# --- UI 関数定義 ---
diversityAnalysisUI <- function(id) {
  ns <- NS(id) # 名前空間を作成
  sidebarLayout(
    sidebarPanel(
      # --- UI要素を配置 ---
      vdjType(ns),
      
      # ★★★ 新しいUI要素: 多様性指数の選択 ★★★
      selectInput(ns("diversity_index"), "多様性指数 (Diversity Index)",
                  choices = c("Shannon" = "shannon",
                              "Simpson" = "simpson",
                              "Inverse Simpson" = "invsimpson",
                              "Chao1 (推定種数)" = "chao1"),
                  selected = "shannon"),
      
      selectInput(ns("target_gene_column"), "クローンID列 (Target Column)",
                  choices = NULL, # Server側で設定
                  selected = NULL),
      groupByInput(ns),
      commonPlotOptions(ns),
    ),
    mainPanel(
      tabsetPanel(
        # ★★★ タブのタイトルを汎用的に変更 ★★★
        tabPanel("Diversity Plot",
                 plotOutput(ns("shannonPlot"), height = "500px")
        ),
        tabPanel("Diversity Table",
                 DTOutput(ns("resultsTable"))
        )
      )
    )
  )
}
# --- サーバー関数定義 ---
diversityAnalysisServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- リアクティブ: データソースの選択 ---
    # (この部分は変更ありません)
    reactive_data <- reactive({
      req(input$vdj_type)
      df <- NULL
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df <- myReactives$bcr_df
      }
      validate(need(!is.null(df) && nrow(df) > 0, paste("選択されたVDJタイプ (", toupper(input$vdj_type), ") のデータが見つからないか、空です。")))
      validate(need("sample" %in% colnames(df), "データに 'sample' 列が必要です。"))
      return(df)
    }) |> bindCache(input$vdj_type, myReactives$tcr_df, myReactives$bcr_df)

    # --- UIの動的更新 ---
    # (この部分は変更ありません)
    observe({
      req(input$vdj_type, nzchar(input$vdj_type))
      df <- reactive_data()
      req(is.data.frame(df) && nrow(df) > 0)
      
      possible_choices <- list()
      selected_choice_val <- NULL
      all_cols <- colnames(df)
      
      if (input$vdj_type == "tcr") {
        possible_choices <- c(
          "Raw Clonotype ID" = "raw_clonotype_id",
          "Exact Clonotype ID" = "exact_subclonotype_id",
          "Pair CDR3 AA" = "TCR_pair_CTaa",
          "Pair CDR3 NT" = "TCR_pair_CTnt",
          "TRA CDR3 AA" = "TCR_TRA_cdr3",
          "TRA CDR3 NT" = "TCR_TRA_cdr3_nt",
          "TRB CDR3 AA" = "TCR_TRB_cdr3",
          "TRB CDR3 NT" = "TCR_TRB_cdr3_nt"
        )
        default_selection_order <- c("raw_clonotype_id", "TCR_pair_CTaa", "TCR_TRB_cdr3")
      } else if (input$vdj_type == "bcr") {
        possible_choices <- c(
          "Raw Clonotype ID" = "raw_clonotype_id",
          "Exact Clonotype ID" = "exact_subclonotype_id",
          "IGH CDR3 AA" = "BCR_IGH_cdr3_aa",
          "IGH CDR3 NT" = "BCR_IGH_cdr3_nt",
          "IGK CDR3 AA" = "BCR_IGK_cdr3_aa",
          "IGK CDR3 NT" = "BCR_IGK_cdr3_nt",
          "IGL CDR3 AA" = "BCR_IGL_cdr3_aa",
          "IGL CDR3 NT" = "BCR_IGL_cdr3_nt"
        )
        default_selection_order <- c("raw_clonotype_id", "BCR_IGH_cdr3_aa")
      } else {
        possible_choices <- c("Select VDJ Type First" = "")
        default_selection_order <- c("")
      }
      
      valid_choices_vals <- unname(possible_choices[possible_choices %in% all_cols])
      
      if(length(valid_choices_vals) == 0){
        final_choices <- c("適切なクローンID列が見つかりません" = "")
        selected_choice_val <- ""
      } else {
        final_choices <- possible_choices[possible_choices %in% valid_choices_vals]
        selected_choice_val <- intersect(default_selection_order, valid_choices_vals)[1]
        if(is.na(selected_choice_val) || is.null(selected_choice_val)) {
          selected_choice_val <- valid_choices_vals[1]
        }
      }
      
      updateSelectInput(session, "target_gene_column",
                        choices = final_choices,
                        selected = selected_choice_val)
    }) |> bindEvent(input$vdj_type, reactive_data(), ignoreNULL = FALSE, ignoreInit = FALSE)

    
    # --- リアクティブ: 解析結果の計算 (★★★ バグを修正したコード ★★★) ---
    analysis_results <- reactive({
      # 必要な入力が存在することを確認
      req(input$target_gene_column, nzchar(input$target_gene_column),
          input$group_by, nzchar(input$group_by),
          input$diversity_index)
      df <- reactive_data()
      req(df)

      clone_id_col <- input$target_gene_column
      group_col <- input$group_by
      index_method <- input$diversity_index # 例: "shannon"

      # --- ★★★ 修正点1: 正しい名前のルックアップ方法 ★★★ ---
      # UIのchoicesで定義した名前と値の対応関係をここで定義します
      index_lookup <- c(
        "Shannon" = "shannon",
        "Simpson" = "simpson",
        "Inverse Simpson" = "invsimpson",
        "Chao1 (推定種数)" = "chao1"
      )
      # サーバーに渡ってきた値（例: "shannon"）から、UI上の表示名（例: "Shannon"）を検索します
      index_display_name <- names(index_lookup)[which(index_lookup == index_method)]
      
      # 安全対策：もし名前が見つからなければ処理を停止
      validate(need(length(index_display_name) > 0, "多様性指数の表示名が見つかりません。"))

      # 正しく取得した表示名で通知を表示
      showNotification(paste(index_display_name, "Index 計算中..."), id="calc_diversity", duration=NULL, type="message")

      # --- 列の存在チェック ---
      validate(need(clone_id_col %in% colnames(df), paste("クローンID列 '", clone_id_col, "' がデータに存在しません。")))
      validate(need(group_col %in% colnames(df), paste("グループ列 '", group_col, "' がデータに存在しません。")))
      validate(need("sample" %in% colnames(df), "'sample' 列がデータに含まれていません。"))

      # --- 計算準備 (abundance_matrix の作成) ---
      sample_metadata <- df %>%
        dplyr::select(sample, !!sym(group_col)) %>%
        mutate(sample = as.character(sample)) %>%
        distinct(sample, .keep_all = TRUE)
      
      count_matrix_long <- df %>%
        mutate(sample = as.character(sample)) %>%
        filter(!is.na(.data[[clone_id_col]])) %>%
        group_by(sample, !!sym(clone_id_col)) %>%
        summarise(count = n(), .groups = 'drop')
      
      count_matrix_wide <- count_matrix_long %>%
        pivot_wider(names_from = !!sym(clone_id_col), values_from = count, values_fill = 0)
      
      abundance_matrix <- count_matrix_wide %>% dplyr::select(-sample) %>% as.matrix()
      rownames(abundance_matrix) <- count_matrix_wide$sample

      # --- 多様性指数の計算 ---
      index_scores <- tryCatch({
        if (index_method %in% c("shannon", "simpson", "invsimpson")) {
          vegan::diversity(abundance_matrix, index = index_method)
        } else if (index_method == "chao1") {
          chao_results <- vegan::estimateR(t(abundance_matrix))
          chao_results[2, ]
        } else {
          stop("無効な多様性指数です。")
        }
      }, error = function(e){
        showNotification(paste(index_display_name, "Index の計算中にエラー:", e$message), type = "error", duration=10)
        return(NULL)
      })

      if(is.null(index_scores)) {
        removeNotification(id="calc_diversity")
        req(FALSE, cancelOutput = TRUE)
      }

      # --- ★★★ 修正点2: 正しい表示名を列名として使用 ★★★ ---
      # 結果をデータフレームに整形します
      results_df <- data.frame(sample = as.character(rownames(abundance_matrix)))
      # ここで正しく取得した表示名を列名として設定します (エラー発生箇所でした)
      results_df[[index_display_name]] <- index_scores
      
      results_df <- results_df %>%
        left_join(sample_metadata, by = "sample")

      removeNotification(id="calc_diversity")
      
      # 結果をリストで返します
      return(list(results = results_df, index_name = index_display_name))
    })

    # --- 出力: 結果テーブル ---
    # (この部分は変更ありません)
    output$resultsTable <- renderDT({
      res_list <- analysis_results()
      req(res_list, res_list$results)

      group_col_name <- input$group_by %||% "Group"
      index_col_name <- res_list$index_name

      display_df <- res_list$results %>%
        rename(!!tools::toTitleCase(gsub("_", " ", group_col_name)) := !!sym(group_col_name))

      datatable(display_df, extensions = 'Buttons',
                options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE,
                               dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                rownames = FALSE,
                caption = paste("各サンプルの", index_col_name, "Index とグループ情報"))
    })

    # --- 出力: 多様性プロット ---
    # (この部分は変更ありません)
    output$shannonPlot <- renderPlot({
      res_list <- analysis_results()
      req(res_list, res_list$results)

      group_col <- input$group_by
      index_col_name <- res_list$index_name
      validate(need(group_col %in% names(res_list$results), paste("グループ列 '", group_col, "' が結果に存在しません。")))
      validate(need(index_col_name %in% names(res_list$results), paste("指数列 '", index_col_name, "' が結果に存在しません。")))

      plot_data <- res_list$results %>%
        filter(!is.na(.data[[index_col_name]]) & !is.na(.data[[group_col]]))
      validate(need(nrow(plot_data) > 0, paste("表示する有効な", index_col_name, "score がありません。")))

      plot_data <- plot_data %>%
        arrange(.data[[group_col]], sample) %>%
        mutate(sample = factor(sample, levels = unique(sample)))

      ggplot(plot_data, aes(x = sample, y = .data[[index_col_name]], fill = .data[[group_col]])) +
        geom_bar(stat = "identity", width = 0.8) +
        labs(
          y = paste("Index Score (", index_col_name, ")"),
          x = "Sample",
          title = paste(index_col_name, "Index per Sample"),
          fill = tools::toTitleCase(gsub("_", " ", group_col))
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "right"
        )
    }, res = 96)

  }) # moduleServer 終了
}