# --- UI 関数定義 (修正箇所) ---
trackClonotypeUI <- function(id) {
  ns <- NS(id) # 名前空間を作成
  sidebarLayout(
    sidebarPanel(
      width = 3,
      vdjType(ns),
      selectInput(ns("clone_identifier_column"), "クローン識別子列", choices = NULL),
      selectInput(ns("target_clonotypes"), "追跡対象クローン", choices = NULL, multiple = TRUE, selectize = TRUE),
      groupByInput(ns),

      # ★★★ プロットタイプ選択を追加 ★★★
      selectInput(ns("plot_viz_type"), "プロットタイプ (Plot Type)",
                   choices = c("積み上げ棒グラフ (Stacked Bar)" = "bar",
                               "Alluvial Plot" = "alluvial"),
                   selected = "bar"), # デフォルトは棒グラフ

      radioButtons(ns("value_type"), "表示値", choices = c("頻度 (%)" = "frequency", "細胞数/リード数" = "count"), selected = "frequency"),
      textInput(ns('order'), "グループ順序 (カンマ区切り)", value = NULL),
      selectInput(ns("legend_position"), "凡例位置", choices = c("右" = "right", "左" = "left", "下" = "bottom", "上" = "top", "なし" = "none"), selected = "right"),
      sliderInput(ns("plot_width"), "プロット幅", min = 200, max = 2000, value = 700, step = 50),
      sliderInput(ns("plot_height"), "プロット高さ", min = 200, max = 2000, value = 500, step = 50),
      downloadButton(ns("download_plot"), "プロットダウンロード (.pdf)"),
      downloadButton(ns("download_table"), "テーブルダウンロード (.csv)")
    ),
    mainPanel(
      width = 9,
      h4("クローン追跡プロット"),
      plotOutput(ns("plot")),
      hr(),
      h4("追跡結果テーブル"),
      DTOutput(ns('table'))
    )
  )
}


# --- サーバー関数定義 (修正版) ---
trackClonotypeServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- リアクティブ: データソースの選択 ---
    reactive_data <- reactive({
      req(input$vdj_type)
      df <- NULL
      expected_cols <- c("sample")

      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df <- myReactives$bcr_df
      }

      validate(need(!is.null(df) && nrow(df) > 0, paste("選択されたVDJタイプ (", toupper(input$vdj_type), ") のデータが見つからないか、空です。")))
      validate(need(all(expected_cols %in% colnames(df)), paste("データに必須列 (", paste(expected_cols, collapse=", "), ") が含まれていません。")))
      return(df)
    })

    # --- オブザーバー: UI要素の更新 ---

    # 1. Group By の選択肢更新
    observeEvent(myReactives$seurat_object, { # トリガーは適宜調整
       req(myReactives$seurat_object)
       update_group_by_select_input(session, myReactives)
    }, ignoreNULL = TRUE)

    # 2. Clone Identifier Column の選択肢更新
    observe({
      req(input$vdj_type, nzchar(input$vdj_type))
      df <- reactive_data()
      req(is.data.frame(df) && nrow(df) > 0)

      possible_choices <- list()
      selected_choice_val <- NULL
      all_cols <- colnames(df)

      if (input$vdj_type == "tcr") {
          possible_choices <- c(
              "Raw Clonotype ID" = "raw_clonotype_id", "Exact Clonotype ID" = "exact_subclonotype_id",
              "Pair CDR3 AA" = "TCR_pair_CTaa", "Pair CDR3 NT" = "TCR_pair_CTnt",
              "TRA CDR3 AA" = "TCR_TRA_cdr3", "TRA CDR3 NT" = "TCR_TRA_cdr3_nt",
              "TRB CDR3 AA" = "TCR_TRB_cdr3", "TRB CDR3 NT" = "TCR_TRB_cdr3_nt"
          )
          default_selection_order <- c("TCR_pair_CTaa", "TCR_TRB_cdr3", "raw_clonotype_id")
      } else if (input$vdj_type == "bcr") {
           possible_choices <- c(
               "Raw Clonotype ID" = "raw_clonotype_id", "Exact Clonotype ID" = "exact_subclonotype_id",
               "IGH CDR3 AA" = "BCR_IGH_cdr3", "IGH CDR3 NT" = "BCR_IGH_cdr3_nt",
               "IGK CDR3 AA" = "BCR_IGK_cdr3", "IGK CDR3 NT" = "BCR_IGK_cdr3_nt",
               "IGL CDR3 AA" = "BCR_IGL_cdr3", "IGL CDR3 NT" = "BCR_IGL_cdr3_nt"
           )
           default_selection_order <- c("BCR_IGH_cdr3", "raw_clonotype_id")
      }

      valid_choices_vals <- unname(possible_choices[possible_choices %in% all_cols])
      if(length(valid_choices_vals) == 0){
          final_choices <- c("適切な列が見つかりません" = "")
          selected_choice_val <- ""
      } else {
          final_choices <- possible_choices[possible_choices %in% valid_choices_vals]
          selected_choice_val <- intersect(default_selection_order, valid_choices_vals)[1]
          if(is.na(selected_choice_val) || is.null(selected_choice_val)) selected_choice_val <- valid_choices_vals[1]
      }
       updateSelectInput(session, "clone_identifier_column", choices = final_choices, selected = selected_choice_val)
    }) |> bindEvent(input$vdj_type, reactive_data(), ignoreNULL = FALSE, ignoreInit = FALSE)


    # 3. Target Clonotypes の選択肢更新 (リスト列対応)
    observe({
      req(input$clone_identifier_column, nzchar(input$clone_identifier_column))
      df <- reactive_data()
      req(is.data.frame(df) && nrow(df) > 0)
      clone_col <- input$clone_identifier_column
      validate(need(clone_col %in% colnames(df), paste("列 '", clone_col, "' がデータに存在しません。")))

      available_clonotypes_char <- tryCatch({
          col_data <- df %>% pull(!!sym(clone_col))
          if (is.list(col_data) && !is.data.frame(col_data)) {
              col_data_char <- vapply(col_data, function(x) { # vapply推奨
                  if (is.null(x) || length(x) == 0 || all(is.na(x))) { NA_character_ }
                  else if (is.atomic(x)) { paste(as.character(x), collapse = ", ") }
                  else { "[COMPLEX_LIST_ELEMENT]" }
              }, FUN.VALUE = character(1), USE.NAMES = FALSE)
          } else { col_data_char <- as.character(col_data) }
          sort(unique(na.omit(col_data_char)))
      }, error = function(e){ character(0) })

      validate(need(length(available_clonotypes_char) > 0, paste("列 '", clone_col, "' に有効なクローン識別子が見つかりません。")))
      # Note: Consider server-side selectize if list is huge
      # selected_targets <- head(available_clonotypes_char, 5) # Default selection can be adjusted
      # updateSelectInput(session, "target_clonotypes", choices = available_clonotypes_char, selected = selected_targets)
      # Server-side selectize update (alternative):
       updateSelectizeInput(session, "target_clonotypes",
                            choices = available_clonotypes_char,
                            selected = head(available_clonotypes_char, 5), # Initial selection
                            server = TRUE) # Enable server-side processing

    }) |> bindEvent(input$clone_identifier_column, reactive_data())


    # 4. Group Order のデフォルト値を設定 (リスト列対応)
    observe({
      req(input$group_by, nzchar(input$group_by))
      df <- reactive_data()
      req(is.data.frame(df) && nrow(df) > 0)
      group_col <- input$group_by
      validate(need(group_col %in% colnames(df), paste("グループ列 '", group_col, "' がデータに存在しません。")))

      group_levels_char <- tryCatch({
          col_data <- df %>% pull(!!sym(group_col))
          if (is.list(col_data) && !is.data.frame(col_data)) {
               col_data_char <- vapply(col_data, function(x) { # vapply推奨
                   if (is.null(x) || length(x) == 0 || all(is.na(x))) { NA_character_ }
                   else if (is.atomic(x)) { paste(as.character(x), collapse=", ") }
                   else { "[COMPLEX_LIST_ELEMENT]" }
               }, FUN.VALUE = character(1), USE.NAMES = FALSE)
          } else { col_data_char <- as.character(col_data) }
          sort(unique(na.omit(col_data_char)))
      }, error = function(e){ character(0) })

      if(length(group_levels_char) > 0) {
          updateTextInput(session, "order", value = paste(group_levels_char, collapse = ", "))
      } else {
          updateTextInput(session, "order", value = "")
      }
    }) |> bindEvent(input$group_by, reactive_data())


    # --- リアクティブ: 追跡結果の計算 (Base R table() 使用版) ---
    tracking_results <- reactive({
      req(input$clone_identifier_column, nzchar(input$clone_identifier_column),
          input$group_by, nzchar(input$group_by),
          input$target_clonotypes); validate(need(length(input$target_clonotypes) > 0, "追跡対象のクローンを少なくとも1つ選択してください。"))

      df_orig <- reactive_data(); req(df_orig)
      clone_col_orig <- input$clone_identifier_column; group_col_orig <- input$group_by; targets <- input$target_clonotypes
      validate(need(clone_col_orig %in% colnames(df_orig), paste("列 '", clone_col_orig, "' が存在しません。")))
      validate(need(group_col_orig %in% colnames(df_orig), paste("列 '", group_col_orig, "' が存在しません。")))

      showNotification("クローン追跡データを計算中...", id="calc_tracking", duration = NULL, type = "message"); on.exit(removeNotification(id="calc_tracking"), add = TRUE)

      # --- 安全な文字列表現を作成 ---
      df <- tryCatch({
          convert_to_char <- function(col_vector) {
              if (is.list(col_vector) && !is.data.frame(col_vector)) {
                  vapply(col_vector, function(x) {
                      if (is.null(x) || length(x) == 0 || all(is.na(x))) { NA_character_ }
                      else if (is.atomic(x)) { paste(as.character(x), collapse = ", ") }
                      else { "[COMPLEX_LIST_ELEMENT]" }
                  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
              } else { as.character(col_vector) }
          }
          df_processed <- df_orig %>%
              mutate(
                  .group_chr = convert_to_char(.data[[group_col_orig]]),
                  .clone_chr = convert_to_char(.data[[clone_col_orig]])
              ) %>%
              filter(!is.na(.group_chr) & !is.na(.clone_chr))
          # stopifnot 確認は省略 (必要なら復活)
          df_processed
          }, error = function(e){showNotification(paste("データの前処理中にエラー:", e$message), type="error"); req(FALSE)})
      validate(need(nrow(df) > 0, "処理可能なデータが見つかりません（NA除去後）。"))

      # --- 1. Total counts (Base R table) ---
      message("--- Preparing for total_counts calculation (using table) ---")
      total_counts <- tryCatch({
          counts_table <- table(df$.group_chr)
          data.frame(.group_chr = names(counts_table), total_n = as.integer(counts_table), stringsAsFactors = FALSE) %>% filter(total_n > 0)
          }, error = function(e){ msg <- paste("Total counts の table() でエラー発生:", e$message); showNotification(msg, type="error"); warning(msg); req(FALSE)})
      message("--- total_counts calculated successfully (using table) ---")
      valid_groups_chr <- total_counts %>% pull(.group_chr); validate(need(length(valid_groups_chr) > 0, "有効なグループが見つかりません。"))

      # --- 2. Target counts (Base R table + 安全な列名割り当て) ---
      message("--- Preparing for target_counts calculation (using table) ---")
      target_counts <- tryCatch({
          df_filtered_for_target <- df %>% filter(.group_chr %in% valid_groups_chr & .clone_chr %in% targets)
          if (nrow(df_filtered_for_target) == 0) {
              message("Target counts: No rows after filtering for targets.")
              data.frame(.group_chr = character(0), .clone_chr = character(0), n = integer(0), stringsAsFactors = FALSE)
          } else {
              counts_table_2d <- table(df_filtered_for_target$.group_chr, df_filtered_for_target$.clone_chr)
              target_df_long <- as.data.frame(counts_table_2d, stringsAsFactors = FALSE)
              message("Target counts: Structure after as.data.frame(table):"); try(print(str(target_df_long)), silent=TRUE); try(print(head(target_df_long)), silent=TRUE)
              message("Actual column names:"); print(colnames(target_df_long))
              if (nrow(target_df_long) > 0 && ncol(target_df_long) == 3) {
                   message("Target counts: Found 3 columns. Assigning names '.group_chr', '.clone_chr', 'n' by position.")
                   colnames(target_df_long) <- c(".group_chr", ".clone_chr", "n")
                   target_df_long %>% filter(n > 0)
              } else {
                   warning("Target counts: Result from as.data.frame(table) did not have 3 columns. Ncol = ", ncol(target_df_long), ". Returning empty result.")
                   data.frame(.group_chr = character(0), .clone_chr = character(0), n = integer(0), stringsAsFactors = FALSE)
              }
          }
          }, error = function(e){ msg <- paste("Target counts の table()/as.data.frame()/colnames() 処理中にエラー発生:", e$message); showNotification(msg, type="error"); warning(msg); req(FALSE)})
      message("--- target_counts calculated successfully (using table) ---")

      # --- 3 & 4 (expand_grid, left_join, mutate) ---
      full_grid <- expand_grid(.group_chr = valid_groups_chr, .clone_chr = targets)
      results <- full_grid %>%
        left_join(target_counts, by = c(".group_chr", ".clone_chr")) %>%
        left_join(total_counts, by = ".group_chr") %>%
        mutate(n = replace_na(n, 0), total_n = replace_na(total_n, 0)) %>%
        mutate(frequency = ifelse(total_n > 0, (n / total_n) * 100, 0)) %>%
        mutate(.clone_fct = factor(.clone_chr, levels = targets),
               .group_fct = factor(.group_chr, levels = unique(str_sort(valid_groups_chr, numeric=TRUE))))
      return(results)
    }) # reactive 終了

    # --- ★★★ 出力: 追跡プロット (修正: プロットタイプ分岐) ★★★ ---
    output$plot <- renderPlot({
      results_df <- tracking_results()
      req(results_df)
      value_col <- input$value_type
      plot_type <- input$plot_viz_type # ★ プロットタイプを取得

      y_col_sym <- if(value_col == "frequency") sym("frequency") else sym("n")
      y_axis_label <- if(value_col == "frequency") "頻度 (%)" else "細胞数 / リード数"

      order_input <- input$order %||% ""
      specified_order_chr <- trimws(strsplit(order_input, ", *")[[1]])
      available_groups_chr <- levels(results_df$.group_fct) # Factorのレベルを使用
      plot_order_chr <- specified_order_chr[specified_order_chr %in% available_groups_chr]
      if(length(plot_order_chr) == 0) plot_order_chr <- available_groups_chr

      # --- プロットタイプに応じて ggplot の中身を分岐 ---
      if (plot_type == "bar") {
          # --- 積み上げ棒グラフ (従来のコード) ---
          p <- ggplot(results_df, aes(x = .group_fct, y = !!y_col_sym, fill = .clone_fct)) +
            geom_col(position = "stack", width=0.8) +
            scale_x_discrete(limits = plot_order_chr, name = toTitleCase(gsub("_", " ", input$group_by))) +
            ylab(y_axis_label) +
            ggtitle(paste("追跡対象クローンの", y_axis_label, "推移 (Stacked Bar)")) +
            labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))) +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
                  legend.position = input$legend_position, plot.title = element_text(hjust = 0.5))

      } else if (plot_type == "alluvial") {
          # --- Alluvial Plot (ggalluvial を使用) ---
           # Alluvial Plot では Factor 型でない方が扱いやすい場合があるので確認
           # results_df$.group_chr <- as.character(results_df$.group_fct)
           # results_df$.clone_chr <- as.character(results_df$.clone_fct)
           # ↑ Factorのままでも ggalluvial は通常動作します

          # 追跡対象クローンが多すぎると見づらくなる警告 (任意)
          if(length(unique(results_df$.clone_fct)) > 15) {
              showNotification("追跡対象クローンが多数選択されています。Alluvial Plot が見づらくなる可能性があります。", type="warning", duration=5)
          }

          p <- ggplot(results_df,
                      aes(x = .group_fct,       # X軸: グループ (Factor型を推奨)
                          y = !!y_col_sym,      # Y軸: 値 (頻度 or カウント)
                          alluvium = .clone_fct,# 流れの要素: クローン
                          stratum = .clone_fct, # 各X軸での層: クローン
                          fill = .clone_fct)) + # 色分け: クローン
                 # 層間の流れを描画 (widthで太さ調整, alphaで透明度)
                 ggalluvial::geom_flow(width = 0.4, alpha = 0.6, na.rm = TRUE) +
                 # 各X軸位置の層（バーに相当）を描画
                 ggalluvial::geom_stratum(width = 0.4, na.rm = TRUE) +
                 # 層の中にラベルを表示（任意、クローン数が多いと重なる）
                 # geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5, na.rm = TRUE) +
                 scale_x_discrete(limits = plot_order_chr, name = toTitleCase(gsub("_", " ", input$group_by))) +
                 ylab(y_axis_label) +
                 ggtitle(paste("追跡対象クローンの", y_axis_label, "推移 (Alluvial Plot)")) +
                 labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))) +
                 theme_bw(base_size = 12) +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
                       legend.position = input$legend_position,
                       plot.title = element_text(hjust = 0.5))
      } else {
           # 未知のプロットタイプの場合はエラー表示またはデフォルト表示
           p <- ggplot() + annotate("text", x=0, y=0, label="無効なプロットタイプです") + theme_void()
      }
      # --- 分岐終了 ---

      print(p)

    }, width = reactive(input$plot_width), height = reactive(input$plot_height))


    # --- ★★★ ダウンロードハンドラ (修正: プロットタイプ分岐) ★★★ ---
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("clonotype_tracking_", input$vdj_type, "_by_", input$group_by, "_", input$value_type, "_", input$plot_viz_type, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        # renderPlot と同様のロジックでプロットオブジェクトを再生成
        results_df <- tracking_results(); req(results_df)
        value_col <- input$value_type; y_col_sym <- if(value_col == "frequency") sym("frequency") else sym("n"); y_axis_label <- if(value_col == "frequency") "頻度 (%)" else "細胞数 / リード数"
        order_input <- input$order %||% ""; specified_order_chr <- trimws(strsplit(order_input, ", *")[[1]]); available_groups_chr <- levels(results_df$.group_fct); plot_order_chr <- specified_order_chr[specified_order_chr %in% available_groups_chr]; if(length(plot_order_chr) == 0) plot_order_chr <- available_groups_chr
        plot_type <- input$plot_viz_type # ★ プロットタイプ取得

        if (plot_type == "bar") {
           p <- ggplot(results_df, aes(x = .group_fct, y = !!y_col_sym, fill = .clone_fct)) + geom_col(position = "stack", width=0.8) # 省略... (renderPlotと同じコード)
                # ... (scale_x_discrete, ylab, ggtitle, labs, theme_bw, theme を追加)
        } else if (plot_type == "alluvial") {
            p <- ggplot(results_df, aes(x = .group_fct, y = !!y_col_sym, alluvium = .clone_fct, stratum = .clone_fct, fill = .clone_fct)) +
                 ggalluvial::geom_flow(width = 0.4, alpha = 0.6, na.rm = TRUE) + ggalluvial::geom_stratum(width = 0.4, na.rm = TRUE) # 省略... (renderPlotと同じコード)
                 # ... (scale_x_discrete, ylab, ggtitle, labs, theme_bw, theme を追加)
        } else {
            p <- ggplot() + annotate("text", x=0, y=0, label="Invalid Plot Type") + theme_void()
        }
         # --- renderPlot内の残りのコードをここにコピー ---
         # (共通部分)
         p <- p + scale_x_discrete(limits = plot_order_chr, name = toTitleCase(gsub("_", " ", input$group_by))) +
              ylab(y_axis_label) +
              ggtitle(paste("追跡対象クローンの", y_axis_label, "推移 (", plot_type, ")")) +
              labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))) +
              theme_bw(base_size = 12) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
                    legend.position = input$legend_position,
                    plot.title = element_text(hjust = 0.5))
        # --- ここまで ---


        plot_width_in <- input$plot_width / 96; plot_height_in <- input$plot_height / 96
        ggsave(file, plot = p, device = "pdf", width = max(plot_width_in, 4), height = max(plot_height_in, 4), units = "in", dpi = 300) # サイズ下限を少し調整
      })


    # --- ダウンロードハンドラ ---
    output$download_plot <- downloadHandler(
      filename = function() { paste0("clonotype_tracking_", input$vdj_type, "_by_", input$group_by, "_", input$value_type, "_", Sys.Date(), ".pdf") },
      content = function(file) {
        results_df <- tracking_results(); req(results_df)
        value_col <- input$value_type; y_col_sym <- if(value_col == "frequency") sym("frequency") else sym("n"); y_axis_label <- if(value_col == "frequency") "頻度 (%)" else "細胞数 / リード数"
        order_input <- input$order %||% ""; specified_order_chr <- trimws(strsplit(order_input, ", *")[[1]]); available_groups_chr <- levels(results_df$.group_fct); plot_order_chr <- specified_order_chr[specified_order_chr %in% available_groups_chr]; if(length(plot_order_chr) == 0) plot_order_chr <- available_groups_chr
        p <- ggplot(results_df, aes(x = .group_fct, y = !!y_col_sym, fill = .clone_fct)) +
             geom_col(position = "stack", width=0.8) + scale_x_discrete(limits = plot_order_chr, name = toTitleCase(gsub("_", " ", input$group_by))) + ylab(y_axis_label) + ggtitle(paste("追跡対象クローンの", y_axis_label, "推移")) + labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10), legend.position = input$legend_position, plot.title = element_text(hjust = 0.5))
        plot_width_in <- input$plot_width / 96; plot_height_in <- input$plot_height / 96
        ggsave(file, plot = p, device = "pdf", width = max(plot_width_in, 3), height = max(plot_height_in, 3), units = "in", dpi = 300)
      })

    output$download_table <- downloadHandler(
      filename = function() { paste0("clonotype_tracking_table_", input$vdj_type, "_by_", input$group_by, "_", Sys.Date(), ".csv") },
      content = function(file) {
        results_df <- tracking_results(); req(results_df)
        download_df <- results_df %>% rename(`Group` = .group_chr, `Clonotype` = .clone_chr, `Count` = n, `Total Count in Group` = total_n, `Frequency (%)` = frequency) %>% mutate(`Frequency (%)` = round(`Frequency (%)`, 5)) %>% select(Group, Clonotype, Count, `Frequency (%)`, `Total Count in Group`)
        write.csv(download_df, file, row.names = FALSE, quote = TRUE, na = "")
      })

  }) # moduleServer 終了
} # trackClonotypeServer 終了

# --- プレースホルダー/既存の関数 (diversityAnalysis から流用・確認) ---
# vdjType, groupByInput, commonPlotOptions, update_group_by_select_input は
# diversityAnalysis モジュールで使用されているものと同じものが利用可能と仮定。
# もし未定義の場合は、前の回答で示したプレースホルダー等を参考に定義してください。