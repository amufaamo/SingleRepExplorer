# --- UI 関数定義 ---
diversityAnalysisUI <- function(id) {
  ns <- NS(id) # 名前空間を作成
  sidebarLayout(
    sidebarPanel(
      # --- UI要素を配置 ---
      vdjType(ns),
      selectInput(ns("target_gene_column"), "クローンID列 (Target Column)",
                  choices = NULL, # Server側で設定
                  selected = NULL),
      groupByInput(ns),
      commonPlotOptions(ns),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Shannon Index Plot",
                 # プロット出力時のエラーメッセージ表示用に uiOutput を使うことも検討可
                 plotOutput(ns("shannonPlot"), height = "500px")
        ),
        tabPanel("Shannon Index Table",
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

    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
      # TODO: target_gene_column の選択肢更新 observeEvent
    })


    # # vdjType 変更時に Target Gene Column の選択肢を更新
    # observeEvent(input$vdj_type,
    #   {
    #     req(input$vdj_type, input$vdj_type != "")

    #     target_choices <- list()
    #     selected_choice <- NULL
    #     # ★ 新しいV-Jペアの列名を定義 ★
    #     vj_pair_names <- list(
    #       tcr_tra = "TCR_TRA_v_j_pair",
    #       tcr_trb = "TCR_TRB_v_j_pair",
    #       bcr_igh = "BCR_IGH_v_j_pair",
    #       bcr_igl = "BCR_IGL_v_j_pair"
    #     )

    #     if (input$vdj_type == "tcr") {
    #       target_choices <- c(
    #         # 単一遺伝子
    #         "Clonotype ID" = 'raw_clonotype_id',
    #         "TRA V Gene" = "TCR_TRA_v_gene",
    #         "TRA J Gene" = "TCR_TRA_j_gene",
    #         "TRB V Gene" = "TCR_TRB_v_gene",
    #         "TRB D Gene" = "TCR_TRB_d_gene",
    #         "TRB J Gene" = "TCR_TRB_j_gene",
    #         # V-J ペア
    #         "TRA V-J Pair" = vj_pair_names$tcr_tra, # 定義した列名を使用
    #         "TRB V-J Pair" = vj_pair_names$tcr_trb # 定義した列名を使用
    #       )
    #       selected_choice <- "TCR_TRB_v_gene"
    #     } else if (input$vdj_type == "bcr") {
    #       target_choices <- c(
    #         # 単一遺伝子
    #         "IGH V Gene" = "BCR_IGH_v_gene",
    #         "IGH D Gene" = "BCR_IGH_d_gene",
    #         "IGH J Gene" = "BCR_IGH_j_gene",
    #         "IGH C Gene" = "BCR_IGH_c_gene",
    #         "IGK V Gene" = "BCR_IGK_v_gene",
    #         "IGK J Gene" = "BCR_IGK_j_gene",
    #         "IGL V Gene" = "BCR_IGL_v_gene",
    #         "IGL J Gene" = "BCR_IGL_j_gene",
    #         # V-J ペア
    #         "IGH V-J Pair" = vj_pair_names$bcr_igh,
    #         "IGK V-J Pair" = vj_pair_names$bcr_igk, # IGK追加
    #         "IGL V-J Pair" = vj_pair_names$bcr_igl
    #       )
    #       selected_choice <- "BCR_IGH_v_gene"
    #     } else {
    #       target_choices <- list("Select VDJ Type First" = "")
    #       selected_choice <- ""
    #     }

    #     # updateSelectInput で選択肢を更新
    #     updateSelectInput(session, "target_gene_column",
    #       choices = target_choices,
    #       selected = selected_choice
    #     )
    #   },
    #   ignoreNULL = FALSE,
    #   ignoreInit = FALSE
    # )

        observe({
        # input$vdj_type が NULL または空文字列でないことを確認
        req(input$vdj_type, nzchar(input$vdj_type))
        df <- reactive_data() # データソースを取得
        # データフレームが NULL でなく、行があることを確認
        req(is.data.frame(df) && nrow(df) > 0)

        possible_choices <- list() # 候補リストを初期化
        selected_choice_val <- NULL # デフォルト選択値を初期化
        all_cols <- colnames(df) # 実際のデータに含まれる列名を取得

        # VDJタイプに応じて候補リストを定義
        if (input$vdj_type == "tcr") {
            # TCR用の候補リスト (表示名 = 実際のカラム名)
            # ユーザー指定のリストに基づいて定義
            # 注: 実際のデータカラム名に合わせて調整済み (例: exact_subclonotype_id, TCR_TRA_cdr3 など)
            possible_choices <- c(
                "Raw Clonotype ID" = "raw_clonotype_id",
                "Exact Clonotype ID" = "exact_subclonotype_id",
                "Pair CDR3 AA" = "TCR_pair_CTaa",
                "Pair CDR3 NT" = "TCR_pair_CTnt",
                "TRA CDR3 AA" = "TCR_TRA_cdr3", # ユーザー指定 TRA CDR3 (TCR_TRA_CTaa) に対応
                "TRA CDR3 NT" = "TCR_TRA_cdr3_nt",# ユーザー指定 TRA CDR3 (TCR_TRA_CTnt) に対応
                "TRB CDR3 AA" = "TCR_TRB_cdr3", # ユーザー指定 TRB CDR3 (TCR_TRB_CTaa) に対応
                "TRB CDR3 NT" = "TCR_TRB_cdr3_nt" # ユーザー指定 TRB CDR3 (TCR_TRB_CTnt) に対応
            )
            # デフォルト選択の優先順位
            default_selection_order <- c("raw_clonotype_id", "TCR_pair_CTaa", "TCR_TRB_cdr3")

        } else if (input$vdj_type == "bcr") {
            # BCR用の候補リスト (表示名 = 想定されるカラム名)
            # ユーザー指定のリストに基づいて定義
            # 注: Pair CDR3 はデータにない可能性が高い。IGH, IGL/IGK の CDR3 カラム名を想定
            possible_choices <- c(
                "Raw Clonotype ID" = "raw_clonotype_id",
                "Exact Clonotype ID" = "exact_subclonotype_id",
                # "Pair CDR3 AA" = "BCR_pair_CTaa", # データにない可能性が高いのでコメントアウト
                # "Pair CDR3 NT" = "BCR_pair_CTnt", # データにない可能性が高いのでコメントアウト
                "IGH CDR3 AA" = "BCR_IGH_cdr3_aa", # ユーザー指定 IGH CDR3 (BCR_IGH_CTaa) に対応する想定
                "IGH CDR3 NT" = "BCR_IGH_cdr3_nt", # ユーザー指定 IGH CDR3 (BCR_IGH_CTnt) に対応する想定
                "IGK CDR3 AA" = "BCR_IGK_cdr3_aa", # IGKも候補に含める
                "IGK CDR3 NT" = "BCR_IGK_cdr3_nt", # IGKも候補に含める
                "IGL CDR3 AA" = "BCR_IGL_cdr3_aa", # ユーザー指定 IGL CDR3 (BCR_IGL_CTaa) に対応する想定
                "IGL CDR3 NT" = "BCR_IGL_cdr3_nt"  # ユーザー指定 IGL CDR3 (BCR_IGL_CTnt) に対応する想定
            )
             # デフォルト選択の優先順位
            default_selection_order <- c("raw_clonotype_id", "BCR_IGH_cdr3_aa")

        } else {
            # VDJタイプが無効な場合
            possible_choices <- c("Select VDJ Type First" = "")
            default_selection_order <- c("")
        }

        # 実際にデータフレームに存在するカラムのみを最終的な選択肢とする
        # possible_choices の値 (カラム名) が all_cols に含まれるものをフィルタリング
        valid_choices_vals <- unname(possible_choices[possible_choices %in% all_cols])

        # 存在する選択肢がない場合の処理
        if(length(valid_choices_vals) == 0){
            final_choices <- c("適切なクローンID列が見つかりません" = "")
            selected_choice_val <- ""
        } else {
            # 表示名と値のペアを保持したまま、存在するカラムのみにフィルタリング
            final_choices <- possible_choices[possible_choices %in% valid_choices_vals]
            # デフォルト選択: 優先順位リストに基づき、最初に存在するものを選ぶ
            selected_choice_val <- intersect(default_selection_order, valid_choices_vals)[1]
            # もし優先順位リストで見つからなければ、存在する選択肢の最初のものを選ぶ
            if(is.na(selected_choice_val) || is.null(selected_choice_val)) {
                selected_choice_val <- valid_choices_vals[1]
            }
        }

        # selectInput を更新
        updateSelectInput(session, "target_gene_column",
                          choices = final_choices,
                          selected = selected_choice_val)

    # VDJタイプまたはデータソースが変更された時にこの observe を実行
    # ignoreNULL=FALSE: vdj_type が NULL になっても実行 (初期状態など)
    # ignoreInit=FALSE: アプリ起動時/モジュール初期化時にも実行
    }) |> bindEvent(input$vdj_type, reactive_data(), ignoreNULL = FALSE, ignoreInit = FALSE)


    # --- リアクティブ: 解析結果の計算 (★ reactive に変更 ★) ---
    analysis_results <- reactive({
      # 必要な入力が存在し、空でないことを確認
      req(input$target_gene_column, nzchar(input$target_gene_column),
          input$group_by, nzchar(input$group_by))
      df <- reactive_data() # データ取得 (エラーは reactive_data 内で処理)
      req(df) # データがない場合はここで停止

      clone_id_col <- input$target_gene_column # UIで選択されたクローンID列
      group_col <- input$group_by           # UIで選択されたグループ列

      # 列の存在を最終確認
      validate(need(clone_id_col %in% colnames(df), paste("クローンID列 '", clone_id_col, "' がデータに存在しません。")))
      validate(need(group_col %in% colnames(df), paste("グループ列 '", group_col, "' がデータに存在しません。")))
      validate(need("sample" %in% colnames(df), "'sample' 列がデータに含まれていません。"))

      # --- 計算処理 (withProgress は reactive 内では直接使えないことが多いので省略するか、工夫が必要) ---
      # showNotificationなどで処理開始を通知する方が良いかも
      showNotification("Shannon Diversity 計算中...", id="calc_shannon", duration=NULL, type="message")

      # サンプルごとのメタデータ
      sample_metadata <- df %>%
        dplyr::select(sample, !!sym(group_col)) %>%
        mutate(sample = as.character(sample)) %>%
        distinct(sample, .keep_all = TRUE)

      # クローン頻度集計
      count_matrix_long <- df %>%
        mutate(sample = as.character(sample)) %>%
        filter(!is.na(.data[[clone_id_col]])) %>%
        group_by(sample, !!sym(clone_id_col)) %>%
        summarise(count = n(), .groups = 'drop')

      # ワイド形式変換
      count_matrix_wide <- count_matrix_long %>%
        pivot_wider(names_from = !!sym(clone_id_col), values_from = count, values_fill = 0)

      # 数値行列化
      abundance_matrix <- count_matrix_wide %>% dplyr::select(-sample) %>% as.matrix()
      rownames(abundance_matrix) <- count_matrix_wide$sample

      # Shannon Diversity 計算
      shannon_index <- tryCatch({
          diversity(abundance_matrix, index = "shannon")
      }, error = function(e){
          showNotification(paste("Shannon Index の計算中にエラー:", e$message), type = "error", duration=10)
          return(NULL) # エラー時はNULLを返す
      })
      # 計算に失敗したら通知を消して停止
      if(is.null(shannon_index)) {
          removeNotification(id="calc_shannon")
          req(FALSE, cancelOutput = TRUE) # 以降の処理をキャンセル
      }

      # 結果整形と結合
      results_df <- data.frame(sample = as.character(rownames(abundance_matrix)), Shannon = shannon_index) %>%
          left_join(sample_metadata, by = "sample")

      removeNotification(id="calc_shannon") # 完了通知

      # 結果を返す
      return(list(results = results_df))

    }) # reactive 終了


    # --- 出力: 結果テーブル ---
    output$resultsTable <- renderDT({
      # analysis_results() を呼び出すと、必要な入力が揃っていれば自動的に計算が実行される
      res_list <- analysis_results()
      # 結果が NULL でないことを確認 (計算エラーの場合など)
      req(res_list, res_list$results)

      group_col_name <- input$group_by %||% "Group"
      display_df <- res_list$results %>%
          rename(!!tools::toTitleCase(gsub("_", " ", group_col_name)) := !!sym(group_col_name))

      datatable(display_df, extensions = 'Buttons',
                options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE,
                               dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                rownames = FALSE, caption = "各サンプルの Shannon Diversity Index とグループ情報")
    })

    # --- 出力: Shannon Diversity ドットプロット ---
    output$shannonPlot <- renderPlot({
      res_list <- analysis_results()
      # 解析結果があり、結果データフレームが存在することを確認
      req(res_list, res_list$results)
      # グループ列名を取得
      group_col <- input$group_by
      # グループ列が結果に含まれているか確認
      validate(need(group_col %in% names(res_list$results), paste("グループ列 '", group_col, "' が結果に存在しません。")))

      # プロット用データを準備 (NAを除外)
      plot_data <- res_list$results %>%
          filter(!is.na(Shannon) & !is.na(.data[[group_col]]))
      # プロットするデータが存在するか確認
      validate(need(nrow(plot_data) > 0, "表示する有効な Shannon score がありません。"))

      # X軸の順序設定: グループでまとめて、その中でサンプル名順にバーが並ぶようにする
      # グループとサンプル名でソートし、その順序で sample を Factor 型に変換する
      plot_data <- plot_data %>%
          arrange(.data[[group_col]], sample) %>%
          mutate(sample = factor(sample, levels = unique(sample)))

      # ggplotでバープロットを作成
      ggplot(plot_data, aes(x = sample, y = Shannon, fill = .data[[group_col]])) +
        # geom_bar(stat = "identity") で、yの値がバーの高さになるようにする
        geom_bar(stat = "identity", width = 0.8) + # バーの幅を少し調整
        # ラベル設定
        labs(
          y = "Index Score (Shannon)", # Y軸ラベル
          x = "Sample",                # X軸ラベル
          title = "Shannon Diversity Index per Sample", # プロットタイトル (任意)
          fill = tools::toTitleCase(gsub("_", " ", group_col)) # 凡例のタイトルを整形
        ) +
        # テーマ設定
        theme_minimal(base_size = 14) + # シンプルなテーマ
        theme(
          # X軸のテキスト（サンプル名）が重ならないように90度回転
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), # 文字サイズを少し小さく
          axis.ticks.x = element_blank(),       # X軸の目盛り線は不要
          panel.grid.major.x = element_blank(), # X軸方向のグリッド線も非表示
          panel.grid.minor.x = element_blank(),
          legend.position = "right"             # 凡例を右側に表示
        )
    }, res = 96) # 解像度

    # --- ダウンロードハンドラは省略 ---

  }) # moduleServer 終了
}