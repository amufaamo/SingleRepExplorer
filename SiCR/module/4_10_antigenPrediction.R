# --- UI 関数定義 (変更なし) ---
antigenPredictionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      vdjType(ns),
      downloadButton(ns("download_table"), "テーブルダウンロード (.csv)")
    ),
    mainPanel(
      h4("抗原予測結果テーブル (Antigen Prediction Table - Matched CDR3s)"),
      DTOutput(ns("table"))
    )
  )
}


antigenPredictionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    # --- 抗原データベースの読み込み (変更なし) ---
    antigen_db_tcr_path <- "data/230323_ruft_TCR_antigen_database.tsv"
    antigen_db_bcr_path <- "data/230323_ruft_BCR_antigen_database.tsv"
    print(antigen_db_tcr_path) # デバッグ用
    print(antigen_db_bcr_path) # デバッグ用

    antigen_db <- reactive({
      req(input$vdj_type) # これを追加
      # ファイルが存在するか確認 (オプションですが推奨)
      tcr_exists <- file.exists(antigen_db_tcr_path)
      bcr_exists <- file.exists(antigen_db_bcr_path)

      data <- NULL # 結果を格納する変数
      if (input$vdj_type == "tcr") {
        validate(need(tcr_exists, paste("TCRデータベースファイルが見つかりません:", antigen_db_tcr_path)))
        tryCatch(
          {
            data <- read.delim(antigen_db_tcr_path, header = TRUE, row.names = 1, sep ='\t' )
            print("TCR data read successfully. Head:") # 読み込み成功時に表示
            print(head(data)) # データフレームの先頭を表示
          },
          error = function(e) {
            msg <- paste("TCRファイル読み込みエラー:", e$message)
            print(msg) # コンソールにエラー表示
            validate(need(FALSE, msg)) # UIにエラーメッセージ表示
          }
        )
      } else if (input$vdj_type == "bcr") {
        validate(need(bcr_exists, paste("BCRデータベースファイルが見つかりません:", antigen_db_bcr_path)))
        tryCatch(
          {
            data <- read.delim(antigen_db_bcr_path, header = TRUE, row.names = 1, sep ='\t' )
            print("BCR data read successfully. Head:") # 読み込み成功時に表示
            print(head(data)) # データフレームの先頭を表示
          },
          error = function(e) {
            msg <- paste("BCRファイル読み込みエラー:", e$message)
            print(msg)
            validate(need(FALSE, msg))
          }
        )
      }

      # 読み込んだデータ（またはエラー時はNULL）を返す
      # data が NULL でないか、行数が0より大きいかなどもチェックするとより安全
      validate(need(!is.null(data) && nrow(data) > 0, "データが見つからないか、空です。"))
      return(data)
    })

    merged_df <- reactive({
      req(input$vdj_type)
      req(antigen_db())
      
      if(input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df <- myReactives$tcr_df
        df <- dplyr::inner_join(df, antigen_db(), 
          by = c("TCR_TRB_cdr3" = "CDR3"))
        df <- df %>% dplyr::select(raw_clonotype_id, exact_subclonotype_id, TCR_TRB_cdr3, Epitope_sequence, Epitope_gene, Epitope_species)
      } else if(input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df <- myReactives$bcr_df
        df <- dplyr::inner_join(df, antigen_db(), 
          by = c("BCR_IGH_cdr3" = "CDR3"))
        df <- df %>% dplyr::select(raw_clonotype_id, exact_subclonotype_id, BCR_IGH_cdr3, Epitope_sequence, Epitope_gene, Epitope_species)
      }
      return(df)
    })

    # renderTable の代わりに DT::renderDataTable を使用
    output$table <- DT::renderDataTable(
      {
        merged_df() # reactiveな値を呼び出す
      },
      options = list(pageLength = 10), # DTのオプション（任意）
      rownames = FALSE # 行名を表示しない（任意）
    )

    # --- ダウンロードハンドラ (もし必要なら) ---
    output$download_table <- downloadHandler(
      filename = function() {
        paste0("antigen_prediction_", selected_vdj_type(), "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        readr::write_csv(antigen_db(), file)
      }
    )
  })
}



# load_antigen_db <- function(filepath_reactiveval) {
#     filepath <- filepath_reactiveval(); shiny::validate(shiny::need(file.exists(filepath), ...)); db <- tryCatch({ readr::read_tsv(filepath, show_col_types = FALSE) }, error = function(e){ NULL }); req(db)
#     required_db_cols <- c("CDR3", "Epitope_sequence", "Epitope_gene", "Epitope_species"); shiny::validate(shiny::need(all(required_db_cols %in% colnames(db)), ...))
#     db %>% dplyr::select(all_of(required_db_cols)) %>% dplyr::mutate(CDR3 = as.character(CDR3)) %>% dplyr::filter(!is.na(CDR3) & CDR3 != "") %>% dplyr::distinct(CDR3, .keep_all = TRUE)
# }
# antigen_db <- reactive({ req(input$vdj_type); if(input$vdj_type == "tcr") { load_antigen_db(antigen_db_tcr_path) } else if (input$vdj_type == "bcr") { load_antigen_db(antigen_db_bcr_path) } else { NULL }; shiny::validate(shiny::need(!is.null(db), ...)); return(db) }) # Added validate & return

# # --- リアクティブ: データソースの選択 (変更なし) ---
# reactive_data <- reactive({
#     req(input$vdj_type); df <- NULL; cdr3_col_name <- NULL
#     message(paste("Selecting data for VDJ type:", input$vdj_type))
#     if (input$vdj_type == "tcr") { req(myReactives$tcr_df); df <- myReactives$tcr_df; cdr3_col_name <- "TCR_TRB_cdr3"; shiny::validate(shiny::need(cdr3_col_name %in% colnames(df), ...)) }
#     else if (input$vdj_type == "bcr") { req(myReactives$bcr_df); df <- myReactives$bcr_df; cdr3_col_name <- "BCR_IGH_cdr3"; shiny::validate(shiny::need(cdr3_col_name %in% colnames(df), ...)) }
#     else { return(NULL) }
#     shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, ...))
#     # ★ barcode 列も必要なのでチェックを追加 ★
#     shiny::validate(shiny::need("barcode" %in% colnames(df), "'barcode' column missing from input data."))
#     return(list(data = df, cdr3_col = cdr3_col_name))
# })

# # --- オブザーバー: GroupBy 更新 (変更なし) ---
# observeEvent(reactive_data(), { req(reactive_data()$data); message("Updating group_by choices..."); if (!is.null(myReactives)) update_group_by_select_input(session, myReactives) }, ignoreNULL = TRUE)
# # filter_groups_ui は一旦コメントアウト (必要なら戻す)
# # output$filter_groups_ui <- renderUI({ ... })


# # --- ★★★ リアクティブ: アノテーション付き結合データ (inner_join 結果を直接返す) ★★★ ---
# joined_data <- reactive({
#     message("--- joined_data reactive started ---")
#     data_list <- reactive_data(); req(data_list); antigen_db_data <- antigen_db(); req(antigen_db_data)
#     df <- data_list$data
#     cdr3_col_name <- data_list$cdr3_col
#     group_col <- input$group_by # グループ列も保持
#     # 元のバーコード、その他の必要な列を保持
#     required_user_cols <- unique(c(group_col, cdr3_col_name, "barcode")) # barcode 列を必須に
#     # TCR/BCR データフレームに含まれる他の有用な列も選択候補に加える (任意)
#     # other_cols_to_keep <- c("他の保持したい列1", "他の保持したい列2")
#     # required_user_cols <- unique(c(required_user_cols, other_cols_to_keep))

#     req(df, antigen_db_data, cdr3_col_name, group_col)
#     shiny::validate(shiny::need(all(required_user_cols %in% colnames(df)),
#                     paste("Required columns missing:", paste(setdiff(required_user_cols, colnames(df)), collapse=", "))))

#     message("Preparing user data for join...")
#     convert_to_char <- function(col_vector) {if (is.list(col_vector) && !is.data.frame(col_vector)) {vapply(col_vector, function(x) {if (is.null(x) || length(x) == 0 || all(is.na(x))) { NA_character_ } else if (is.atomic(x)) { paste(as.character(x), collapse = ", ") } else { "[COMPLEX]" }}, FUN.VALUE = character(1), USE.NAMES = FALSE)} else { as.character(col_vector) }}

#     df_prep <- df %>%
#         # ★ 必要な列を選択 (all_of で安全に) ★
#         dplyr::select(all_of(required_user_cols)) %>%
#         dplyr::mutate(
#             cdr3_for_join = convert_to_char(.data[[cdr3_col_name]])
#             # グループ列は型変換不要かも (元の型を保持)
#             # group_col_str = as.character(.data[[group_col]])
#             ) %>%
#          dplyr::filter(!is.na(cdr3_for_join) & cdr3_for_join != "")

#     db_prep <- antigen_db_data # Already prepared

#     message("Joining data with antigen database using INNER JOIN...")
#     df_joined <- tryCatch({
#          # inner_join を使用
#          dplyr::inner_join(df_prep, db_prep, by = c("cdr3_for_join" = "CDR3"), relationship = "many-to-one")
#         }, error = function(e) { warning("Join error: ", e$message); NULL })

#     req(df_joined)
#     message(paste("Rows after INNER join:", nrow(df_joined)))

#     # 結合後の行数が0の場合のバリデーション
#     shiny::validate(shiny::need(nrow(df_joined) > 0, "Antigen Database にマッチする CDR3 が見つかりませんでした。"))

#     # 不要になった結合キー列を削除し、列順序を調整
#     df_final <- df_joined %>%
#         dplyr::select(-cdr3_for_join) %>%
#         # 元のCDR3列は df_prep で選択済みなのでそのまま残る
#         # グループ列、barcode、元のCDR3列、抗原情報の順に並べる
#         dplyr::relocate(all_of(group_col), barcode, all_of(cdr3_col_name), starts_with("Epitope"))

#     message("--- joined_data reactive finished ---")
#     return(df_final)
# })


# # --- リアクティブ: 集計データ (一旦コメントアウト) ---
# # summary_data <- reactive({ ... })


# # --- ★★★ 出力: テーブル (joined_data を表示) ★★★ ---
# output$table <- renderDT({
#     message("--- renderDT started ---")
#     # joined_data() を直接呼び出す
#     joined_df <- joined_data()
#     req(joined_df)
#     shiny::validate(shiny::need(nrow(joined_df) > 0, "No matched data to display in table."))
#     message("Rendering table...")

#     # グループ列名を "Group" に変更（任意）
#     display_df <- joined_df %>%
#          dplyr::rename(Group = all_of(input$group_by)) # all_of を使う

#     datatable(display_df,
#               rownames = FALSE, filter = 'top', extensions="Buttons",
#               options = list(pageLength = 10, scrollX = TRUE, autoWidth = FALSE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
#               caption = "Antigen Prediction Results (Matched Cells)")
# })


# # --- 出力: プロット (無効化) ---
# output$plot <- renderPlot({
#     # プロットは後回し
#     ggplot() + theme_void() + labs(title="Plot generation is currently disabled.")
# })

# # --- ダウンロードハンドラ (テーブルのみ実装) ---
# output$download_plot <- downloadHandler(
#   filename = function() { "plot_disabled.pdf" },
#   content = function(file) { pdf(file); plot.new(); text(0.5,0.5,"Plot disabled"); dev.off() }
# )
# output$download_table <- downloadHandler(
#   filename = function() { paste0("antigen_prediction_matched_data_", input$vdj_type %||% "vdj", "_", input$group_by %||% "group", ".csv") },
#   content = function(file) {
#     # ★ joined_data をダウンロード ★
#     data_to_write <- joined_data(); req(data_to_write)
#     write.csv(data_to_write, file, row.names = FALSE, quote = TRUE, na = "")
#   })

#   }) # moduleServer 終了
# } # antigenPredictionServer 終了
