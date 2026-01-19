#source("../utils.R")
source("utils.R")
# --- UI 関数定義 (デバッグ用・最小版) ---
clonalOverlapUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            vdjType(ns),
            selectInput(ns("clone_identifier_column"), "クローン識別子列", choices = NULL),
            groupByInput(ns),
            commonPlotOptions(ns) # 共通プロットオプション
        ),
        mainPanel(
            h3("Plot"),
            downloadButton(ns("download_plot"), "Download plot (.pdf)"),
            plotOutput(ns("plot")),
            h3("Table"),
            DTOutput(ns("table")) # テーブルID変更
        ),
    )
}

# --- サーバー関数定義 (デバッグ用・最小版) ---
clonalOverlapServer <- function(id, myReactives) {
    moduleServer(id, function(input, output, session) {
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

    observeEvent(myReactives$seurat_object, {
      update_group_by_select_input(session, myReactives)
    })

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

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



        # --- ★ リアクティブ: 重複度計算 (修正: shiny::validate, shiny::need) ★ ---
        overlap_matrix_debug <- reactive({
            # データ取得と固定設定
            req(reactive_data())
            df <- reactive_data()
            group_col <- input$group_by
            clone_col <- input$clone_identifier_column
            overlap_method <- "overlap"

            # ★★★ shiny::validate と shiny::need を使用 ★★★
            shiny::validate(
                shiny::need(group_col %in% colnames(df), paste("Group column '", group_col, "' not found in tcr_df.")),
                shiny::need(clone_col %in% colnames(df), paste("Clonotype column '", clone_col, "' not found in tcr_df."))
            )
            # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★

            message("--- Starting overlap calculation (Debug Mode) ---")
            showNotification("Calculating overlap (Debug Mode)...", duration = NULL, id = "overlap_calc_debug")
            on.exit(removeNotification("overlap_calc_debug"), add = TRUE)

            # --- データ前処理 ---
            df_processed <- tryCatch(
                {
                    message("Preprocessing data...")
                    convert_to_char <- function(col_vector) {
                        if (is.list(col_vector) && !is.data.frame(col_vector)) {
                            vapply(col_vector, function(x) {
                                if (is.null(x) || length(x) == 0 || all(is.na(x))) {
                                    NA_character_
                                } else if (is.atomic(x)) {
                                    paste(as.character(x), collapse = ", ")
                                } else {
                                    "[COMPLEX]"
                                }
                            }, FUN.VALUE = character(1), USE.NAMES = FALSE)
                        } else {
                            as.character(col_vector)
                        }
                    }
                    temp_df <- df %>%
                        mutate(
                            .group_str = as.character(.data[[group_col]]),
                            effective_clonotype = as.character(.data[[clone_col]])
                        ) %>%
                        filter(!is.na(.group_str) & .group_str != "" &
                            !is.na(effective_clonotype) & effective_clonotype != "") %>%
                        dplyr::select(.group_str, effective_clonotype) %>%
                        distinct()
                    message("Preprocessing finished.")
                    if (!all(c(".group_str", "effective_clonotype") %in% colnames(temp_df))) {
                        stop("Preprocessing failed to create required columns.")
                    }
                    temp_df
                },
                error = function(e) {
                    err_msg <- paste("Data processing error:", e$message)
                    warning(err_msg)
                    showNotification(err_msg, type = "error", duration = 10)
                    NULL
                }
            )

            # ★★★ shiny::validate と shiny::need を使用 ★★★
            req(df_processed, cancelOutput = TRUE)
            shiny::validate(
                shiny::need(is.data.frame(df_processed), "Preprocessing result ('df_processed') is not a data frame."),
                shiny::need(nrow(df_processed) > 0, "No valid group/clonotype pairs found after filtering (df_processed is empty)."),
                shiny::need(".group_str" %in% colnames(df_processed), "Required column '.group_str' is missing after preprocessing."),
                shiny::need("effective_clonotype" %in% colnames(df_processed), "Required column 'effective_clonotype' is missing after preprocessing.")
            )
            # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
            message("--- df_processed successfully created ---")
            message("Structure:")
            try(print(str(df_processed)), silent = TRUE)
            message("Head:")
            try(print(head(df_processed)), silent = TRUE)
            message("-------------------------------------------")

            # --- グループリスト取得 ---
            groups <- tryCatch(
                {
                    unique(df_processed$.group_str)
                },
                error = function(e) {
                    warning("Could not get unique groups. Error: ", e$message)
                    NULL
                }
            )
            req(groups, cancelOutput = TRUE)
            # ★★★ shiny::validate と shiny::need を使用 ★★★
            shiny::validate(shiny::need(length(groups) >= 2, "Need at least two groups to calculate overlap."))
            # ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
            groups <- sort(groups)
            message(paste("Groups found:", paste(groups, collapse = ", ")))

            # --- ペアワイズ計算関数 ---
            calc_overlap_simple <- function(group1, group2, data_grouped) {
                clones1 <- data_grouped %>%
                    filter(.group_str == group1) %>%
                    pull(effective_clonotype)
                clones2 <- data_grouped %>%
                    filter(.group_str == group2) %>%
                    pull(effective_clonotype)
                if (length(clones1) == 0 || length(clones2) == 0) {
                    return(0)
                }
                intersect_size <- length(intersect(clones1, clones2))
                min_size <- min(length(clones1), length(clones2))
                value <- if (min_size > 0) intersect_size / min_size else 0
                return(value)
            }

            # --- 行列計算ループ ---
            mat <- matrix(NA_real_, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
            message("Calculating matrix...")
            tryCatch(
                {
                    for (i in 1:length(groups)) {
                        for (j in 1:length(groups)) {
                            g1 <- groups[i]
                            g2 <- groups[j]
                            if (i == j) {
                                mat[g1, g2] <- 1.0
                            } else {
                                mat[g1, g2] <- calc_overlap_simple(g1, g2, df_processed)
                            }
                        }
                    }
                    message("Matrix calculation complete.")
                },
                error = function(e) {
                    err_msg <- paste("Error during matrix calculation loop:", e$message)
                    warning(err_msg)
                    showNotification(err_msg, type = "error", duration = 10)
                    return(matrix(NA_real_, nrow = 0, ncol = 0))
                }
            )
            req(mat, cancelOutput = TRUE)

            return(mat) # 行列のみを返す
        }) # reactive 終了

        # --- プロットオブジェクトの生成 (reactive) ---
        plot_obj <- reactive({
            mat <- overlap_matrix_debug()
            req(mat)
            shiny::validate(shiny::need(is.matrix(mat) && all(dim(mat) > 0), "Overlap matrix is empty or invalid for plotting."))

            df_plot_long <- tryCatch(
                {
                    as.data.frame(as.table(mat), stringsAsFactors = FALSE)
                },
                error = function(e) {
                    showNotification(paste("Error converting matrix for plot:", e$message), type = "error")
                    NULL
                }
            )
            req(df_plot_long, cancelOutput = TRUE)

            if (ncol(df_plot_long) == 3) {
                colnames(df_plot_long) <- c("Group1", "Group2", "Value")
            } else {
                shiny::validate(shiny::need(FALSE, "Cannot create plot data from matrix conversion result."))
            }

            df_plot <- df_plot_long %>% mutate(Group1 = factor(Group1, levels = rownames(mat)), Group2 = factor(Group2, levels = colnames(mat)))
            fill_limits <- c(0, 1)
            p <- ggplot(df_plot, aes(x = Group1, y = Group2, fill = Value)) +
                geom_tile(color = "grey50") +
                scale_fill_viridis_c(option = "plasma", limits = fill_limits, na.value = "grey80", name = "Overlap") +
                geom_text(aes(label = round(Value, 2)), color = "white", size = 3.5, na.rm = TRUE, check_overlap = TRUE) +
                coord_fixed() +
                labs(x = NULL, y = NULL, title = "Clonal Overlap Heatmap") +
                theme_minimal(base_size = 11) +
                theme(
                    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                    panel.grid = element_blank(),
                    legend.position = input$legend %||% "right"
                )
            return(p)
        })

        # --- プロット出力 ---
        output$plot <- renderPlot(
            { plot_obj() },
            width = reactive(input$plot_width %||% 600),
            height = reactive(input$plot_height %||% 500)
        )

        # --- テーブル出力 ---
        output$table <- renderDT({
            mat <- overlap_matrix_debug()
            req(mat)
            shiny::validate(shiny::need(is.matrix(mat) || is.data.frame(mat), "Overlap result is not a matrix/dataframe for table."))
            display_df <- as.data.frame(mat) %>%
                rownames_to_column("Group") %>%
                mutate(across(where(is.numeric), ~ round(., 3)))
            datatable(display_df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = min(10, nrow(mat))))
        })

        # --- PDFダウンロード機能 ---
        output$download_plot <- downloadHandler(
            filename = function() {
                paste0("clonal_overlap_", input$group_by, "_", input$clone_identifier_column, ".pdf")
            },
            content = function(file) {
                p <- plot_obj()
                req(p)
                ggsave(file, plot = p, width = (input$plot_width %||% 600) / 72, height = (input$plot_height %||% 500) / 72, device = "pdf", dpi = 300)
            }
        )
    }) # moduleServer 終了
} # clonalOverlapServer_debug 終了
