#source("../utils.R")
source("utils.R")
# --- UI Definition ---
clonalOverlapUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            vdjType(ns),
            selectInput(ns("clone_identifier_column"), "Clonotype Column", choices = NULL),
            groupByInput(ns),
            hr(),
            h4("Analysis Settings"),
            selectInput(ns("overlap_method"), "Overlap Method (Beta Diversity):", 
                        choices = c("Simple Overlap (Intersect/Min)" = "overlap", 
                                    "Jaccard Index" = "jaccard", 
                                    "Morisita-Horn Index" = "morisita")),
            numericInput(ns("top_n_clones"), "Top N Clones for Alluvial:", value = 15, min = 5, max = 100),
            hr(),
            commonPlotOptions(ns) # Common Plot options
        ),
        mainPanel(
            tabsetPanel(
              tabPanel("Overlap Heatmap", 
                 br(),
                 h3("Beta Diversity Heatmap"),
                 downloadButton(ns("download_plot"), "Download plot (.pdf)"),
                 plotOutput(ns("plot"), height = "500px"),
                 hr(),
                 h3("Overlap Matrix"),
                 DTOutput(ns("table"))
              ),
              tabPanel("Alluvial Flow Map",
                 br(),
                 h3("Clonotype Flow across Groups"),
                 downloadButton(ns("download_alluvial"), "Download Map (.pdf)"),
                 plotOutput(ns("alluvial_plot"), height = "600px")
              )
            )
        )
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

      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, paste("選択されたVDJタイプ (", toupper(input$vdj_type %||% ""), ") のデータが見つからないか、空です。")))
      shiny::validate(shiny::need(all(expected_cols %in% colnames(df)), paste("データに必須列 (", paste(expected_cols, collapse=", "), ") が含まれていません。")))
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
    }) |> bindEvent(input$vdj_type, reactive_data, ignoreNULL = FALSE, ignoreInit = FALSE)



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

            # --- Data Preprocessing for counts ---
            df_counts <- tryCatch({
              df %>%
                mutate(.group_str = as.character(.data[[group_col]]),
                       effective_clonotype = as.character(.data[[clone_col]])) %>%
                filter(!is.na(.group_str) & .group_str != "" &
                       !is.na(effective_clonotype) & effective_clonotype != "") %>%
                group_by(.group_str, effective_clonotype) %>%
                summarise(count = n(), .groups = 'drop')
            }, error = function(e) { NULL })

            req(df_counts, cancelOutput = TRUE)

            groups <- sort(unique(df_counts$.group_str))
            shiny::validate(shiny::need(length(groups) >= 2, "Need at least two groups to calculate overlap."))

            # --- Pairwise Calculation Function ---
            calc_overlap_index <- function(g1, g2, data_counts, method) {
                d1 <- data_counts %>% filter(.group_str == g1)
                d2 <- data_counts %>% filter(.group_str == g2)
                
                clones1 <- d1$effective_clonotype
                clones2 <- d2$effective_clonotype
                
                if (method == "overlap") {
                    intersect_size <- length(intersect(clones1, clones2))
                    min_size <- min(length(clones1), length(clones2))
                    return(if (min_size > 0) intersect_size / min_size else 0)
                } 
                else if (method == "jaccard") {
                    intersect_size <- length(intersect(clones1, clones2))
                    union_size <- length(union(clones1, clones2))
                    return(if (union_size > 0) intersect_size / union_size else 0)
                }
                else if (method == "morisita") {
                    common_clones <- intersect(clones1, clones2)
                    if(length(common_clones) == 0) return(0)
                    
                    x_i <- d1 %>% filter(effective_clonotype %in% common_clones) %>% arrange(effective_clonotype) %>% pull(count)
                    y_i <- d2 %>% filter(effective_clonotype %in% common_clones) %>% arrange(effective_clonotype) %>% pull(count)
                    X <- sum(d1$count)
                    Y <- sum(d2$count)
                    sum_x_sq <- sum(d1$count^2)
                    sum_y_sq <- sum(d2$count^2)
                    
                    numerator <- 2 * sum(x_i * y_i)
                    denominator <- (sum_x_sq / X^2 + sum_y_sq / Y^2) * (X * Y)
                    return(if (denominator > 0) numerator / denominator else 0)
                }
                return(0)
            }

            # --- Matrix Loop ---
            mat <- matrix(NA_real_, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
            method_used <- input$overlap_method
            
            tryCatch({
                for (i in seq_along(groups)) {
                    for (j in seq_along(groups)) {
                        g1 <- groups[i]
                        g2 <- groups[j]
                        if (i == j) {
                            mat[g1, g2] <- 1.0
                        } else {
                            mat[g1, g2] <- calc_overlap_index(g1, g2, df_counts, method_used)
                        }
                    }
                }
            }, error = function(e){ NULL })

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
                labs(x = NULL, y = NULL, title = paste("Beta Diversity:", tools::toTitleCase(input$overlap_method))) +

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
        # --- Alluvial Plot ---
        alluvial_obj <- reactive({
            req(reactive_data(), input$group_by, input$clone_identifier_column)
            df <- reactive_data()
            group_col <- input$group_by
            clone_col <- input$clone_identifier_column
            
            shiny::validate(
                shiny::need(group_col %in% colnames(df), paste("Group column missing.")),
                shiny::need(clone_col %in% colnames(df), paste("Clonotype column missing."))
            )
            
            # Count clone occurrences
            df_alluvial <- df %>%
                filter(!is.na(.data[[group_col]]) & !is.na(.data[[clone_col]]) & 
                       .data[[group_col]] != "" & .data[[clone_col]] != "") %>%
                group_by(.data[[group_col]], .data[[clone_col]]) %>%
                summarise(Freq = n(), .groups = 'drop')
                
            # Filter to top N clones based on TOTAL frequency across all groups
            top_clones <- df_alluvial %>%
                group_by(.data[[clone_col]]) %>%
                summarise(TotalFreq = sum(Freq)) %>%
                top_n(input$top_n_clones, TotalFreq) %>%
                pull(.data[[clone_col]])
                
            df_plot <- df_alluvial %>% filter(.data[[clone_col]] %in% top_clones)
            
            shiny::validate(shiny::need(nrow(df_plot) > 0, "No data available for plotting top clones."))
            
            # Use ggalluvial
            p <- ggplot(df_plot,
                   aes(x = .data[[group_col]], stratum = .data[[clone_col]], alluvium = .data[[clone_col]],
                       y = Freq, fill = .data[[clone_col]], label = .data[[clone_col]])) +
                scale_x_discrete(expand = c(.1, .1)) +
                ggalluvial::geom_flow(alpha = 0.5) +
                ggalluvial::geom_stratum(alpha = 0.8) +
                theme_pubr() +
                theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(y = "Frequency", x = tools::toTitleCase(gsub("_", " ", group_col)), 
                     title = paste("Top", input$top_n_clones, "Clonal Transitions"))
            
            return(p)
        })

        output$alluvial_plot <- renderPlot({
            alluvial_obj()
        })
        
        output$download_alluvial <- downloadHandler(
            filename = function() { paste0("alluvial_", input$group_by, ".pdf") },
            content = function(file) {
                 ggsave(file, plot = alluvial_obj(), width = 10, height = 7, device = "pdf")
            }
        )

    }) # moduleServer 終了
}

