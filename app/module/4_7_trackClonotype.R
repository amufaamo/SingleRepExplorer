# Clonotype tracking module
# --- UI 関数定義 (修正箇所) ---
trackClonotypeUI <- function(id) {
  ns <- NS(id) # Create namespace
  sidebarLayout(
    sidebarPanel(
      width = 3,
      vdjType(ns),
      selectInput(ns("clone_identifier_column"), "Clone Identifier Column", choices = NULL),
      selectInput(ns("target_clonotypes"), "Target Clonotypes to Track", choices = NULL, multiple = TRUE, selectize = TRUE),
      groupByInput(ns),

      selectInput(ns("plot_viz_type"), "Plot Type",
                   choices = c("Stacked Bar Plot" = "bar",
                               "Alluvial Plot" = "alluvial"),
                   selected = "bar"),

      radioButtons(ns("value_type"), "Display Value", choices = c("Frequency (%)" = "frequency", "Cell/Read Count" = "count"), selected = "frequency"),
      textInput(ns('order'), "Group Order (comma-separated)", value = NULL),
      selectInput(ns("legend_position"), "Legend Position", choices = c("Right" = "right", "Left" = "left", "Bottom" = "bottom", "Top" = "top", "None" = "none"), selected = "right"),
      div(style = "display: flex; gap: 10px;",
          numericInput(ns("plot_width"), "Plot Width (px)", min = 100, max = 2000, value = 750, step = 100),
          numericInput(ns("plot_height"), "Plot Height (px)", min = 100, max = 2000, value = 500, step = 100)
      )
    ),
    mainPanel(
      width = 9,
      h3("Clonotype Tracking Plot"),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns("plot")),
      hr(),
      h4("Tracking Results Table"),
      downloadButton(ns("download_table"), "Download Table (.xlsx)"),
      br(), br(),
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

      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, paste("No data found for the selected VDJ type (", toupper(input$vdj_type %||% ""), ").")))
      shiny::validate(shiny::need(all(expected_cols %in% colnames(df)), paste("Mandatory columns (", paste(expected_cols, collapse=", "), ") missing from data.")))
      return(df)
    })

    # --- オブザーバー: UI要素の更新 ---

    # 1. Group By の選択肢更新
    observeEvent(myReactives$seurat_object, { # トリガーは適宜調整
       req(myReactives$seurat_object)
       update_group_by_select_input(session, myReactives)
    }, ignoreNULL = TRUE)

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    # 2. Clone Identifier Column の選択肢更新（タイミング問題修正）
    observe({
      req(input$vdj_type, nzchar(input$vdj_type))

      vdj_df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(!is.null(vdj_df) && nrow(vdj_df) > 0)

      possible_choices <- list()
      selected_choice_val <- NULL
      all_cols <- colnames(vdj_df)

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
               "IGH CDR3 AA" = "BCR_IGH_cdr3_aa", "IGH CDR3 NT" = "BCR_IGH_cdr3_nt",
               "IGK CDR3 AA" = "BCR_IGK_cdr3_aa", "IGK CDR3 NT" = "BCR_IGK_cdr3_nt",
               "IGL CDR3 AA" = "BCR_IGL_cdr3_aa", "IGL CDR3 NT" = "BCR_IGL_cdr3_nt"
           )
           default_selection_order <- c("BCR_IGH_cdr3_aa", "raw_clonotype_id")
      }

      valid_choices_vals <- unname(possible_choices[possible_choices %in% all_cols])
      if(length(valid_choices_vals) == 0){
          final_choices <- c("No suitable columns found" = "")
          selected_choice_val <- ""
      } else {
          final_choices <- possible_choices[possible_choices %in% valid_choices_vals]
          selected_choice_val <- intersect(default_selection_order, valid_choices_vals)[1]
          if(is.na(selected_choice_val) || is.null(selected_choice_val)) selected_choice_val <- valid_choices_vals[1]
      }
       updateSelectInput(session, "clone_identifier_column", choices = final_choices, selected = selected_choice_val)
    }) |> bindEvent(input$vdj_type, myReactives$tcr_df, myReactives$bcr_df, ignoreNULL = FALSE, ignoreInit = FALSE)


    # 3. Target Clonotypes の選択肢更新 (リスト列対応)
    observe({
      req(input$clone_identifier_column, nzchar(input$clone_identifier_column))
      df <- reactive_data()
      req(is.data.frame(df) && nrow(df) > 0)
      clone_col <- input$clone_identifier_column
      validate(need(clone_col %in% colnames(df), paste("Column '", clone_col, "' does not exist in data.")))

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

      validate(need(length(available_clonotypes_char) > 0, paste("No valid clone identifiers found in column '", clone_col, "'.")))
      # Note: Consider server-side selectize if list is huge
      # selected_targets <- head(available_clonotypes_char, 5) # Default selection can be adjusted
      # updateSelectInput(session, "target_clonotypes", choices = available_clonotypes_char, selected = selected_targets)
      # Server-side selectize update (alternative):
       updateSelectizeInput(session, "target_clonotypes",
                            choices = available_clonotypes_char,
                            selected = head(available_clonotypes_char, 5), # Initial selection
                            server = TRUE) # Enable server-side processing

    }) |> bindEvent(input$clone_identifier_column, reactive_data)


    # 4. Group Order のデフォルト値を設定 (リスト列対応)
    observe({
      req(input$group_by, nzchar(input$group_by))
      df <- reactive_data()
      req(is.data.frame(df) && nrow(df) > 0)
      group_col <- input$group_by
      validate(need(group_col %in% colnames(df), paste("Group column '", group_col, "' does not exist in data.")))

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
    }) |> bindEvent(input$group_by, reactive_data)


    # --- リアクティブ: 追跡結果の計算 (Base R table() 使用版) ---
    tracking_results <- reactive({
      req(input$clone_identifier_column, nzchar(input$clone_identifier_column),
          input$group_by, nzchar(input$group_by),
          input$target_clonotypes); validate(need(length(input$target_clonotypes) > 0, "Please select at least one clone to track."))

      df_orig <- reactive_data(); req(df_orig)
      clone_col_orig <- input$clone_identifier_column; group_col_orig <- input$group_by; targets <- input$target_clonotypes
      shiny::validate(shiny::need(group_col_orig %in% colnames(df_orig), paste("Column '", group_col_orig, "' does not exist.")))

      showNotification("Computing clonotype tracking data...", id="calc_tracking", duration = NULL, type = "message"); on.exit(removeNotification(id="calc_tracking"), add = TRUE)

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
          }, error = function(e){showNotification(paste("Error during data preprocessing:", e$message), type="error"); req(FALSE)})
      validate(need(nrow(df) > 0, "No processable data found after NA removal."))

      # --- 1. Total counts (Base R table) ---
      message("--- Preparing for total_counts calculation (using table) ---")
      total_counts <- tryCatch({
          counts_table <- table(df$.group_chr)
          data.frame(.group_chr = names(counts_table), total_n = as.integer(counts_table), stringsAsFactors = FALSE) %>% filter(total_n > 0)
          }, error = function(e){ msg <- paste("Error while computing total counts:", e$message); showNotification(msg, type="error"); warning(msg); req(FALSE)})
      message("--- total_counts calculated successfully (using table) ---")
      valid_groups_chr <- total_counts %>% pull(.group_chr); validate(need(length(valid_groups_chr) > 0, "No valid groups found."))

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
          }, error = function(e){ msg <- paste("Error while computing target counts:", e$message); showNotification(msg, type="error"); warning(msg); req(FALSE)})
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
      y_axis_label <- if(value_col == "frequency") "Frequency (%)" else "Cell/Read Count"

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
            ggtitle(paste("Tracking of target clones (Stacked Bar)")) +
            labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))) +
            theme_classic(base_size = 14) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
                  legend.position = input$legend_position, plot.title = element_text(hjust = 0.5))

      } else if (plot_type == "alluvial") {
          # --- Alluvial Plot (ggalluvial を使用) ---
           # Alluvial Plot では Factor 型でない方が扱いやすい場合があるので確認
           # results_df$.group_chr <- as.character(results_df$.group_fct)
           # results_df$.clone_chr <- as.character(results_df$.clone_fct)
           # ↑ Factorのままでも ggalluvial は通常動作します

          # Notice if too many clones
          if(length(unique(results_df$.clone_fct)) > 15) {
              showNotification("Many clones selected. Alluvial Plot may become cluttered.", type="warning", duration=5)
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
                 ggtitle(paste("Tracking of target clones (Alluvial Plot)")) +
                 labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))) +
                 theme_classic(base_size = 14) +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
                       legend.position = input$legend_position,
                       plot.title = element_text(hjust = 0.5))
      } else {
           # Invalid plot type
           p <- ggplot() + annotate("text", x=0, y=0, label="Invalid Plot Type") + theme_void()
      }
      # --- 分岐終了 ---

      print(p)

    }, width = reactive(input$plot_width), height = reactive(input$plot_height))


    # --- ダウンロードハンドラ (PPTXに変更) ---
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("clonotype_tracking_", input$vdj_type, "_by_", input$group_by,
               "_", input$value_type, "_", input$plot_viz_type, "_", Sys.Date(), ".pptx")
      },
      content = function(file) {
        results_df <- tracking_results()
        req(results_df)
        value_col   <- input$value_type
        y_col_sym   <- if (value_col == "frequency") sym("frequency") else sym("n")
        y_axis_label <- if (value_col == "frequency") "Frequency (%)" else "Cell/Read Count"
        plot_type   <- input$plot_viz_type

        order_input       <- input$order %||% ""
        specified_order   <- trimws(strsplit(order_input, ", *")[[1]])
        available_groups  <- levels(results_df$.group_fct)
        plot_order        <- specified_order[specified_order %in% available_groups]
        if (length(plot_order) == 0) plot_order <- available_groups

        common_theme <- list(
          scale_x_discrete(limits = plot_order,
                           name = toTitleCase(gsub("_", " ", input$group_by))),
          ylab(y_axis_label),
          ggtitle(paste("Clonotype tracking (", plot_type, ")")),
          labs(fill = toTitleCase(gsub("_", " ", input$clone_identifier_column))),
          theme_classic(base_size = 14),
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
                legend.position = input$legend_position,
                plot.title = element_text(hjust = 0.5))
        )

        p <- if (plot_type == "bar") {
          ggplot(results_df, aes(x = .group_fct, y = !!y_col_sym, fill = .clone_fct)) +
            geom_col(position = "stack", width = 0.8)
        } else if (plot_type == "alluvial") {
          ggplot(results_df,
                 aes(x = .group_fct, y = !!y_col_sym,
                     alluvium = .clone_fct, stratum = .clone_fct, fill = .clone_fct)) +
            ggalluvial::geom_flow(width = 0.4, alpha = 0.6, na.rm = TRUE) +
            ggalluvial::geom_stratum(width = 0.4, na.rm = TRUE)
        } else {
          ggplot() + annotate("text", x = 0, y = 0, label = "Invalid Plot Type") + theme_void()
        }
        p <- p + common_theme

        save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
      })

    output$download_table <- downloadHandler(
      filename = function() { paste0("clonotype_tracking_table_", input$vdj_type, "_by_", input$group_by, "_", Sys.Date(), ".xlsx") },
      content = function(file) {
        results_df <- tracking_results(); req(results_df)
        download_df <- results_df %>% rename(`Group` = .group_chr, `Clonotype` = .clone_chr, `Count` = n, `Total Count in Group` = total_n, `Frequency (%)` = frequency) %>% mutate(`Frequency (%)` = round(`Frequency (%)`, 5)) %>% select(Group, Clonotype, Count, `Frequency (%)`, `Total Count in Group`)
        openxlsx::write.xlsx(download_df, file)
      })

  }) # moduleServer 終了
} # trackClonotypeServer 終了

# --- プレースホルダー/既存の関数 (diversityAnalysis から流用・確認) ---
# vdjType, groupByInput, commonPlotOptions, update_group_by_select_input は
# diversityAnalysis モジュールで使用されているものと同じものが利用可能と仮定。
# もし未定義の場合は、前の回答で示したプレースホルダー等を参考に定義してください。
