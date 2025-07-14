# --- UI Definition ---
cdrLengthDistUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      vdjType(ns),
      selectInput(ns("target_cdr_column"), "Target Column (for Length Calculation)",
        choices = c("Select VDJ type first..." = ""),
        selected = ""
      ),
      groupByInput(ns),
      uiOutput(ns("filter_groups_ui")), # 重複していたものを1つに修正

      selectInput(ns("plot_type"), "Plot Type",
        choices = c(
          "Dodged Bars" = "dodge",
          "Stacked Bars" = "stack",
          "Boxplot" = "boxplot"
        ),
        selected = "dodge"
      ),
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] != 'boxplot'"),
        selectInput(ns("display_value"), "Display Value",
          choices = c("Count" = "count", "Percentage (%)" = "percentage"),
          selected = "count"
        )
      ),
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] != 'boxplot'"),
        checkboxInput(ns("horizontal"), "Horizontal Barplot", value = FALSE)
      ),
      hr(),
      commonPlotOptions(ns)
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download Plot (.pdf)"),
      plotOutput(ns('plot')),
      hr(),
      h3("Table"),
      downloadButton(ns("download_table"), "Download Table (.csv)"),
      DTOutput(ns("table")),
    )
  )
}

# --- Server Definition ---
cdrLengthDistServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    observeEvent(myReactives$seurat_object, {
      update_group_by_select_input(session, myReactives)
    })

    observeEvent(input$vdj_type, {
      req(input$vdj_type) # VDJタイプが選択されていることを確認
      
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        target_choices_all <- c(
          "TCR Paired CTnt" = "TCR_pair_CTnt", "TCR Paired CTaa" = "TCR_pair_CTaa",
          "TCR TRA CDR3 aa" = "TCR_TRA_cdr3", "TCR TRA CDR3 nt" = "TCR_TRA_cdr3_nt",
          "TCR TRB CDR3 aa" = "TCR_TRB_cdr3", "TCR TRB CDR3 nt" = "TCR_TRB_cdr3_nt"
        )
        selected_choice <- "TCR_TRB_cdr3"
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        target_choices_all <- c(
          "BCR Paired CTnt" = "BCR_pair_CTnt", "BCR Paired CTaa" = "BCR_pair_CTaa",
          "BCR IGH CDR3 aa" = "BCR_IGH_cdr3", "BCR IGH CDR3 nt" = "BCR_IGH_cdr3_nt",
          "BCR IGK CDR3 aa" = "BCR_IGK_cdr3", "BCR IGK CDR3 nt" = "BCR_IGK_cdr3_nt",
          "BCR IGL CDR3 aa" = "BCR_IGL_cdr3", "BCR IGL CDR3 nt" = "BCR_IGL_cdr3_nt"
        )
        selected_choice <- "BCR_IGH_cdr3"
      } else {
        target_choices_all <- c("Select VDJ type first..." = "")
        selected_choice <- ""
      }
      updateSelectInput(session, "target_cdr_column", choices = target_choices_all, selected = selected_choice)
    })

    reactive_df_raw <- reactive({
      df <- NULL # 未定義エラーを防ぐためにNULLで初期化
      req(input$vdj_type)
      
      if (input$vdj_type == "tcr" && !is.null(myReactives$tcr_df)) {
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr" && !is.null(myReactives$bcr_df)) {
        df <- myReactives$bcr_df
      }
      return(df)
    })

    table <- reactive({
      # ★★★★★ 修正点 ★★★★★
      # 計算に必要なデータと入力が揃っていることを確認
      req(reactive_df_raw(), input$target_cdr_column)
      
      target_col_sym <- sym(input$target_cdr_column)

      df <- reactive_df_raw() %>%
        # 対象列にNAや空文字が含まれているとncharでエラーになることがあるため、除外する
        filter(!is.na(!!target_col_sym) & !!target_col_sym != "") %>%
        mutate(
          length = nchar(as.character(!!target_col_sym))
        )
      
      df <- df %>% dplyr::select(input$group_by, length)
      df <- df %>%
        group_by(.data[[input$group_by]], length) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(.data[[input$group_by]]) %>%
        mutate(percentage = count / sum(count) * 100) %>% # パーセンテージ計算を修正
        ungroup()
      return(df)
    })

    output$filter_groups_ui <- renderUI({
      df <- reactive_df_raw()
      req(df, input$group_by)
      grouping_var <- input$group_by

      available_groups <- df[[grouping_var]] %>% unique() %>% sort()

      validate(
        need(length(available_groups) > 0, paste("No available groups in column:", shQuote(grouping_var)))
      )

      checkboxGroupInput(session$ns("filter_groups"),
        label = paste("Filter", tools::toTitleCase(gsub("_", " ", grouping_var)), ":"),
        choices = available_groups,
        selected = available_groups,
        inline = TRUE
      )
    })

    table_filtered <- reactive({
      df <- table()
      req(df, input$group_by, input$filter_groups)

      df_filtered <- df %>%
        dplyr::filter(.data[[input$group_by]] %in% input$filter_groups)

      validate(
        need(nrow(df_filtered) > 0, "No data remaining after filtering. Please adjust filter.")
      )
      
      return(df_filtered)
    })

    boxplot_df <- reactive({
      # ★★★★★ 修正点 ★★★★★
      req(reactive_df_raw(), input$target_cdr_column)
      
      target_col_sym <- sym(input$target_cdr_column)

      df <- reactive_df_raw() %>%
        filter(!is.na(!!target_col_sym) & !!target_col_sym != "") %>%
        mutate(
          length = nchar(as.character(!!target_col_sym))
        )
      df <- df %>% dplyr::select(input$group_by, length)
      
      req(df, input$group_by, input$filter_groups)

      df_filtered <- df %>%
        dplyr::filter(.data[[input$group_by]] %in% input$filter_groups)

      validate(
        need(nrow(df_filtered) > 0, "No data remaining after filtering. Please adjust filter.")
      )
      
      return(df_filtered)
    })

    output$plot <- renderPlot({
      req(table_filtered())
      
      p <- if (input$plot_type == "boxplot") {
        req(boxplot_df())
        boxplot_df() %>% ggplot(aes(x = .data[[input$group_by]], y = length, fill = .data[[input$group_by]])) +
          geom_boxplot() +
          labs(x = input$group_by, y = "Length", fill = input$group_by)
      } else { # dodge or stack
        y_val <- if(input$display_value == "count") "count" else "percentage"
        p <- table_filtered() %>% ggplot(aes(x = length, y = .data[[y_val]], fill = .data[[input$group_by]])) +
          geom_bar(stat = 'identity', position = input$plot_type) +
          labs(x = "Length", y = tools::toTitleCase(y_val), fill = input$group_by) +
          scale_y_continuous(expand = c(0, 0))

        if (input$horizontal) {
          p <- p + coord_flip()
        }
        p
      }
      
      p + theme_classic()
    })

    output$table <- renderDT({
      table_filtered()
    })
  })
}