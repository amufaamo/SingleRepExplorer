# クローンサイズ割合のUIとサーバーモジュール

# --- UI 関数定義 ---
clonalProportionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      vdjType(ns),
      groupByInput(ns),
      radioButtons(ns('clone_call'), label = "クローン定義方法 (Clone Definition)",
                   choices = list("Clonotype ID" = "raw_clonotype_id",
                                  "CDR3 nucleotide" = "nt",
                                  "CDR3 amino acid" = "aa"),
                   selected = "raw_clonotype_id"),
      uiOutput(ns("chain_ui_wrapper")),

      h4("Clonal Size Categories (Upper Limits):"),
      numericInput(ns("split_1"), label = "1. Unique Max Count (e.g., 1)", value = 1, min = 1),
      numericInput(ns("split_2"), label = "2. Rare Max Count (e.g., 5)", value = 5, min = 1),
      numericInput(ns("split_3"), label = "3. Small Max Count (e.g., 20)", value = 20, min = 1),
      numericInput(ns("split_4"), label = "4. Medium Max Count (e.g., 100)", value = 100, min = 1),
      numericInput(ns("split_5"), label = "5. Large Max Count (e.g., 1000)", value = 1000, min = 1),
      numericInput(ns("split_6"), label = "6. Hyperexpanded Max Count (e.g., 10000)", value = 10000, min = 1),

      selectInput(ns("legend_position"), "凡例位置", choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
      sliderInput(ns("plot_width"), "プロット幅", min = 200, max = 2000, value = 700, step = 50),
      sliderInput(ns("plot_height"), "プロット高さ", min = 200, max = 2000, value = 500, step = 50),

      downloadButton(ns("download_plot"), "プロットダウンロード (.pdf)"),
      downloadButton(ns("download_table"), "テーブルダウンロード (.csv)")
    ),
    mainPanel(
      width = 9,
      h4("グループ別クローンサイズ割合 (Clonal Proportion by Size Category)"),
      plotOutput(ns("plot")),
      hr(),
      h4("割合データテーブル (Proportion Data Table)"),
      DTOutput(ns('table'))
    )
  )
}

# --- サーバー関数定義 ---
clonalProportionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    debug_mode <- TRUE

    reactive_data <- reactive({
      req(input$vdj_type)
      vdj_type <- input$vdj_type
      message("Using vdj_type from input: ", vdj_type)
      df <- NULL
      if (vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df <- myReactives$tcr_df
      } else if (vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df <- myReactives$bcr_df
      } else {
        stop("Invalid vdj_type: ", vdj_type)
      }
      validate(shiny::need(!is.null(df) && nrow(df) > 0, "Data is NULL or empty."))
      validate(shiny::need("sample" %in% colnames(df), "'sample' column not found."))
      return(df)
    })

    observeEvent(reactive_data(), { 
       req(reactive_data())
    }, ignoreNULL = TRUE)

    output$chain_ui_wrapper <- renderUI({
      choices <- list("Both" = "both", "TRA" = "TRA", "TRB" = "TRB")
      radioButtons(ns('chain'), label = "対象チェーン (Chain)", choices = choices, selected = 'both', inline = TRUE)
    })

    proportion_data <- reactive({
      message("--- proportion_data calculation started ---")
      req(input$chain, nzchar(input$chain))
      req(is.numeric(input$split_1), is.numeric(input$split_2), is.numeric(input$split_3),
          is.numeric(input$split_4), is.numeric(input$split_5), is.numeric(input$split_6),
          reactive_data())
      df_orig <- reactive_data(); req(df_orig)

      group_col <- input$group_by
      clone_call_method <- input$clone_call
      chain_selection <- input$chain

      split_values <- sort(c(input$split_1, input$split_2, input$split_3, input$split_4, input$split_5, input$split_6))
      shiny::validate(
        shiny::need(all(!is.na(split_values) & split_values >= 1), "Clonal split values must be >= 1."),
        shiny::need(length(unique(split_values)) == 6, "Split values must be unique."),
        shiny::need(all(diff(split_values) >= 0), "Split values must be in increasing order.")
      )

      category_labels <- c(
        paste0("Unique (<= ", split_values[1], ")"),
        paste0("(", split_values[1], " < Ct <= ", split_values[2], ")"),
        paste0("(", split_values[2], " < Ct <= ", split_values[3], ")"),
        paste0("(", split_values[3], " < Ct <= ", split_values[4], ")"),
        paste0("(", split_values[4], " < Ct <= ", split_values[5], ")"),
        paste0("(", split_values[5], " < Ct <= ", split_values[6], ")"),
        paste0("(Ct > ", split_values[6], ")")
      )

      validate(need(group_col %in% colnames(df_orig), paste("Group column '", group_col, "' not found.")))

      df_processed <- df_orig %>% 
        dplyr::filter(!is.na(.data[[clone_call_method]])) %>%
        dplyr::group_by(.data[[group_col]], .data[[clone_call_method]]) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(.data[[group_col]]) %>%
        dplyr::mutate(category = cut(count,
                                     breaks = c(0, split_values, Inf),
                                     labels = category_labels,
                                     right = TRUE)) %>%
        dplyr::count(.data[[group_col]], category, name = "n") %>%
        dplyr::group_by(.data[[group_col]]) %>%
        dplyr::mutate(percentage = n / sum(n) * 100)

      return(df_processed)
    })

    output$plot <- renderPlot({
      df <- proportion_data(); req(df)
      ggplot(df, aes(x = .data[[input$group_by]], y = percentage, fill = category)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = input$group_by, y = "Percentage (%)") +
        theme_minimal() +
        theme(legend.position = input$legend_position)
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))

    output$table <- renderDT({
      df <- proportion_data(); req(df)
      datatable(df)
    })

    output$download_plot <- downloadHandler(
      filename = function() { "clonal_proportion_plot.pdf" },
      content = function(file) {
        pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
        print(
          ggplot(proportion_data(), aes(x = .data[[input$group_by]], y = percentage, fill = category)) +
            geom_bar(stat = "identity", position = "stack") +
            labs(x = input$group_by, y = "Percentage (%)") +
            theme_minimal() +
            theme(legend.position = input$legend_position)
        )
        dev.off()
      }
    )

    output$download_table <- downloadHandler(
      filename = function() { "clonal_proportion_data.csv" },
      content = function(file) {
        write.csv(proportion_data(), file, row.names = FALSE)
      }
    )
  })
}
