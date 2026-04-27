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
      uiOutput(ns("filter_groups_ui")),

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
      h5("Plot Options"),
      commonPlotOptions(ns, legend_selected = "right", width_value = 700, height_value = 500)
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      plotOutput(ns('plot')),
      hr(),
      h3("Table"),
      downloadButton(ns("download_table"), "Download Table (.xlsx)"),
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

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    observe({
      req(input$vdj_type, nzchar(input$vdj_type))

      vdj_df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(!is.null(vdj_df) && nrow(vdj_df) > 0)
      all_cols <- colnames(vdj_df)

      if (input$vdj_type == "tcr") {
        target_choices_all <- c(
          "TCR Paired CTnt" = "TCR_pair_CTnt", "TCR Paired CTaa" = "TCR_pair_CTaa",
          "TCR TRA CDR3 aa" = "TCR_TRA_cdr3", "TCR TRA CDR3 nt" = "TCR_TRA_cdr3_nt",
          "TCR TRB CDR3 aa" = "TCR_TRB_cdr3", "TCR TRB CDR3 nt" = "TCR_TRB_cdr3_nt"
        )
        selected_choice <- "TCR_TRB_cdr3"
      } else if (input$vdj_type == "bcr") {
        target_choices_all <- c(
          "BCR Paired CTnt" = "BCR_pair_CTnt", "BCR Paired CTaa" = "BCR_pair_CTaa",
          "BCR IGH CDR3 aa" = "BCR_IGH_cdr3_aa", "BCR IGH CDR3 nt" = "BCR_IGH_cdr3_nt",
          "BCR IGK CDR3 aa" = "BCR_IGK_cdr3_aa", "BCR IGK CDR3 nt" = "BCR_IGK_cdr3_nt",
          "BCR IGL CDR3 aa" = "BCR_IGL_cdr3_aa", "BCR IGL CDR3 nt" = "BCR_IGL_cdr3_nt"
        )
        selected_choice <- "BCR_IGH_cdr3_aa"
      } else {
        target_choices_all <- c("Select VDJ type first..." = "")
        selected_choice <- ""
      }

      valid_choices <- target_choices_all[target_choices_all %in% all_cols]
      if (length(valid_choices) == 0) {
        valid_choices <- c("適切な列が見つかりません" = "")
        selected_choice <- ""
      } else if (!selected_choice %in% valid_choices) {
        selected_choice <- valid_choices[1]
      }
      updateSelectInput(session, "target_cdr_column", choices = valid_choices, selected = selected_choice)
    }) |> bindEvent(input$vdj_type, myReactives$tcr_df, myReactives$bcr_df, ignoreNULL = FALSE, ignoreInit = FALSE)

    reactive_df_raw <- reactive({
      df <- NULL
      req(input$vdj_type)

      if (input$vdj_type == "tcr" && !is.null(myReactives$tcr_df)) {
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr" && !is.null(myReactives$bcr_df)) {
        df <- myReactives$bcr_df
      }
      req(df)

      if (!is.null(input$group_by) && nzchar(input$group_by) &&
          !input$group_by %in% names(df) && !is.null(myReactives$seurat_object)) {
        so_meta <- myReactives$seurat_object@meta.data
        if (input$group_by %in% colnames(so_meta)) {
          so_meta$barcode <- rownames(so_meta)
          meta_join <- so_meta[, c("barcode", input$group_by), drop = FALSE]
          df <- dplyr::left_join(df, meta_join, by = "barcode")
        }
      }

      return(df)
    })

    table <- reactive({
      req(reactive_df_raw(), input$target_cdr_column)

      target_col_sym <- sym(input$target_cdr_column)

      df <- reactive_df_raw() %>%
        filter(!is.na(!!target_col_sym) & !!target_col_sym != "") %>%
        mutate(
          length = nchar(as.character(!!target_col_sym))
        )

      df <- df %>% dplyr::select(input$group_by, length)
      df <- df %>%
        group_by(.data[[input$group_by]], length) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(.data[[input$group_by]]) %>%
        mutate(percentage = count / sum(count) * 100) %>%
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

    build_plot <- function() {
      req(table_filtered())
      angle     <- as.numeric(input$x_axis_angle %||% "45")
      hjust_val <- if (angle == 0) 0.5 else 1
      base_font <- input$base_font_size %||% 12
      leg_pos   <- input$legend %||% "right"

      p <- if (input$plot_type == "boxplot") {
        req(boxplot_df())
        boxplot_df() %>% ggplot(aes(x = .data[[input$group_by]], y = length, fill = .data[[input$group_by]])) +
          geom_boxplot() +
          labs(x = input$group_by, y = "Length", fill = input$group_by)
      } else {
        y_val <- if (input$display_value == "count") "count" else "percentage"
        p2 <- table_filtered() %>% ggplot(aes(x = length, y = .data[[y_val]], fill = .data[[input$group_by]])) +
          geom_bar(stat = 'identity', position = input$plot_type) +
          labs(x = "Length", y = tools::toTitleCase(y_val), fill = input$group_by) +
          scale_y_continuous(expand = c(0, 0))
        if (input$horizontal) p2 <- p2 + coord_flip()
        p2
      }

      p + theme_classic(base_size = base_font) +
        theme(
          axis.text.x    = element_text(angle = angle, hjust = hjust_val),
          legend.position = leg_pos
        )
    }

    output$plot <- renderPlot({
      build_plot()
    }, width = reactive(input$plot_width %||% 700), height = reactive(input$plot_height %||% 500))

    output$table <- renderDT({
      table_filtered()
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("cdr_length_", input$vdj_type, "_", input$target_cdr_column, "_", Sys.Date(), ".pptx")
      },
      content = function(file) {
        p <- build_plot()
        save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
      }
    )

    output$download_table <- downloadHandler(
      filename = function() {
        paste0("cdr_length_table_", input$vdj_type, "_", input$target_cdr_column, "_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        openxlsx::write.xlsx(table_filtered(), file)
      }
    )
  })
}
