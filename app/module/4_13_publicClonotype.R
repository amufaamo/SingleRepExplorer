# Public clonotype analysis module

# --- UI ---
publicClonotypeUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      vdjType(ns),
      selectInput(ns("clonotype_col"), "Clonotype Column",
                  choices = NULL, selected = NULL),
      groupByInput(ns),

      h4("Venn Diagram Settings"),
      uiOutput(ns("group_selector_ui")),

      radioButtons(ns("table_unit"), "Table Value Unit",
                   choices = c("Count & Proportion" = "both", "Count" = "count", "Proportion" = "proportion"),
                   selected = "both"),
      hr(),
      actionButton(ns("run_venn"), "Run Venn Analysis", class = "btn-primary", width = "100%", icon = icon("play")),
      hr(),
      h5("Plot Options"),
      commonPlotOptions(ns, legend_selected = "right", width_value = 600, height_value = 600),
      # downloadButton moved to mainPanel
    ),
    mainPanel(
      width = 9,
      h3("Public Clonotypes Venn Diagram"),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns("venn_plot")),

      hr(),

      h3("Clonotype Details"),
      uiOutput(ns("intersection_selector_ui")),
      DTOutput(ns("clonotype_table"))
    )
  )
}

# --- Server ---
publicClonotypeServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    # Dynamically populate clonotype_col based on actual available columns
    observe({
      req(input$vdj_type)
      df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(df)
      prefix <- if (input$vdj_type == "tcr") "TCR_pair_" else "BCR_pair_"
      candidates <- c("raw_clonotype_id",
                      paste0(prefix, "CTaa"),
                      paste0(prefix, "CTgene"),
                      paste0(prefix, "CTnt"))
      available <- intersect(candidates, names(df))
      if (length(available) == 0) available <- names(df)[1]
      selected <- if ("raw_clonotype_id" %in% available) "raw_clonotype_id" else available[1]
      updateSelectInput(session, "clonotype_col", choices = available, selected = selected)
    })

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    reactive_df_raw <- reactive({
      req(input$vdj_type, input$group_by)
      df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      clono_col <- input$clonotype_col %||% "raw_clonotype_id"
      req(df, clono_col %in% names(df))

      # Join Seurat metadata if group_by column is not in the repertoire df
      group_col <- input$group_by
      if (!group_col %in% names(df)) {
        so <- myReactives$seurat_object
        if (!is.null(so) && group_col %in% names(so@meta.data) && "barcode" %in% names(df)) {
          meta_sub <- data.frame(barcode = rownames(so@meta.data), stringsAsFactors = FALSE)
          meta_sub[[group_col]] <- so@meta.data[[group_col]]
          df <- dplyr::left_join(df, meta_sub, by = "barcode")
        }
      }
      shiny::validate(shiny::need(group_col %in% names(df),
        paste("Column", group_col, "not found in repertoire data or Seurat metadata.")))

      df_filtered <- df %>%
        dplyr::filter(!is.na(.data[[group_col]]) & .data[[group_col]] != "")
      shiny::validate(shiny::need(nrow(df_filtered) > 0, paste("No data after removing NA/empty from", group_col)))
      df_filtered
    })

    output$group_selector_ui <- renderUI({
      df <- reactive_df_raw()
      req(df, input$group_by)
      available_groups <- sort(unique(df[[input$group_by]]))
      default_selected <- head(available_groups, 4)

      selectizeInput(session$ns("selected_groups"),
                     label = "Select groups to compare (2-4 recommended):",
                     choices = available_groups,
                     selected = default_selected,
                     multiple = TRUE,
                     options = list(maxItems = 4))
    })

    venn_input_list <- eventReactive(input$run_venn, {
      req(reactive_df_raw(), input$selected_groups, length(input$selected_groups) >= 2)

      df_venn <- reactive_df_raw() %>%
        dplyr::filter(.data[[input$group_by]] %in% input$selected_groups)

      clono_col <- input$clonotype_col %||% "raw_clonotype_id"
      split(df_venn[[clono_col]], df_venn[[input$group_by]]) %>%
        lapply(unique)
    })

    make_venn_plot <- reactive({
      req(venn_input_list())
      set_name_sz <- input$set_name_size %||% 4
      txt_sz      <- input$text_size %||% 5

      ggvenn::ggvenn(
        venn_input_list(),
        fill_color  = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
        stroke_size = 0.5,
        set_name_size = set_name_sz,
        text_size   = txt_sz
      ) +
      labs(title = paste("Clonotype Overlap between", tools::toTitleCase(gsub("_", " ", input$group_by)))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })

    output$venn_plot <- renderPlot({
      make_venn_plot()
    }, width  = reactive(input$plot_width  %||% 600),
       height = reactive(input$plot_height %||% 600))

    intersections <- eventReactive(input$run_venn, {
        req(venn_input_list())
        sets <- venn_input_list()
        set_names <- names(sets)

        combinations <- unlist(lapply(1:length(sets), function(i) combn(set_names, i, simplify = FALSE)), recursive = FALSE)

        intersection_list <- list()
        for (combo in combinations) {
            intersect_clones <- Reduce(intersect, sets[combo])
            other_sets <- set_names[!set_names %in% combo]
            if (length(other_sets) > 0) {
                clones_in_others <- unlist(sets[other_sets])
                final_clones <- setdiff(intersect_clones, clones_in_others)
            } else {
                final_clones <- intersect_clones
            }
            combo_name <- paste(combo, collapse = " & ")
            if (length(final_clones) > 0) {
              intersection_list[[combo_name]] <- final_clones
            }
        }
        return(intersection_list)
    })

    output$intersection_selector_ui <- renderUI({
      req(intersections())
      choices_with_counts <- purrr::map2_chr(names(intersections()), intersections(), ~ paste0(.x, " (", length(.y), " clones)"))
      selectInput(session$ns("selected_intersection"),
                  "Select intersection to inspect:",
                  choices = choices_with_counts)
    })

    output$clonotype_table <- renderDT({
      req(reactive_df_raw(), input$selected_intersection, intersections())
      selected_name <- sub(" \\(.*\\)$", "", input$selected_intersection)
      target_clonotypes <- intersections()[[selected_name]]

      clono_col <- input$clonotype_col %||% "raw_clonotype_id"
      table_data <- reactive_df_raw() %>%
        dplyr::filter(
          .data[[input$group_by]] %in% input$selected_groups,
          .data[[clono_col]] %in% target_clonotypes
        ) %>%
        dplyr::count(.data[[input$group_by]], .data[[clono_col]], name = "count") %>%
        dplyr::group_by(.data[[input$group_by]]) %>%
        dplyr::mutate(proportion = count / sum(count)) %>%
        dplyr::ungroup() %>%
        dplyr::rename(group = .data[[input$group_by]]) %>%
        tidyr::pivot_wider(
            names_from  = group,
            values_from = c(count, proportion),
            values_fill = 0
        )

      datatable(table_data, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    })

    output$download_plot <- downloadHandler(
      filename = function() { paste0("venn_clonotype_", input$group_by, "_", Sys.Date(), ".pptx") },
      content = function(file) {
        p <- make_venn_plot()
        save_plot_as_pptx(file, p, input$plot_width %||% 600, input$plot_height %||% 600)
      }
    )
  })
}
