library(shiny)
library(Seurat)
library(DT)
library(dplyr)
library(plotly)
library(ggplot2)
library(shinyjs)



differentialGeneExpressionUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      selectInput(
        ns("analysis_type"),
        label = "Comparison Type",
        choices = c(
          "To entire dataset" = "all",
          "Between selected group(s)" = "two"
        ),
        selected = "all"
      ),
      selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
      conditionalPanel(
        condition = "input.analysis_type == 'two'",
        ns = ns,
        selectInput(ns("target_cluster"), label = "Target cluster(s)", choices = NULL, selected = NULL, multiple = TRUE),
        selectInput(ns("reference_cluster"), label = "Reference cluster(s)", choices = NULL, selected = NULL, multiple = TRUE),
      ),
      numericInput(ns("logfc"), "Threshold of Log fold change (default: 0.1)", min = 0, max = 1, value = 0.1, step = 0.01),
      numericInput(ns("minpct"), "Minimum percentage of cells expressing a gene (%) (default: 0.01)", min = 0, max = 1, value = 0.01, step = 0.01),
      conditionalPanel(
        condition = "input.analysis_type == 'all'",
        ns = ns,
        selectInput(
          ns("choice"),
          "Show Features",
          choices = list(
            "Only Up-regulated" = "positive",
            "Only Down-regulated" = "negative",
            "Both" = "both"
          ),
          selected = "positive"
        ),
      ),
      actionButton(ns("run"), "Run"),
      conditionalPanel(
        condition = "input.analysis_type == 'all' && output.table_markers",
        ns = ns,
        numericInput(ns("num"), label = "Number of markers", value = 10),
        numericInput(ns("p_val_adj_marker"), "Threshold of adjusted p-value (default: 0.05)", min = 0, max = 1, value = 0.05, step = 0.01),
      ),
      conditionalPanel(
        condition = "input.analysis_type == 'two' && output.table_two",
        ns = ns,
        textInput(ns("gene_input"), "Highlight Genes (comma separated)", value = ""),
        numericInput(ns("p_val_adj"), "Threshold of adjusted p-value (default: 0.05)", min = 0, max = 1, value = 0.05, step = 0.01),
        numericInput(ns("highlight_size"), "Highlight Size", min = 0.001, max = 10, value = 2),
        numericInput(ns("other_size"), "Other Points Size", min = 0.001, max = 10, value = 1),
        numericInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 400, step = 100),
        numericInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 400, step = 100),
      ),
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.analysis_type == 'all'",
        ns = ns,
        h3("Most Differentially Expressed Features (markers) for Each Group"),
        downloadButton(ns("download_marker_table"), "Download marker table (.csv)"),
        DTOutput(ns("table_markers")),
        hr(),
        h3("Differentially Expressed Features across Groups"),
        downloadButton(ns("download_table_all"), "Download all table (.csv)"),
        DTOutput(ns("table_all")),
      ),
      conditionalPanel(
        condition = "input.analysis_type == 'two'",
        ns = ns,
        h3("Volcano Plot"),
        downloadButton(ns("download_plot"), "Download Plot (PDF)"),
        plotlyOutput(ns("volcanoPlot")),
        hr(),
        h3("Differentially Expressed Features Between Selected Group(s)"),
        downloadButton(ns("download_table_two"), "Download all table (.csv)"), # Keep using all
        DTOutput(ns("table_two")),
      )
    )
  )
}



differentialGeneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    # Initialize disabled state
    shinyjs::disable("download_marker_table")
    shinyjs::disable("download_table_all")
    shinyjs::disable("download_table_two")


    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_for_marker(session, input, myReactives)
    })


    # twoのときに、クラスターの種類
    observe({
      req(myReactives$seurat_object, input$group_by, input$analysis_type == "two")
      all_clusters <- unique(myReactives$seurat_object@meta.data[[input$group_by]])

      updated_target_choices <- setdiff(all_clusters, input$reference_cluster)
      numeric_choices_target <- as.numeric(updated_target_choices)
      sorted_choices_target <- if (all(!is.na(numeric_choices_target))) {
        sort(numeric_choices_target)
      } else {
        sort(updated_target_choices)
      }
      updateSelectInput(session, "target_cluster", choices = sorted_choices_target, selected = input$target_cluster)

      updated_reference_choices <- setdiff(all_clusters, input$target_cluster)
      numeric_choices_reference <- as.numeric(updated_reference_choices)
      sorted_choices_reference <- if (all(!is.na(numeric_choices_reference))) {
        sort(numeric_choices_reference)
      } else {
        sort(updated_reference_choices)
      }
      updateSelectInput(session, "reference_cluster", choices = sorted_choices_reference, selected = input$reference_cluster)
    })


    observeEvent(input$run, {
      req(myReactives$seurat_object, input$group_by)
      so <- myReactives$seurat_object
      Idents(so) <- input$group_by

      if (input$analysis_type == "all") {
        withProgress(message = "Calculating...", value = 0, {
          clusters <- unique(Idents(so))
          n_clusters <- length(clusters)
          all_data_list <- list()

          for (i in seq_along(clusters)) {
            cluster <- clusters[i]
            incProgress(1 / n_clusters, detail = paste("Processing cluster", cluster))

            # Use FindMarkers with only.pos based on input$choice
            markers <- FindMarkers(
              so,
              ident.1 = cluster,
              logfc.threshold = input$logfc,
              min.pct = input$minpct,
              only.pos = (input$choice == "positive"), # Correct only.pos handling
              verbose = FALSE
            )
            if (input$choice == "negative") {
              markers <- markers[markers$avg_log2FC < 0, ]
            } else if (input$choice == "both") {
              # No filter
            } else {
              markers <- markers[markers$avg_log2FC > 0, ]
            }
            markers$cluster <- cluster
            markers$gene <- rownames(markers)
            all_data_list[[i]] <- markers
          }
          myReactives$all_data_all <- do.call(rbind, all_data_list)
        })
      } else if (input$analysis_type == "two") {
        req(input$target_cluster, input$reference_cluster)
        withProgress(message = "Calculating...", value = 0, {
          myReactives$all_data_two <- FindMarkers(
            so,
            ident.1 = input$target_cluster,
            ident.2 = input$reference_cluster,
            logfc.threshold = input$logfc,
            min.pct = input$minpct
          )
          myReactives$all_data_two$gene <- rownames(myReactives$all_data_two) # Add gene names

          incProgress(1) # Full progress for single comparison
        })
      }
    })


    # Marker Table (analysis_type == "all")
    observe({
      req(input$analysis_type == "all")
      req(myReactives$all_data_all)

      myReactives$marker_table <- myReactives$all_data_all %>%
        dplyr::filter(p_val_adj <= input$p_val_adj_marker) %>%
        group_by(cluster) %>%
        arrange(desc(abs(avg_log2FC))) %>%
        slice_head(n = input$num) %>%
        summarise(genes = paste(gene, collapse = ", "), .groups = "drop")
    })

    output$table_markers <- renderDT({
      req(myReactives$marker_table)
      datatable(
        myReactives$marker_table,
        options = list(paging = FALSE)
      )
    })

    # All Data Table (both analysis types)
    output$table_all <- renderDT({
      req(myReactives$all_data_all) # Use all_data
      datatable(myReactives$all_data_all,
        filter = "top",
        options = list(pageLength = 10)
      )
    })

    # All Data Table (both analysis types)
    output$table_two <- renderDT({
      req(myReactives$all_data_two) # Use all_data
      datatable(myReactives$all_data_two,
        filter = "top",
        options = list(pageLength = 10)
      )
    })

    volcano_plot_data <- reactive({
      req(myReactives$all_data_two)

      filtered_data <- myReactives$all_data_two %>%
        dplyr::filter(p_val_adj <= input$p_val_adj)

      volcano_data <- filtered_data %>%
        mutate(log_p_value = -log10(p_val_adj))

      highlighted <- if (!is.null(input$gene_input) && nzchar(input$gene_input)) {
        genes <- strsplit(input$gene_input, ",\\s*")[[1]]
        tolower(volcano_data$gene) %in% tolower(genes) # Use the gene column
      } else {
        rep(FALSE, nrow(volcano_data))
      }

      num_highlighted <- sum(highlighted)

      p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log_p_value, text = gene)) + # Use the gene column
        geom_point(
          aes(
            color = ifelse(highlighted, "Highlighted", "Other"),
            size = ifelse(highlighted, input$highlight_size, input$other_size)
          ),
          alpha = 0.7
        ) +
        scale_color_manual(
          values = c("Highlighted" = "red", "Other" = "gray"),
          labels = c("Highlighted" = paste0("Highlighted (", num_highlighted, ")"), "Other" = "Other")
        ) +
        scale_size_identity() +
        geom_text(
          data = subset(volcano_data, highlighted),
          aes(label = gene), # Use the gene column
          size = 3,
          vjust = -0.5,
          hjust = -0.1,
          check_overlap = TRUE
        ) +
        labs(
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted p-value"
        ) +
        theme_classic() +
        theme(legend.position = "none")

      return(p)
    })

    output$volcanoPlot <- renderPlotly({
      volcano_plot_data() %>%
        ggplotly(tooltip = "text", width = input$plot_width, height = input$plot_height) %>%
        layout(
          plot_bgcolor = "white",
          panel_grid = list(visible = FALSE),
          xaxis = list(showline = TRUE, linecolor = "black"),
          yaxis = list(showline = TRUE, linecolor = "black"),
          showlegend = FALSE
        )
    })

    output$download_marker_table <- downloadHandler(
      filename = function() {
        "marker_table.csv"
      },
      content = function(file) {
        write.csv(myReactives$marker_table, file, row.names = FALSE)
      }
    )

    output$download_table_all <- downloadHandler(
      filename = function() {
        "all_table.csv"
      },
      content = function(file) {
        write.csv(myReactives$all_data_all, file, row.names = TRUE) # Use all_data
      }
    )

    output$download_table_two <- downloadHandler(
      filename = function() {
        "all_table.csv"
      },
      content = function(file) {
        write.csv(myReactives$all_data_two, file, row.names = TRUE) # Use all_data
      }
    )

   # PDF Download Handler
    output$download_plot <- downloadHandler(
      filename = function() {
        "Volcano_plot.pdf"
      },
      content = function(file) {
        ggsave(file, plot = volcano_plot_data(), device = "pdf", width = input$plot_width / 72, height = input$plot_height / 72, units = "in", dpi = 300)
      }
    )


    # Enable/disable download buttons based on table content
    observe({
      if (input$analysis_type == "all") {
        if (is.null(myReactives$marker_table) || nrow(myReactives$marker_table) == 0) {
          shinyjs::disable("download_marker_table")
        } else {
          shinyjs::enable("download_marker_table")
        }
        if (is.null(myReactives$all_data_all) || nrow(myReactives$all_data_all) == 0) {
          shinyjs::disable("download_table_all")
        } else {
          shinyjs::enable("download_table_all")
        }
      } else if (input$analysis_type == "two") {
        if (is.null(myReactives$all_data_two) || nrow(myReactives$all_data_two) == 0) {
          shinyjs::disable("download_table_two")
          shinyjs::disable("download_table_all")
        } else {
          shinyjs::enable("download_table_two")
          shinyjs::enable("download_table_all")
        }
      }
    })
  })
}
