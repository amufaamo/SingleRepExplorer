# findMarkersUI <- function(id) {
#   ns <- NS(id)
#   sidebarLayout(
#     sidebarPanel(
#       radioButtons(ns("analysis_type"),
#         label = "Analysis Type",
#         choices = c("Find All Markers" = "all", "Find Markers (Two Groups)" = "two"),
#         selected = "all"
#       ),
#       radioButtons(ns("group_by"),
#         label = "Group By",
#         choices = c("sample", "seurat_clusters"),
#         selected = "seurat_clusters"
#       ),
#       uiOutput(ns("dynamic_ui")),
#       actionButton(ns("run"), "Calculate Markers"),
#     ),
#     mainPanel(
#       plotOutput(ns("volcano_plot")),
#       DTOutput(ns("table_markers")),
#       numericInput(ns("num"), label = "Number of markers", value = 10),
#       radioButtons(ns("choice"), "Positive or negative markers", choices = list("positive", "negative", "both"), selected = "positive"),
#       sliderInput(ns("logfc"),
#         "Threshold of Log fold change (default: 0.1)",
#         min = 0,
#         max = 1,
#         value = 0.1,
#         step = 0.01
#       ),
#       sliderInput(ns("minpct"),
#         "Minimum percentage of cells expressing a gene (%) (default: 0.01)",
#         min = 0,
#         max = 1,
#         value = 0.01,
#         step = 0.01
#       ),
#       sliderInput(ns("p_val_adj"),
#         "Threshold of adjusted p-value (default: 0.05)",
#         min = 0,
#         max = 1,
#         value = 0.05,
#         step = 0.01
#       ),
#       downloadButton(ns("download_table"), "Download table (.csv)"),
#       DTOutput(ns("table"))
#     )
#   )
# }

# findMarkersServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {
#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       update_group_by_for_marker(session, input, myReactives) # input を渡す
#     })
#     output$dynamic_ui <- renderUI({
#       ns <- session$ns
#       if (input$analysis_type == "two") {
#         tagList(
#           selectInput(ns("target_cluster"), "Target Cluster", choices = NULL),
#           selectInput(ns("reference_cluster"), "Reference Cluster", choices = NULL)
#         )
#       }
#     })

#     observeEvent(input$run, {
#       req(myReactives$seurat_object)
#       so <- myReactives$seurat_object
#       Idents(so) <- input$group_by

#       if (input$analysis_type == "all") {
#         myReactives$markers_table <- FindAllMarkers(so, logfc.threshold = 0, min.pct = 0)
#       } else {
#         req(input$target_cluster, input$reference_cluster)
#         myReactives$markers_table <- FindMarkers(so, ident.1 = input$target_cluster, ident.2 = input$reference_cluster)
#       }
#     })

#     output$volcano_plot <- renderPlot({
#       req(myReactives$markers_table)
#       ggplot(myReactives$markers_table, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
#         geom_point(aes(color = p_val_adj < input$p_val_adj & abs(avg_log2FC) >= input$logfc)) +
#         scale_color_manual(values = c("black", "red")) +
#         labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
#         theme_minimal()
#     })

#     output$table_markers <- renderDT({
#       req(myReactives$markers_table)
#       filtered_markers <- switch(input$choice,
#         positive = myReactives$markers_table %>% filter(avg_log2FC > 0),
#         negative = myReactives$markers_table %>% filter(avg_log2FC < 0),
#         both = myReactives$markers_table
#       )
#       filtered_markers <- filtered_markers %>%
#         dplyr::filter(
#           abs(avg_log2FC) >= input$logfc,
#           pct.1 >= input$minpct,
#           pct.2 >= input$minpct,
#           p_val_adj <= input$p_val_adj
#         )
#       table <- filtered_markers %>%
#         group_by(cluster) %>%
#         arrange(desc(abs(avg_log2FC))) %>%
#         slice_head(n = input$num) %>%
#         summarise(genes = paste(gene, collapse = ", "), .groups = "drop")
#       return(table)
#     })

#     output$download_table <- downloadHandler(
#       filename = function() {
#         "table.csv"
#       },
#       content = function(file) {
#         write.csv(myReactives$markers_table, file)
#       }
#     )
#   })
# }
