# # UI部分
# highlightPlotUI <- function(id) {
#   ns <- NS(id)
#   sidebarLayout(
#     sidebarPanel(
#       selectInput(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne"), selected = "umap"),
#       selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
#       checkboxGroupInput(ns("unique_group"), "Select Group", choices = NULL, inline = TRUE),
#       sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
#       sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
#       downloadButton(ns("download_plot"), "Download plot (.pdf)")
#     ),
#     mainPanel(
#       plotOutput(ns("plot"))
#     )
#   )
# }
# 
# highlightPlotServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {
# 
#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       update_group_by_for_dimplot(session, myReactives)
#     })
# 
#     observeEvent(input$group_by, {
#       req(myReactives$seurat_object, input$group_by)
#       update_unique_group_choices(session, myReactives, input$group_by)
#     })
# 
#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       # group_by の選択肢を更新
#       update_group_by_for_dimplot(session, myReactives)
#       # 初回の unique_group の選択肢を更新（デフォルトの "sample" を基準に）
#       update_unique_group_choices(session, myReactives, "sample")
#     })
# 
#     plot <- reactive({
#       req(myReactives$seurat_object)
#       seurat_object <- myReactives$seurat_object
#       Idents(seurat_object) <- input$group_by
#       highlights <- WhichCells(seurat_object, idents = input$unique_group)
# 
#       DimPlot(
#         seurat_object,
#         reduction = input$reduction,
#         group.by = input$group_by,
#         cells.highlight = highlights,
#         cols.highlight = 'darkblue',
#         cols = 'grey'
#       )
#     })
# 
#     # プロットのレンダリング
#     output$plot <- renderPlot({
#       plot()
#     }, width = reactive(input$plot_width), height = reactive(input$plot_height))
# 
#     # PDFダウンロードハンドラー
#     output$download_plot <- downloadHandler(
#       filename = function() {
#         "UMAP_plot.pdf"
#       },
#       content = function(file) {
#         ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
#       }
#     )
#   })
# }
# 
