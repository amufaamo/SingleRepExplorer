# qualityControlUI <- function(id) {
#   ns <- NS(id)
#   sidebarLayout(
#     sidebarPanel(
#       sliderInput(ns("slider_nCount_RNA"), label = "cCount", min = 0, max = 100000, value = c(0, 100000)),
#       sliderInput(ns("slider_nFeatures_RNA"), label = "nFeature", min = 0, max = 10000, value = c(0, 10000)),
#       sliderInput(ns("slider_percent_mt"), label = "percent mt", min = 0, max = 100, value = 100),
#       sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 1, value = 0.1, step = 0.01),
#       sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 900, step = 100),
#       sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 900, step = 100),
#       actionButton(ns('rerun'), label = "Rerun"),
# #      textOutput(ns("filtered_cell_ratio"))
#     ),
#     mainPanel(
#       div(style = "font-size: 30px;",
#         textOutput(ns("filtered_cell_ratio"))
#       ),
#       plotOutput(ns("plot")),
#     )
#   )
# }

# qualityControlServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {

#     # Update input
#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       seurat_object <- myReactives$seurat_object
#       updateSliderInput(session, "slider_nCount_RNA", min = min(seurat_object@meta.data$nCount_RNA), max = max(seurat_object@meta.data$nCount_RNA), value = c(min(seurat_object@meta.data$nCount_RNA), max(seurat_object@meta.data$nCount_RNA)))
#       updateSliderInput(session, "slider_nFeatures_RNA", min = min(seurat_object@meta.data$nFeature_RNA), max = max(seurat_object@meta.data$nFeature_RNA), value = c(min(seurat_object@meta.data$nFeature_RNA), max(seurat_object@meta.data$nFeature_RNA)))
#       updateSliderInput(session, "slider_percent_mt", min = 0, max = floor(max(seurat_object@meta.data$percent.mt)) + 1, value = floor(max(seurat_object@meta.data$percent.mt)) + 1)
#     })

#      observeEvent({
#       input$slider_nCount_RNA
#       input$slider_nFeatures_RNA
#       input$slider_percent_mt
#     }, {
#       req(myReactives$seurat_object)
#       seurat_object <- myReactives$seurat_object
      
#       total_cells <- nrow(seurat_object@meta.data)
      
#       filtered_seurat <- subset(seurat_object, 
#                                 subset = nFeature_RNA >= input$slider_nFeatures_RNA[1] & 
#                                          nFeature_RNA <= input$slider_nFeatures_RNA[2] & 
#                                          nCount_RNA >= input$slider_nCount_RNA[1] & 
#                                          nCount_RNA <= input$slider_nCount_RNA[2] & 
#                                          percent.mt <= input$slider_percent_mt)
      
#       filtered_cells <- nrow(filtered_seurat@meta.data)
#       ratio <- round(filtered_cells / total_cells, 3)
      
#       output$filtered_cell_ratio <- renderText({
#         paste0("Filtered Cells: ", filtered_cells, " / ", total_cells, " (", ratio * 100, "%)")
#       })
#     })


      
#     plot <- reactive({
#       req(myReactives$seurat_object)
#       myReactives$quality_sca <- myquality_control_scatterplot(myReactives$seurat_object, input$slider_nFeatures_RNA[1], input$slider_nFeatures_RNA[2], input$slider_nCount_RNA[1], input$slider_nCount_RNA[2], input$slider_percent_mt, input$point_size)
#       myReactives$quality_vln <- myquality_control_violinplot(myReactives$seurat_object, input$slider_nFeatures_RNA[1], input$slider_nFeatures_RNA[2], input$slider_nCount_RNA[1], input$slider_nCount_RNA[2], input$slider_percent_mt)
#       ggarrange(myReactives$quality_sca, myReactives$quality_vln, ncol = 1) + theme(legend.position = input$legend)
#       })
    
#     # プロットのレンダリング
#     output$plot <- renderPlot({
#       plot()
#     }, width = reactive(input$plot_width), height = reactive(input$plot_height))
    
#     # PDFダウンロードハンドラー
#     output$download_plot <- downloadHandler(
#       filename = function() { "UMAP_plot.pdf" },
#       content = function(file) {
#         ggsave(file, plot = renderCustomPlot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
#       }
#     )
  

#     observeEvent(input$rerun, {
#       myReactives <- myseurat_normalize_umap_rerun(myReactives, input$slider_nFeatures_RNA[1], input$slider_nFeatures_RNA[2], input$slider_nCount_RNA[1], input$slider_nCount_RNA[2], input$slider_percent_mt)
# #      myReactives <- run_tcr(myReactives)
# #      myReactives <- run_bcr(myReactives)
      
#     })
#   })
# }


# myquality_control_scatterplot <- function(seurat_object, feature_low, feature_high, count_low, count_high, mito, plotsize){
  
#   plot1 <- seurat_object@meta.data %>% 
#     ggplot(mapping = aes(x = nCount_RNA, y = percent.mt)) + 
#     geom_pointdensity(size = plotsize) +
#     scale_color_viridis(guide = 'none') +
#     theme(legend.position = 'none') +
#     theme_classic() +
#     geom_vline(xintercept=count_low, col="red") +
#     geom_vline(xintercept=count_high, col="red") +
#     geom_hline(yintercept=mito, col="red")
  
#   plot1 <- ggMarginal(plot1, 
#                       type = 'density', 
#                       margins = 'both', 
#                       fill = 'red')
  
#   plot2 <- seurat_object@meta.data %>% 
#     ggplot(mapping = aes(x = nCount_RNA, y = nFeature_RNA)) + 
#     geom_pointdensity(size = plotsize) +
#     scale_color_viridis(guide = 'none') +
#     theme(legend.position = 'none') +
#     theme_classic() +
#     geom_vline(xintercept=count_low, col="red") +
#     geom_vline(xintercept=count_high, col="red") +
#     geom_hline(yintercept=feature_low, col="red") +
#     geom_hline(yintercept=feature_high, col="red")
  
#   plot2 <- ggMarginal(plot2, 
#                       type = 'density', 
#                       margins = 'both',
#                       fill = 'red')
  
  
#   plot <- ggarrange(plot1, plot2, ncol = 2)
#   return(plot)
# }

# myquality_control_violinplot <- function(seurat_object, feature_low, feature_high, count_low, count_high, mito){
#   vln1 <- VlnPlot(seurat_object, features = "nFeature_RNA", group.by='orig.ident') + theme(legend.position = "none") + geom_hline(yintercept=feature_low, col="red") + geom_hline(yintercept=feature_high, col="red")
#   vln2 <- VlnPlot(seurat_object, features = "nCount_RNA", group.by='orig.ident') + theme(legend.position = "none") + geom_hline(yintercept=count_low, col="red") + geom_hline(yintercept=count_high, col="red")
#   vln3 <- VlnPlot(seurat_object, features = "percent.mt", group.by='orig.ident') + theme(legend.position = "none") + geom_hline(yintercept=mito, col="red")
#   plot <- ggarrange(vln1, vln2, vln3, ncol = 3)
#   return(plot)
# }


# myseurat_normalize_umap_rerun <- function(myReactives, low_nFeature, high_nFeature, low_Count, high_Count, mito){
#   seurat_object <- myReactives$seurat_object
#   seurat_object <- subset(seurat_object, subset = nFeature_RNA > low_nFeature)
#   seurat_object <- subset(seurat_object, subset = nFeature_RNA < high_nFeature)
#   seurat_object <- subset(seurat_object, subset = nCount_RNA > low_Count)
#   seurat_object <- subset(seurat_object, subset = nCount_RNA < high_Count)
#   seurat_object <- subset(seurat_object, subset = percent.mt < mito)
#   seurat_object <- NormalizeData(seurat_object)
#   seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
#   all.genes <- rownames(seurat_object)
#   seurat_object <- ScaleData(seurat_object, features = all.genes)
#   seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), npcs = 50)
#   seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
#   seurat_object <- FindClusters(seurat_object, resolution = 0.5)
#   seurat_object <- RunUMAP(seurat_object, dims = 1:10)
#   seurat_object <- RunTSNE(seurat_object, dims = 1:10)
#   myReactives$seurat_object <- seurat_object
#   return(myReactives)
  
# }



