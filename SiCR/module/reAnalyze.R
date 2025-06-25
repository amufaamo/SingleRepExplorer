# # 統合されたUI (単一画面レイアウト)
# integratedQCSubsetUI <- function(id) {
#   ns <- NS(id)

#   sidebarLayout(
#     sidebarPanel(
#       selectInput(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne"), selected = "umap"),
#       selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
#       checkboxGroupInput(ns("unique_group"), "Select Group", choices = NULL, inline = TRUE),
#       sliderInput(ns("slider_nCount_RNA"), label = "nCount_RNA", min = 0, max = 100000, value = c(0, 100000)),
#       sliderInput(ns("slider_nFeatures_RNA"), label = "nFeature_RNA", min = 0, max = 10000, value = c(0, 10000)),
#       sliderInput(ns("slider_percent_mt"), label = "percent.mt", min = 0, max = 100, value = 100),
#       sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 10, value = 0.1, step = 0.01),
#       sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
#       sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100)
#     ),
#     mainPanel(
#       div(
#         style = "font-size: 30px;",
#         textOutput(ns("filtered_cell_ratio")),
#         plotOutput(ns("umap_tsne_plot")), # UMAP/t-SNE plot
#         plotOutput(ns("qc_plot")) # QC plots
#       ),
#     ),
#   )
# }
# # # 統合されたUI (単一画面レイアウト)
# # integratedQCSubsetUI <- function(id) {
# #   ns <- NS(id)

# #   fluidPage(
# #     # Subset Section
# #     fluidRow(
# #       column(
# #         12,
# #         div(
# #           style = "font-size: 30px;",
# #           textOutput(ns("filtered_cell_ratio"))
# #         ),
# #       ),
# #     ),
# #     fluidRow(
# #       div(
# #         style = "font-size: 30px;",
# #         textOutput(ns("filtered_cell_ratio")) # Filtered cell ratio above the plots
# #       ),
# #     ),
# #     fluidRow(
# #       column(
# #         4, # Subset Sidebar
# #         selectInput(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne"), selected = "umap"),
# #         selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
# #         checkboxGroupInput(ns("unique_group"), "Select Group", choices = NULL, inline = TRUE),
# #       ),
# #       column(
# #         4, # Subset Main Panel
# #         sliderInput(ns("slider_nCount_RNA"), label = "nCount_RNA", min = 0, max = 100000, value = c(0, 100000)),
# #         sliderInput(ns("slider_nFeatures_RNA"), label = "nFeature_RNA", min = 0, max = 10000, value = c(0, 10000)),
# #         sliderInput(ns("slider_percent_mt"), label = "percent.mt", min = 0, max = 100, value = 100),
# #       ),
# #       column(
# #         4,
# #         sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 10, value = 0.1, step = 0.01),
# #         sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
# #         sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
# #       ),
# #     ),
# #     fluidRow(
# #       column(
# #         6,
# #         plotOutput(ns("umap_tsne_plot")), # UMAP/t-SNE plot
# #       ),
# #       column(
# #         6,
# #         plotOutput(ns("qc_plot")) # QC plots
# #       ),
# #     ),
# #   )
# # }

# # 統合されたServer (UMAP/t-SNE plotのwidth/height reactives 削除, QC plot width/height reactives 削除)
# integratedQCSubsetServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {
#     # --- Subset Logic ---

#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       update_group_by_for_dimplot(session, myReactives)
#     })

#     observeEvent(input$group_by, {
#       req(myReactives$seurat_object, input$group_by)
#       update_unique_group_choices(session, myReactives, input$group_by)
#     })

#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       # group_by の選択肢を更新
#       update_group_by_for_dimplot(session, myReactives)
#       # 初回の unique_group の選択肢を更新（デフォルトの "sample" を基準に）
#       update_unique_group_choices(session, myReactives, "sample")
#     })

#     # Reactive expression to filter the Seurat object
#     filtered_seurat <- reactive({
#       req(myReactives$seurat_object, input$unique_group, input$group_by)
#       so <- myReactives$seurat_object
#       # フィルタリング対象のセルを取得
#       cells_to_keep <- rownames(so@meta.data)[so@meta.data[[input$group_by]] %in% input$unique_group]
#       # フィルタされたオブジェクトを作成
#       so_filtered <- subset(so, cells = cells_to_keep)
#       return(so_filtered)
#     })

#     umap_tsne_plot <- reactive({
#       req(filtered_seurat()) # filtered_seurat を正しく取得
#       so <- filtered_seurat()
#       print(dim(so@meta.data)) # デバッグ: メタデータの次元を確認
#       DimPlot(
#         so,
#         reduction = input$reduction,
#         group.by = input$group_by,
#         pt.size = input$point_size
#       ) + theme(legend.position = "none")
#     })

#     output$umap_tsne_plot <- renderPlot(
#       {
#         umap_tsne_plot()
#       },
#       width = reactive(input$plot_width),
#       height = reactive(input$plot_height)
#     )

#     # Update input
#     observeEvent(filtered_seurat(), {
#       req(filtered_seurat())
#       seurat_object <- filtered_seurat()
#       updateSliderInput(session, "slider_nCount_RNA", min = min(seurat_object@meta.data$nCount_RNA), max = max(seurat_object@meta.data$nCount_RNA), value = c(min(seurat_object@meta.data$nCount_RNA), max(seurat_object@meta.data$nCount_RNA)))
#       updateSliderInput(session, "slider_nFeatures_RNA", min = min(seurat_object@meta.data$nFeature_RNA), max = max(seurat_object@meta.data$nFeature_RNA), value = c(min(seurat_object@meta.data$nFeature_RNA), max(seurat_object@meta.data$nFeature_RNA)))
#       updateSliderInput(session, "slider_percent_mt", min = 0, max = floor(max(seurat_object@meta.data$percent.mt)) + 1, value = floor(max(seurat_object@meta.data$percent.mt)) + 1)
#     })

#     quality_plot <- reactive({
#       req(filtered_seurat())
#       myReactives$quality_sca <- myquality_control_scatterplot(filtered_seurat(), input$slider_nFeatures_RNA[1], input$slider_nFeatures_RNA[2], input$slider_nCount_RNA[1], input$slider_nCount_RNA[2], input$slider_percent_mt, input$point_size, input$group_by)
#       myReactives$quality_vln <- myquality_control_violinplot(filtered_seurat(), input$slider_nFeatures_RNA[1], input$slider_nFeatures_RNA[2], input$slider_nCount_RNA[1], input$slider_nCount_RNA[2], input$slider_percent_mt, input$point_size)
#       ggarrange(myReactives$quality_sca, myReactives$quality_vln, ncol = 1)
#     })

#     # プロットのレンダリング
#     output$qc_plot <- renderPlot(
#       {
#         quality_plot()
#       },
#       width = reactive(input$plot_width),
#       height = reactive(input$plot_height)
#     )
#   })


#   # observeEvent(
#   #   {
#   #     input$slider_nCount_RNA
#   #     input$slider_nFeatures_RNA
#   #     input$slider_percent_mt
#   #   },
#   #   {
#   #     req(filtered_seurat())
#   #     seurat_object <- filtered_seurat()

#   #     total_cells <- nrow(seurat_object@meta.data)

#   #     filtered_seurat <- subset(seurat_object,
#   #       subset = nFeature_RNA >= input$slider_nFeatures_RNA[1] &
#   #         nFeature_RNA <= input$slider_nFeatures_RNA[2] &
#   #         nCount_RNA >= input$slider_nCount_RNA[1] &
#   #         nCount_RNA <= input$slider_nCount_RNA[2] &
#   #         percent.mt <= input$slider_percent_mt
#   #     )

#   #     filtered_cells <- nrow(filtered_seurat@meta.data)
#   #     ratio <- round(filtered_cells / total_cells, 3)

#   #     output$filtered_cell_ratio <- renderText({
#   #       paste0("Filtered Cells: ", filtered_cells, " / ", total_cells, " (", ratio * 100, "%)")
#   #     })
#   #   }
#   # )
# }


# myquality_control_scatterplot <- function(seurat_object, feature_low, feature_high, count_low, count_high, mito, plotsize, group_by) {
#   plot1 <- seurat_object@meta.data %>%
#     ggplot(mapping = aes(x = nCount_RNA, y = percent.mt, color = .data[[group_by]])) +
#     geom_point(size = plotsize) +
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_vline(xintercept = count_low, col = "red") +
#     geom_vline(xintercept = count_high, col = "red") +
#     geom_hline(yintercept = mito, col = "red")

#   plot1 <- ggMarginal(plot1,
#     type = "density",
#     margins = "both",
#     fill = "red"
#   )

#   plot2 <- seurat_object@meta.data %>%
#     ggplot(mapping = aes(x = nCount_RNA, y = nFeature_RNA, color = .data[[group_by]])) +
#     geom_point(size = plotsize) +
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_vline(xintercept = count_low, col = "red") +
#     geom_vline(xintercept = count_high, col = "red") +
#     geom_hline(yintercept = feature_low, col = "red") +
#     geom_hline(yintercept = feature_high, col = "red")

#   plot2 <- ggMarginal(plot2,
#     type = "density",
#     margins = "both",
#     fill = "red"
#   )

#   plot <- ggarrange(plot1, plot2, ncol = 2)
#   return(plot)
# }

# myquality_control_violinplot <- function(seurat_object, feature_low, feature_high, count_low, count_high, mito, plot_size) {
#   vln1 <- VlnPlot(seurat_object, features = "nFeature_RNA", group.by = "orig.ident", pt.size = plot_size) + theme(legend.position = "none") + geom_hline(yintercept = feature_low, col = "red") + geom_hline(yintercept = feature_high, col = "red")
#   vln2 <- VlnPlot(seurat_object, features = "nCount_RNA", group.by = "orig.ident", pt.size = plot_size) + theme(legend.position = "none") + geom_hline(yintercept = count_low, col = "red") + geom_hline(yintercept = count_high, col = "red")
#   vln3 <- VlnPlot(seurat_object, features = "percent.mt", group.by = "orig.ident", pt.size = plot_size) + theme(legend.position = "none") + geom_hline(yintercept = mito, col = "red")
#   plot <- ggarrange(vln1, vln2, vln3, ncol = 3)
#   return(plot)
# }


# myseurat_normalize_umap_rerun <- function(myReactives, low_nFeature, high_nFeature, low_Count, high_Count, mito) {
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
