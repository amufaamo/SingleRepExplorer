# # UI部分
# subsetUI <- function(id) {
#   ns <- NS(id)
#   sidebarLayout(
#     sidebarPanel(
#       #      radioButtons(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
#       #      radioButtons(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne", "PCA" = "pca"), selected = "umap"),
#       selectInput(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne"), selected = "umap"),
#       selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
#       checkboxGroupInput(ns("unique_group"), "Select Group", choices = NULL), # checkboxGroupInputに変更
#       #      selectInput(ns("split_by"), "Split by", choices = c("sample", "seurat_clusters"), selected = "sample"),
#       sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 10, value = 0.1, step = 0.01),
#       sliderInput(ns("label_size"), "Size of labels", min = 0, max = 20, value = 10, step = 1),
#       radioButtons(ns("legend"), "Legend", choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
#       sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
#       sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
#       downloadButton(ns("download_plot"), "Download plot (.pdf)"),
#       downloadButton(ns("download_table"), "Download Table (.csv)") # ダウンロードボタン
#     ),
#     mainPanel(
#       plotOutput(ns("plot")),
#       actionButton(ns('rerun'), "Re-run Clustering"),
#       DT::dataTableOutput(ns("coordinates_table")) # 座標テーブル表示
#     )
#   )
# }
# 
# # Server部分
# subsetServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {
#     
#     observeEvent(myReactives$seurat_object,{
#       req(myReactives$seurat_object)
#       update_group_by_for_dimplot(session, myReactives)
#     })
#     
#     # group_byの変更を監視し、unique_groupの選択肢を更新
#     observeEvent(input$group_by, {
#       req(myReactives$seurat_object, input$group_by)
#       update_unique_group_choices(session, myReactives, input$group_by) # unique_groupの選択肢を更新
#     })
#     
#     observeEvent(myReactives$seurat_object, {
#       req(myReactives$seurat_object)
#       
#       # group_by の選択肢を更新
#       update_group_by_for_dimplot(session, myReactives)
#       
#       # 初回の unique_group の選択肢を更新（デフォルトの "sample" を基準に）
#       update_unique_group_choices(session, myReactives, "sample")
#     })
#     
#     
#     # Reactive expression to filter the Seurat object
#     filtered_seurat <- reactive({
#       req(myReactives$seurat_object, input$unique_group, input$group_by)
#       so <- myReactives$seurat_object
#       
#       # フィルタリング対象のセルを取得
#       cells_to_keep <- rownames(so@meta.data)[so@meta.data[[input$group_by]] %in% input$unique_group]
#       
#       # フィルタされたオブジェクトを作成
#       so_filtered <- subset(so, cells = cells_to_keep)
#       
#       # デバッグ用に選択された細胞数を出力
#       print(paste("Filtered cells:", length(cells_to_keep)))
#       
#       return(so_filtered)
#     })
#     
#     observeEvent(input$rerun, {
#       req(filtered_seurat())
#       myReactives <- myseurat_normalize_umap_rerun_subset(myReactives, filtered_seurat())
#     })
#     
#     plot <- reactive({
#       req(filtered_seurat()) # filtered_seurat を正しく取得
#       so <- filtered_seurat()
#       
#       print(dim(so@meta.data))  # デバッグ: メタデータの次元を確認
#       
#       DimPlot(
#         so,
#         reduction = input$reduction,
#         label = TRUE,
#         pt.size = input$point_size,
#         group.by = input$group_by,
#         split.by = NULL,
#         label.size = input$label_size
#       ) +
#         theme(legend.position = input$legend)
#     })
#     
#     output$plot <- renderPlot({
#       print("Rendering plot")  # デバッグ: プロットの更新確認
#       plot()
#     }, width = reactive(input$plot_width), height = reactive(input$plot_height))
#     
#     # 座標テーブルの作成
#     output$coordinates_table <- DT::renderDataTable({
#       req(myReactives$seurat_object)
#       so <- myReactives$seurat_object
#       
#       # reductionに応じた座標を取得
#       reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
#       
#       # metadataと座標を結合
#       table_data <- cbind(so@meta.data, reduction_data)
#       
#       # 表示する列を選択（必要に応じて調整）
#       table_data <- table_data %>% select(input$group_by, starts_with(paste0(toupper(input$reduction), "_")))
#       
#       table_data
#     })
#     
#     
#     # PDFダウンロードハンドラー
#     output$download_plot <- downloadHandler(
#       filename = function() { "UMAP_plot.pdf" },
#       content = function(file) {
#         ggsave(file, plot = renderCustomPlot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
#       }
#     )
#     
#     # CSVダウンロードハンドラー
#     output$download_table <- downloadHandler(
#       filename = function() {
#         paste0("coordinates_", input$reduction, "_", Sys.Date(), ".csv")
#       },
#       content = function(file) {
#         req(myReactives$seurat_object)
#         so <- myReactives$seurat_object
#         
#         # reductionに応じた座標を取得
#         reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
#         
#         # metadataと座標を結合
#         table_data <- cbind(so@meta.data, reduction_data)
#         
#         # 表示する列を選択（必要に応じて調整）
#         table_data <- table_data %>% select(input$group_by, starts_with(paste0(toupper(input$reduction), "_")))
#         
#         write.csv(table_data, file, row.names = FALSE)
#       }
#     )
#     
#   })
# }
# 
# 
# update_group_by_for_dimplot <- function(session, myReactives) {
#   # 除外する列名
#   minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")
#   
#   # メタデータ列を選択
#   metadatas <- myReactives$seurat_object@meta.data %>%
#     select(-all_of(minus_column)) %>%
#     select(-starts_with("TCR"), -starts_with("BCR"))
#   metadata_cols <- names(metadatas)
#   
#   # TCR/BCRデータが存在する場合の追加
#   if (!is.null(myReactives$tcr_df)) {
#     metadata_cols <- append(metadata_cols, c("TCR", "TCR_clonalFrequency", "TCR_cloneSize"))
#   }
#   if (!is.null(myReactives$bcr_df)) {
#     metadata_cols <- append(metadata_cols, c("BCR", "BCR_clonalFrequency", "BCR_cloneSize"))
#   }
#   
#   # グループ列をリスト形式で作成
#   group_cols <- setNames(metadata_cols, metadata_cols)
#   
#   # デフォルトの選択肢を設定
#   default_selection <- if ("sample" %in% metadata_cols) "sample" else metadata_cols[1]
#   # ラジオボタンを更新
#   updateSelectInput(session, "group_by", choices = group_cols, selected = default_selection)
#   updateSelectInput(session, "split_by", choices = group_cols, selected = default_selection)
#   #  updateSelectInput(session, "split_by", choices = c(group_cols, 'none' = NULL), selected = NULL)
# }
# 
# 
# # 
# # # unique_group の選択肢を更新する関数
# # update_unique_group_choices <- function(session, myReactives, group_by_col) {
# #   req(myReactives$seurat_object, group_by_col)
# #   
# #   unique_groups <- unique(myReactives$seurat_object@meta.data[[group_by_col]])
# #   
# #   # 数値順 or 文字列順に並べる
# #   if (all(grepl("^\\d+$", unique_groups))) {
# #     unique_groups <- sort(as.numeric(unique_groups))  # 数値順
# #   } else {
# #     unique_groups <- sort(unique_groups)  # 文字列順
# #   }
# #   
# #   updateCheckboxGroupInput(session, "unique_group", choices = unique_groups, selected = unique_groups) # デフォルトで全て選択
# # }
# 
# myseurat_normalize_umap_rerun_subset <- function(myReactives, seurat_object){
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
#   
# }
