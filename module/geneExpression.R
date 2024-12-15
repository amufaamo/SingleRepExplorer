# UI部分
geneExpressionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns('plot_type'), "Plot type", choices = c('Feature plot' = 'feature_plot', "Violin plot" = 'violin_plot', "Dot plot" = "dot_plot", "Heatmap" = "heatmap"), selected = "feature_plot"),
      textInput(ns("gene"), "Enter feature (gene) names (ex. CD3E, CD19, CD14):", value = "CD3E, CD19"),
      downloadButton(ns('available_feature'), "You can download available gene name"),
      radioButtons(ns('gene'), "Plot type", choices = c('Feature plot' = 'feature_plot', "Violin plot" = 'violin_plot', "Dot plot" = "dot_plot", "Heatmap" = "heatmap"), selected = "feature_plot"),
      radioButtons(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
      radioButtons(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne", "PCA" = "pca"), selected = "umap"),
      sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 10, value = 0.1, step = 0.01),
      sliderInput(ns("label_size"), "Size of labels", min = 0, max = 20, value = 10, step = 1),
      radioButtons(ns("legend"), "Legend", choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
      sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
      sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
      downloadButton(ns("download_plot"), "Download plot (.pdf)")
    ),
    mainPanel(
      plotOutput(ns("plot"))
    )
  )
}

# Server部分
geneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(myReactives$seurat_object,{
      req(myReactives$seurat_object)
      update_group_by(session, myReactives)
    })
    
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])
    })
    
    observeEvent(input$gene, {
      req(input$gene)
      req(myReactives$available_genes)
      text_list <- unlist(strsplit(input$gene, ",\\s*"))
      myReactives$valid_genes <- text_list[text_list %in% myReactives$available_genes]
    })
    
    # プロットを生成する関数（他のプロットにも再利用可能）
    plot <- reactive({
      req(myReactives$seurat_object)
      
      if (input$plot_type == 'feature_plot'){
        FeaturePlot(
          myReactives$seurat_object,
          reduction = input$reduction,
          features = myReactives$valid_genes,
          pt.size = input$point_size,
        ) + theme(legend.position = input$legend)
      } else if (input$plot_type == 'violin_plot'){
        VlnPlot(
          myReactives$seurat_object,
          features = myReactives$valid_genes,
          pt.size = input$point_size,
          group.by = input$group_by
        ) + theme(legend.position = input$legend)
      } else if (input$plot_type == 'dot_plot'){
        DotPlot(
          myReactives$seurat_object,
          features = myReactives$valid_genes,
          group.by = input$group_by
        ) + theme(legend.position = input$legend)
      } else if (input$plot_type == 'heatmap'){
        DoHeatmap(
          myReactives$seurat_object,
          features = myReactives$valid_genes,
          group.by = input$group_by
        )
      }
      
      
    })
    
    # プロットのレンダリング
    output$plot <- renderPlot({
      plot()
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))
    
    # PDFダウンロードハンドラー
    output$download_plot <- downloadHandler(
      filename = function() { "UMAP_plot.pdf" },
      content = function(file) {
        ggsave(file, plot = renderCustomPlot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
      }
    )
  })
}

update_group_by <- function(session, myReactives) {
  print('update_group_by')
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")
  metadatas <- myReactives$seurat_object@meta.data %>%
    select(-minus_column) %>%
    select(!starts_with("TCR")) %>%
    select(!starts_with("BCR"))
  metadata_cols <- names(metadatas)
  
  group_cols <- list()
  for (col in metadata_cols) {
    group_cols[[col]] <- col
  }
  updateRadioButtons(session, "group_by", choices = group_cols, selected = "sample")
}

update_group_by_for_dimplot <- function(session, myReactives) {
  print('update_group_by')
  
  # 除外する列名
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")
  
  # メタデータ列を選択
  metadatas <- myReactives$seurat_object@meta.data %>%
    select(-all_of(minus_column)) %>%
    select(-starts_with("TCR"), -starts_with("BCR"))
  metadata_cols <- names(metadatas)
  
  # TCR/BCRデータが存在する場合の追加
  if (!is.null(myReactives$tcr_df)) {
    metadata_cols <- append(metadata_cols, c("TCR", "TCR_clonalFrequency", "TCR_cloneSize"))
  }
  if (!is.null(myReactives$bcr_df)) {
    metadata_cols <- append(metadata_cols, c("BCR", "BCR_clonalFrequency", "BCR_cloneSize"))
  }
  
  # グループ列をリスト形式で作成
  group_cols <- setNames(metadata_cols, metadata_cols)
  
  # デフォルトの選択肢を設定
  default_selection <- if ("sample" %in% metadata_cols) "sample" else metadata_cols[1]
  
  # ラジオボタンを更新
  updateRadioButtons(session, "group_by", choices = group_cols, selected = default_selection)
}
