# UI
geneExpressionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      selectInput(ns("plot_type"), "Plot type", choices = c("Feature plot" = "feature_plot", "Violin plot" = "violin_plot", "Dot plot" = "dot_plot", "Heatmap" = "heatmap"), selected = "feature_plot"),
      textInput(ns("gene"), "Enter feature (gene) names (ex. CD3E, CD19, CD14):", value = "CD3E, CD19"),
      downloadButton(ns("available_feature"), "Download Available Feature Name"),
      conditionalPanel(
        condition = sprintf("input['%s'] != 'feature_plot'", ns("plot_type")),
        groupByInput(ns),
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'feature_plot'", ns("plot_type")),
        reductionInput(ns),
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'feature_plot' || input['%s'] == 'violin_plot'", ns("plot_type"), ns("plot_type")),
        pointSizeInput(ns),
      ),
      commonPlotOptions(ns),
    ),
    mainPanel(
      h3('Plot'),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot"))
    )
  )
}

geneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)

      # Seuratオブジェクトがロードされたら、利用可能な遺伝子リストを更新し、初期値を検証する
      myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])

      # 初期値を検証し、myReactives$valid_genesを設定
      initial_genes <- unlist(strsplit(input$gene, ",\\s*"))
      myReactives$valid_genes <- trimws(initial_genes)[trimws(initial_genes) %in% myReactives$available_genes]
    })

    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      # 上記のobserveEventで処理するため、ここでは不要
      # myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])
    })

    observeEvent(input$gene, {
      req(input$gene)
      req(myReactives$available_genes)
      text_list <- unlist(strsplit(input$gene, ",\\s*"))
      myReactives$valid_genes <- trimws(text_list)[trimws(text_list) %in% myReactives$available_genes]
    })

    # プロットを生成する関数（他のプロットにも再利用可能）
    plot <- reactive({
      req(myReactives$seurat_object)
      req(myReactives$valid_genes) # valid_genesが存在することを確認

      if (input$plot_type == "feature_plot") {
        FeaturePlot(
          myReactives$seurat_object,
          reduction = input$reduction,
          features = myReactives$valid_genes,
          pt.size = input$point_size,
        ) + theme(legend.position = input$legend)
      } else if (input$plot_type == "violin_plot") {
        VlnPlot(
          myReactives$seurat_object,
          features = myReactives$valid_genes,
          pt.size = input$point_size,
          group.by = input$group_by
        ) + theme(legend.position = input$legend)
      } else if (input$plot_type == "dot_plot") {
        DotPlot(
          myReactives$seurat_object,
          features = myReactives$valid_genes,
          group.by = input$group_by
        ) + theme(legend.position = input$legend)
      } else if (input$plot_type == "heatmap") {
        DoHeatmap(
          myReactives$seurat_object,
          features = myReactives$valid_genes,
          group.by = input$group_by
        )
      }
    })

    # プロットのレンダリング
    output$plot <- renderPlot(
      {
        plot()
      },
      width = reactive(input$plot_width),
      height = reactive(input$plot_height)
    )

    # PDFダウンロードハンドラー
    output$download_plot <- downloadHandler(
      filename = function() {
        "gene_expression_plot.pdf"
      }, # ファイル名を修正
      content = function(file) {
        ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300) # renderCustomPlot() を plot() に修正
      }
    )

    output$available_feature <- downloadHandler(
      filename = function() {
        "available_features.txt"
      },
      content = function(file) {
        # 遺伝子名をソートしてから書き出す
        sorted_genes <- sort(myReactives$available_genes)
        write.table(sorted_genes, file, row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
    )
  })
}
