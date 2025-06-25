# UI
geneExpressionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      selectInput(ns("plot_type"), "Plot type", choices = c("Feature plot" = "feature_plot", "Violin plot" = "violin_plot", "Dot plot" = "dot_plot", "Heatmap" = "heatmap"), selected = "feature_plot"),
      # UIのラベルをより分かりやすく変更
      selectizeInput(ns("gene"), 
        label = "Enter or select gene names:", 
        choices = NULL, # 初期値は空にしておき、サーバー側で更新します
        multiple = TRUE, # 複数選択を許可
        options = list(
          placeholder = 'e.g., CD3E, CD19', # 入力欄のヒント
          # 以下の設定で、リストにない遺伝子も手入力できるようになります
          create = TRUE, 
          # 選択した項目を消すボタンを追加
          plugins = list('remove_button') 
        )
      ),
      # ★ Runボタンを追加！
      actionButton(ns("run_plot"), "Run", icon = icon("play")),
      downloadButton(ns("available_feature"), "Download Available Feature List"),
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
      plotOutput(ns('plot'))
    )
  )
}

# Server
geneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {

    # ---- 1. Seuratオブジェクトが読み込まれたときの初期設定 ----
    observeEvent(myReactives$seurat_object, {
      myReactives$available_genes <- c()
      myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])

      # `server = TRUE` は、大量の選択肢でもパフォーマンスを維持するための重要なおまじないです！
      updateSelectizeInput(
        session,
        "gene",
        choices = myReactives$available_genes,
        selected = c('CD3E', 'CD19'), # 初期選択値
        server = TRUE # とっても大事！
      )
      update_group_by_select_input(session, myReactives)
    })

    # ---- 2. Runボタンが押されたら、プロット計算を実行 ----
    plot <- eventReactive(input$run_plot, {
      # 入力値のバリデーション
      req(myReactives$seurat_object, input$gene)
      validate(
        need(length(input$gene) > 0, "Please enter at least one gene name.")
      )

      # Seuratオブジェクトに存在する有効な遺伝子のみを抽出
      valid_genes <- intersect(input$gene, myReactives$available_genes)

      # 有効な遺伝子がない場合はエラーメッセージを表示
      validate(
        need(length(valid_genes) > 0,
             "No valid genes found. Please check the entered gene names or download the available feature list.")
      )

      # プロットタイプに応じてプロットを生成
      tryCatch({
        switch(input$plot_type,
          "feature_plot" = FeaturePlot(myReactives$seurat_object, features = valid_genes, reduction = input$reduction, pt.size = input$point_size) + theme(legend.position = input$legend),
          "violin_plot"  = VlnPlot(myReactives$seurat_object, features = valid_genes, group.by = input$group_by, pt.size = input$point_size) + theme(legend.position = input$legend),
          "dot_plot"     = DotPlot(myReactives$seurat_object, features = valid_genes, group.by = input$group_by) + theme(legend.position = input$legend),
          "heatmap"      = DoHeatmap(myReactives$seurat_object, features = valid_genes, group.by = input$group_by)
        )
      }, error = function(e) {
        validate(
          need(FALSE, paste("An error occurred while generating the plot:", e$message))
        )
      })
    })


    output$plot <- renderPlot({
      plot()
    })

    # ---- 3. ダウンロード処理 ----
    output$download_plot <- downloadHandler(
      filename = function() {
        # isolate() を使って、入力が変更されてもファイル名が再計算されないようにする
        plot_type <- isolate(input$plot_type)
        genes <- isolate(input$gene)
        valid_genes <- intersect(genes, myReactives$available_genes)

        # ファイル名が長くなりすぎないように調整
        gene_string <- if(length(valid_genes) > 3) {
          paste0(paste(head(valid_genes, 3), collapse="_"), "_etc")
        } else {
          paste(valid_genes, collapse="_")
        }

        sprintf("%s_%s.pdf", plot_type, gene_string)
      },
      content = function(file) {
        # plot() は eventReactive なので、これを呼び出すとプロットが再生成される
        p <- plot()
        req(p)

        # isolateで囲むことで、ダウンロード時にwidth/heightが変更されてもエラーにならないようにする
        plot_width <- isolate(input$plot_width %||% 800)
        plot_height <- isolate(input$plot_height %||% 600)

        # ggsaveでPDFとして保存
        ggsave(file, plot = p, width = plot_width / 72, height = plot_height / 72, dpi = 300)
      }
    )

    output$available_feature <- downloadHandler(
      filename = function() {"available_features.txt"},
      content = function(file) {
        req(myReactives$available_genes)
        write.table(sort(myReactives$available_genes), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
    )
  })
}