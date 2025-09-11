#source("../utils.R")
source("utils.R")
# UI
geneExpressionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      selectInput(ns("plot_type"), "Plot type", choices = c("Feature plot" = "feature_plot", "Violin plot" = "violin_plot", "Dot plot" = "dot_plot", "Heatmap" = "heatmap"), selected = "feature_plot"),
      selectizeInput(ns("gene"),
        label = "Enter or select gene names:",
        choices = NULL, # 初期値は空にしておき、サーバー側で更新します
        multiple = TRUE, # 複数選択を許可
        options = list(
          placeholder = 'e.g., CD3E, CD19', # 入力欄のヒント
          create = TRUE,
          plugins = list('remove_button')
        )
      ),
      # ★ Runボタンを削除しました！
      # actionButton(ns("run_plot"), "Run", icon = icon("play")),
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
      # commonPlotOptionsでプロットの幅(plot_width)と高さ(plot_height)のUIが定義されていると仮定します
      commonPlotOptions(ns),
    ),
    mainPanel(
      h3('Plot'),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns('plot'))
    )
  )
}

# # Server
# geneExpressionServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {

#     # ---- 1. Seuratオブジェクトが読み込まれたときの初期設定 ----
#     observeEvent(myReactives$seurat_object, {
#       myReactives$available_genes <- c()
#       myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])

#       # `server = TRUE` は、大量の選択肢でもパフォーマンスを維持するための重要なおまじないです！
#       updateSelectizeInput(
#         session,
#         "gene",
#         choices = myReactives$available_genes,
#         selected = c('CD3E', 'CD19'), # 初期選択値
#         server = TRUE # とっても大事！
#       )
#       update_group_by_select_input(session, myReactives)
#     })

#     # ---- 2. Runボタンが押されたら、プロット計算を実行 ----
#     plot <- eventReactive(input$run_plot, {
#       # 入力値のバリデーション
#       req(myReactives$seurat_object, input$gene)
#       validate(
#         need(length(input$gene) > 0, "Please enter at least one gene name.")
#       )

#       # Seuratオブジェクトに存在する有効な遺伝子のみを抽出
#       valid_genes <- intersect(input$gene, myReactives$available_genes)

#       # 有効な遺伝子がない場合はエラーメッセージを表示
#       validate(
#         need(length(valid_genes) > 0,
#              "No valid genes found. Please check the entered gene names or download the available feature list.")
#       )

#       # プロットタイプに応じてプロットを生成
#       tryCatch({
#         switch(input$plot_type,
#           "feature_plot" = FeaturePlot(myReactives$seurat_object, features = valid_genes, reduction = input$reduction, pt.size = input$point_size) + theme(legend.position = input$legend),
#           "violin_plot"  = VlnPlot(myReactives$seurat_object, features = valid_genes, group.by = input$group_by, pt.size = input$point_size) + theme(legend.position = input$legend),
#           "dot_plot"     = DotPlot(myReactives$seurat_object, features = valid_genes, group.by = input$group_by) + theme(legend.position = input$legend),
#           "heatmap"      = DoHeatmap(myReactives$seurat_object, features = valid_genes, group.by = input$group_by)
#         )
#       }, error = function(e) {
#         validate(
#           need(FALSE, paste("An error occurred while generating the plot:", e$message))
#         )
#       })
#     })


#     output$plot <- renderPlot({
#       plot()
#     })

#     # ---- 3. ダウンロード処理 ----
#     output$download_plot <- downloadHandler(
#       filename = function() {
#         # isolate() を使って、入力が変更されてもファイル名が再計算されないようにする
#         plot_type <- isolate(input$plot_type)
#         genes <- isolate(input$gene)
#         valid_genes <- intersect(genes, myReactives$available_genes)

#         # ファイル名が長くなりすぎないように調整
#         gene_string <- if(length(valid_genes) > 3) {
#           paste0(paste(head(valid_genes, 3), collapse="_"), "_etc")
#         } else {
#           paste(valid_genes, collapse="_")
#         }

#         sprintf("%s_%s.pdf", plot_type, gene_string)
#       },
#       content = function(file) {
#         # plot() は eventReactive なので、これを呼び出すとプロットが再生成される
#         p <- plot()
#         req(p)

#         # isolateで囲むことで、ダウンロード時にwidth/heightが変更されてもエラーにならないようにする
#         plot_width <- isolate(input$plot_width %||% 800)
#         plot_height <- isolate(input$plot_height %||% 600)

#         # ggsaveでPDFとして保存
#         ggsave(file, plot = p, width = plot_width / 72, height = plot_height / 72, dpi = 300)
#       }
#     )

#     output$available_feature <- downloadHandler(
#       filename = function() {"available_features.txt"},
#       content = function(file) {
#         req(myReactives$available_genes)
#         write.table(sort(myReactives$available_genes), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
#       }
#     )
#   })
# }

# Server
# Server
geneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {

    # ---- 1. Seuratオブジェクトが読み込まれたときの初期設定 ----
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])

      updateSelectizeInput(
        session,
        "gene",
        choices = myReactives$available_genes,
        selected = c('CD3E', 'CD19'),
        server = TRUE
      )
      update_group_by_select_input(session, myReactives)
      update_reduction_choices(session, myReactives)
    })

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    # ---- 2. Runボタンなしで、プロットをリアクティブに生成 ----
    plot <- reactive({
      # 入力値のバリデーション
      req(myReactives$seurat_object, input$gene, myReactives$available_genes)
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
      
      # プロットタイプに応じて、必要な入力値が揃うまで待機
      req(
        # ★ 修正: !is.null() を外して直接評価するように変更
        if (input$plot_type != "feature_plot") input$group_by else TRUE,
        if (input$plot_type == "feature_plot") input$reduction else TRUE
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

    # ---- 3. プロットのレンダリング（サイズの動的変更に対応） ----
    output$plot <- renderPlot({
      req(plot())
      plot()
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))

    # ---- 4. ダウンロード処理（動的なサイズで保存） ----
    output$download_plot <- downloadHandler(
      filename = function() {
        plot_type <- isolate(input$plot_type)
        genes <- isolate(input$gene)
        
        req(myReactives$available_genes)
        valid_genes <- intersect(genes, myReactives$available_genes)

        gene_string <- if(length(valid_genes) > 3) {
          paste0(paste(head(valid_genes, 3), collapse="_"), "_etc")
        } else {
          paste(valid_genes, collapse="_")
        }
        if(gene_string == "") gene_string <- "no_valid_genes"

        sprintf("%s_%s.pdf", plot_type, gene_string)
      },
      content = function(file) {
        p <- plot()
        req(p)

        # UIで指定された幅と高さを取得
        plot_width <- req(input$plot_width)
        plot_height <- req(input$plot_height)

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


# reductionPlotServerから持ってきたヘルパー関数
# Seuratオブジェクトに存在するリダクション名を取得し、UIの選択肢を動的に更新する
update_reduction_choices <- function(session, myReactives) {
  ns <- session$ns
  req(myReactives$seurat_object)

  # Seuratオブジェクトからリダクション名（'pca', 'umap', 'tsne'など）を取得
  reduction_names <- names(myReactives$seurat_object@reductions)

  # UIで表示する名前（例: 'umap' -> 'UMAP'）と、サーバー側で使う値（'umap'）のペアを作成
  choices <- stats::setNames(reduction_names, toupper(reduction_names))

  # デフォルトの選択肢を決定（'umap'があればそれを、なければ最初のものを選択）
  default_selection <- if ("umap" %in% reduction_names) "umap" else reduction_names[1]

  # selectInputを更新
  updateSelectInput(session, "reduction", choices = choices, selected = default_selection)
}