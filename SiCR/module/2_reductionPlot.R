# UI部分
reductionPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      reductionInput(ns),
      groupByInput(ns),

#      selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
      checkboxInput(ns("highlight_toggle"), "Highlight Cells", value = FALSE), # Toggle switch
      conditionalPanel( # Conditional panel for highlight options
        condition = sprintf("input['%s'] == true", ns("highlight_toggle")),  # Corrected condition
        checkboxGroupInput(ns("unique_group"), "Select Group", choices = NULL, inline = TRUE)
      ),
      pointSizeInput(ns),
      labelSizeInput(ns),
#      sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 10, value = 0.1, step = 0.01),
#      sliderInput(ns("label_size"), "Size of labels", min = 0, max = 20, value = 10, step = 1),
      commonPlotOptions(ns),
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      hr(),
      h3("Coordinates Table"),
      downloadButton(ns("download_table"), "Download Table (.csv)"), # ダウンロードボタン
      DT::dataTableOutput(ns("coordinates_table")) # 座標テーブル表示
    )
  )
}

# Server部分
reductionPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(myReactives$seurat_object,{
      req(myReactives$seurat_object)
      # Seuratオブジェクトから利用可能なリダクションを動的に読み込みUIを更新
      update_reduction_choices(session, myReactives)
      update_group_by_for_dimplot(session, myReactives)
      # Initial update of unique_group choices based on default "sample"
      update_unique_group_choices(session, myReactives, "sample")
    })
    
    observeEvent(input$group_by, {
      req(myReactives$seurat_object, input$group_by)
      update_unique_group_choices(session, myReactives, input$group_by)
      # Reset highlight_toggle to FALSE when group_by changes
      updateCheckboxInput(session, "highlight_toggle", value = FALSE)
    })
    
    
    # プロットを生成する関数（他のプロットにも再利用可能）
    plot <- reactive({
      req(myReactives$seurat_object, input$reduction) # FIX: Ensure input$reduction is not empty
      so <- myReactives$seurat_object
      so@meta.data$orig.ident <- NULL
      
      if (input$highlight_toggle) {
        # Highlight Plot Logic
        Idents(so) <- input$group_by
        highlights <- WhichCells(so, idents = input$unique_group)
        
        p <- DimPlot(
          so,
          reduction = input$reduction, # UIから渡される値は既に正しい小文字の名前です
          group.by = input$group_by,
          cells.highlight = highlights,
          cols.highlight = 'darkblue',
          cols = 'grey',
          pt.size = input$point_size, # Add point size
          label.size = input$label_size, #add label size
          label = TRUE
          
        ) + theme(legend.position = input$legend)
        
      } else {
        # Standard DimPlot Logic
        p <- DimPlot(
          so,
          reduction = input$reduction, # UIから渡される値は既に正しい小文字の名前です
          label = TRUE,
          pt.size = input$point_size,
          group.by = input$group_by,
          split.by = NULL,
          label.size = input$label_size
        ) +
          theme(legend.position = input$legend)
      }
      return(p)
    })
    
    # プロットのレンダリング
    output$plot <- renderPlot({
      plot()
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))
    
    # 座標テーブルの作成
    output$coordinates_table <- DT::renderDataTable({
      req(myReactives$seurat_object, input$reduction) # FIX: Ensure input$reduction is not empty
      so <- myReactives$seurat_object
      
      # reductionに応じた座標を取得
      reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
      
      # metadataと座標を結合
      table_data <- cbind(so@meta.data, reduction_data)
      
      # 表示する列を選択（必要に応じて調整）
      table_data <- table_data %>% select(input$group_by, starts_with(paste0(input$reduction, "_")))
      
      table_data
    })
    
    
    # PDFダウンロードハンドラー
    output$download_plot <- downloadHandler(
      filename = function() { "UMAP_plot.pdf" },
      content = function(file) {
        ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
      }
    )
    
    # CSVダウンロードハンドラー
    output$download_table <- downloadHandler(
      filename = function() {
        paste0("coordinates_", input$reduction, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(myReactives$seurat_object, input$reduction) # FIX: Ensure input$reduction is not empty
        so <- myReactives$seurat_object
        
        # reductionに応じた座標を取得
        reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
        
        # metadataと座標を結合
        table_data <- cbind(so@meta.data, reduction_data)
        
        # 表示する列を選択（必要に応じて調整）
        table_data <- table_data %>% select(input$group_by, starts_with(paste0(input$reduction, "_")))
        
        write.csv(table_data, file, row.names = FALSE)
      }
    )
    
  })
}

# Seuratオブジェクトに存在するリダクション名を取得し、UIの選択肢を動的に更新する関数
update_reduction_choices <- function(session, myReactives) {
  ns <- session$ns
  req(myReactives$seurat_object)
  
  # Seuratオブジェクトからリダクション名（'pca', 'umap', 'tsne'など）を取得
  reduction_names <- names(myReactives$seurat_object@reductions)
  
  # UIで表示する名前（例: 'umap' -> 'UMAP'）と、サーバー側で使う値（'umap'）のペアを作成
  # toupperは 'tsne' を 'TSNE' に変換します。'T-SNE' のようにしたい場合は個別の処理が必要です。
  choices <- stats::setNames(reduction_names, toupper(reduction_names))
  
  # デフォルトの選択肢を決定（'umap'があればそれを、なければ最初のものを選択）
  default_selection <- if ("umap" %in% reduction_names) "umap" else reduction_names[1]
  
  # selectInputを更新
  updateSelectInput(session, "reduction", choices = choices, selected = default_selection)
}