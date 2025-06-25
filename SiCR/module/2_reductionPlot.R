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
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      so@meta.data$orig.ident <- NULL
      
      if (input$highlight_toggle) {
        # Highlight Plot Logic
        Idents(so) <- input$group_by
        highlights <- WhichCells(so, idents = input$unique_group)
        
        p <- DimPlot(
          so,
          reduction = input$reduction,
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
          reduction = input$reduction,
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
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      
      # reductionに応じた座標を取得
      reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
      
      # metadataと座標を結合
      table_data <- cbind(so@meta.data, reduction_data)
      
      # 表示する列を選択（必要に応じて調整）
      table_data <- table_data %>% select(input$group_by, starts_with(paste0(toupper(input$reduction), "_")))
      
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
        req(myReactives$seurat_object)
        so <- myReactives$seurat_object
        
        # reductionに応じた座標を取得
        reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
        
        # metadataと座標を結合
        table_data <- cbind(so@meta.data, reduction_data)
        
        # 表示する列を選択（必要に応じて調整）
        table_data <- table_data %>% select(input$group_by, starts_with(paste0(toupper(input$reduction), "_")))
        
        write.csv(table_data, file, row.names = FALSE)
      }
    )
    
  })
}


update_group_by_for_dimplot <- function(session, myReactives) {
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
  updateSelectInput(session, "group_by", choices = group_cols, selected = default_selection)
  # updateSelectInput(session, "split_by", choices = group_cols, selected = default_selection) # Remove split_by
  #  updateSelectInput(session, "split_by", choices = c(group_cols, 'none' = NULL), selected = NULL)
  
  
}

# Function to update choices for unique_group
update_unique_group_choices <- function(session, myReactives, group_by_value) {
  req(myReactives$seurat_object)
  unique_groups <- unique(myReactives$seurat_object@meta.data[[group_by_value]])
  updateCheckboxGroupInput(session, "unique_group", choices = unique_groups, selected = NULL, inline = TRUE)
}