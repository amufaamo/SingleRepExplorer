# UI部分
highlightPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns("reduction"), "Reduction", choices = c("UMAP" = "umap", "T-SNE" = "tsne", "PCA" = "pca"), selected = "umap"),
      radioButtons(ns('vdj'), "VDJ", choices = c("TCR", "BCR"), selected = "TCR"),
      radioButtons(ns('clone_call'),
                   label = "Clone call",
                   choices = list("Clonotype ID" = "raw_clonotype_id",
                                  "Sub clonotype ID" = "exact_subclonotype_id",
                                  "VDJC + CDR3 nucleotide" = "CTstrict",
                                  "VDJC" = "CTgene",
                                  "CDR3 nucleotide" = "CTnt",
                                  "CDR3 amino acid" = "CTaa"),
                   selected = "raw_clonotype_id"),
#      checkboxGroupInput(ns('selected_values'), label = "Value", choices = NULL, selected = NULL),
      selectInput(ns('selected_values'), label = "Value", choices = NULL, selected = NULL, multiple = TRUE),
#      uiOutput(ns('selectUI')),
      sliderInput(ns("highlight_point_size"), "Size of highlight points", min = 0.01, max = 10, value = 0.1, step = 0.01),
      sliderInput(ns("un_highlight_point_size"), "Size of NA points", min = 0.01, max = 10, value = 0.1, step = 0.01),
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

highlightPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    
    selected_column <- reactive({
      req(myReactives$seurat_object)
      # selected_columnを作成
      selected_column <- paste0(input$vdj, "_", input$clone_call)
      return(selected_column)
    })
    
    observeEvent(selected_column(), {
      req(myReactives$seurat_object)
      new_choices <- unique(myReactives$seurat_object@meta.data[[selected_column()]])
      new_choices <- str_sort(new_choices, numeric = TRUE)
#      new_choices <- sort(unique(myReactives$seurat_object@meta.data[[selected_column()]]))
      updateSelectInput(session, "selected_values", choices = new_choices, selected = character(0))
#      updateCheckboxGroupInput(session, "selected_values", choices = new_choices, selected = character(0))
    })
    

    observe({
      req(input$selected_values)
      print(input$selected_values)
    })
    
#    DimPlot(integrated, label=T, group.by="Treat", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")
                                                    
    
    plot <- reactive({
      req(myReactives$seurat_object)
      so <- highlightClones(myReactives$seurat_object,
                             cloneCall = selected_column(),
                             sequence = input$selected_values)
      # DimPlotの作成
      DimPlot(
        so,
        reduction = input$reduction,
        group.by = 'highlight',
        pt.size = ifelse(is.na(so@meta.data$'highlight'), input$un_highlight_point_size, input$highlight_point_size)
#        pt.size = input$point_size
      ) +
        theme(legend.position = input$legend)
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

