# UI部分
clonalProportionUI <- function(id, vdj = 'tcr') {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns('group_by'), 
                   label = 'Group', 
                   choices = c('sample', 'seurat_clusters'),
                   selected = 'sample'),
                   radioButtons(ns('clone_call'),
                                label = "Clone call",
                                choices = list("Clonotype ID" = "raw_clonotype_id",
                                               "Sub clonotype ID" = "exact_subclonotype_id",
                                               "VDJC + CDR3 nucleotide" = "strict",
                                               "VDJC" = "gene",
                                               "CDR3 nucleotide" = "nt",
                                               "CDR3 amino acid" = "aa"),
                                selected = "raw_clonotype_id"),
      if (vdj == 'tcr'){
        radioButtons(ns('chain'), label = "Chain", 
                     choices = list("both", "TRA", "TRB"),
                     selected = 'both')
      } else if (vdj == 'bcr'){
        radioButtons(ns('chain'), label = "Chain", 
                     choices = list("both", "IGH", "IGL"),
                     selected = 'both')
      },
      numericInput(ns("unique"), label = "Unique clones (default: 1)", value = 1),
      numericInput(ns("rare"), label = "Rare clones (default: 5)", value = 5),
      numericInput(ns("uncommon"), label = "Uncommon clones (default: 10)", value = 10),
      numericInput(ns("intermediate"), label = "Intermediate clones (default: 100)", value = 100),
      numericInput(ns("expanded"), label = "Expanded clones (default: 1000)", value = 1000),
      numericInput(ns("dominant"), label = "Dominant clones (default: 10000)", value = 10000),
      radioButtons(ns("legend"), "Legend", choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
      sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
      sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      downloadButton(ns("download_table"), "Download table (.csv)") # 追加
    ),
    mainPanel(
      plotOutput(ns("plot")),
      DTOutput(ns('table'))
    )
  )
}

# Server部分

clonalProportionServer <- function(id, myReactives, vdj = 'tcr') {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(myReactives$seurat_object,{
      req(myReactives$seurat_object)
      update_group_by(session, myReactives)
    })
    
    
    plot <- reactive({
      if (vdj == 'tcr'){
        df <- myReactives$tcr_df
      } else if (vdj == 'bcr'){
        df <- myReactives$bcr_df
      }
      req(df)
      

      clonalProportion(df,
                  group.by = input$group_by,
                  cloneCall = input$clone_call,
                  chain = input$chain,
                  clonalSplit = c(input$unique,
                                  input$rare,
                                  input$uncommon,
                                  input$intermediate,
                                  input$expanded,
                                  input$dominant)
                  ) +
        scale_y_continuous(expand = c( 0, 0 )) +
        theme(legend.position = input$legend)
    })
    
    table <- reactive({
      if (vdj == 'tcr'){
        df <- myReactives$tcr_df
      } else if (vdj == 'bcr'){
        df <- myReactives$bcr_df
      }
      req(df)
      
      
      clonalProportion(df,
                       group.by = input$group_by,
                  cloneCall = input$clone_call,
                  chain = input$chain,
                  clonalSplit = c(input$unique,
                                  input$rare,
                                  input$uncommon,
                                  input$intermediate,
                                  input$expanded,
                                  input$dominant),
                  exportTable = TRUE
      )
    })
    
        
    output$plot <- renderPlot({
      plot()
    }, width = reactive(input$plot_width), height = reactive(input$plot_height)
    )
    
    output$table <- renderDT({
      table()
    })
    
    # PDFダウンロードハンドラー
    output$download_plot <- downloadHandler(
      filename = function() { "plot.pdf" },
      content = function(file) {
        ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
      }
    )
    
    
    # CSVダウンロードハンドラー
    output$download_table <- downloadHandler(
      filename = function() { "table.csv" },
      content = function(file) {
        write.csv(table(), file)  # table()で生成したデータをCSVとして保存
      }
    )
    

  })
}
