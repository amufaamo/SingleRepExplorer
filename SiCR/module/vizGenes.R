# UI部分
vizGenesUI <- function(id, vdj = 'tcr') {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns('group_by'), 
                   label = 'Group', 
                   choices = c('sample', 'seurat_clusters'),
                   selected = 'sample'),
      radioButtons(ns('plot'), "Plot", choices = list('heatmap', 'barplot'), selected = 'heatmap'),
      radioButtons(ns('order'), "Order", choices = list('gene', 'variance'), selected = 'gene'),
      radioButtons(ns('scale'), label = "Scale",
                   choices = c('Relative' = TRUE, "Absolute" = FALSE),
                   selected = TRUE),
      if (vdj == 'tcr'){
        radioButtons(ns('x_axis'), label = "X axis", 
                     choices = list("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"),
                     selected = 'TRBV')
      } else if (vdj == 'bcr'){
        radioButtons(ns('x_axis'), label = "X axis", 
                     choices = list("IGHV", "IGHD", "IGHJ", "IGLV", "IGLJ", "IGKV", "IGKJ"),
                     selected = 'IGHV')
      },
      if (vdj == 'tcr'){
        radioButtons(ns('y_axis'), label = "Y axis", 
                     choices = list("none", "TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"),
                     selected = "none")
      } else if (vdj == 'bcr'){
        radioButtons(ns('y_axis'), label = "Y axis", 
                     choices = list("none", "IGHV", "IGHD", "IGHJ", "IGLV", "IGLJ", "IGKV", "IGKJ"),
                     selected = "none")
      },
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

vizGenesServer <- function(id, myReactives, vdj = 'tcr') {
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
      y_axis_value <- if (input$y_axis == "none") NULL else input$y_axis
      vizGenes(df,
               group.by = input$group_by,
               plot = input$plot,
               order = input$order,
#               scale = input$scale,
               x.axis = input$x_axis,
               y.axis = y_axis_value
              ) + theme(legend.position = input$legend)
    })
    
    table <- reactive({
      if (vdj == 'tcr'){
        df <- myReactives$tcr_df
      } else if (vdj == 'bcr'){
        df <- myReactives$bcr_df
      }
      req(df)
      
      vizGenes(df,
               group.by = input$group_by,
               plot = input$plot,
               order = input$order,
#               scale = input$scale,
               x.axis = input$x_axis,
               y.axis = input$y_axis,
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
