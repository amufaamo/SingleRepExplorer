# UI部分
clonalAbundanceUI <- function(id, vdj = 'tcr') {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns('group_by'), 
                   label = 'Group', 
                   choices = c('sample', 'seurat_clusters'),
                   selected = 'sample'),
      radioButtons(ns('scale'), label = "Scale",
                   choices = c('Relative' = TRUE, "Absolute" = FALSE),
                   selected = TRUE),
      if (vdj == 'tcr'){
        radioButtons(ns('chain'), label = "Chain", 
                     choices = list("both", "TRA", "TRB"),
                     selected = 'both')
      } else if (vdj == 'bcr'){
        radioButtons(ns('chain'), label = "Chain", 
                     choices = list("both", "IGH", "IGL"),
                     selected = 'both')
      },
      radioButtons(ns('clone_call'),
                   label = "Clone call",
                   choices = list("Clonotype ID" = "raw_clonotype_id",
                                  "Sub clonotype ID" = "exact_subclonotype_id",
                                  "VDJC + CDR3 nucleotide" = "strict",
                                  "VDJC" = "gene",
                                  "CDR3 nucleotide" = "nt",
                                  "CDR3 amino acid" = "aa"),
                   selected = "raw_clonotype_id"),
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

clonalAbundanceServer <- function(id, myReactives, vdj = 'tcr') {
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
      sample_var <- input$group_by
      df <- df$S1 %>% group_by(!!sym(sample_var)) %>%
        group_split() %>%
        setNames(unique(df$S1[[sample_var]]))
      
      
      clonalAbundance(df,
                  group.by = input$group_by,
                  scale = input$scale,
                  cloneCall = input$clone_call,
                  chain = input$chain,
                  palette = 'Dynamic'
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
      
      sample_var <- input$group_by
      df <- df$S1 %>% group_by(!!sym(sample_var)) %>%
        group_split() %>%
        setNames(unique(df$S1[[sample_var]]))
      
      clonalAbundance(df,
                  group.by = input$group_by,
                  scale = input$scale,
                  cloneCall = input$clone_call,
                  chain = input$chain,
                  palette = 'Dynamic',
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
