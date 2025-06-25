# UI部分
clonalCountUI <- function(id, vdj = 'tcr') {
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
      numericInput(ns("rare"), label = "Rare clones (default: 2)", value = 2),
      numericInput(ns("uncommon"), label = "Uncommon clones (default: 5)", value = 5),
#      numericInput(ns("intermediate"), label = "Intermediate clones (default: 10)", value = 10),
      numericInput(ns("expanded"), label = "Expanded clones (default: 20)", value = 20),
      numericInput(ns("dominant"), label = "Dominant clones (default: 100)", value = 100),
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

clonalCountServer <- function(id, myReactives, vdj = 'tcr') {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(myReactives$seurat_object,{
      req(myReactives$seurat_object)
      update_group_by(session, myReactives)
    })
    
    table <- reactive({
#      req(myReactives$tcr_path)
#      req(myReactives$bcr_path)
      imm <- NULL
      imm <- immunarch_data(myReactives, input, vdj)
      repClonality(imm$data, 
                   .method = 'rare',
                   .bound = c(input$unique,
                              input$rare,
                              input$uncommon,
                              input$expanded) #,
#                              input$dominant)
      ) 
    })
    
    plot <- reactive({
      req(table())
      vis(table()) + scale_y_continuous(expand = c(0, 0))
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



immunarch_data <- function(myReactives, input, vdj){
  unlink('immunarch', recursive = TRUE)
  dir.create("immunarch", showWarnings = FALSE)
  if(vdj == 'tcr'){
    req(myReactives$tcr_path)
    dt <- read.csv(myReactives$tcr_path)
  } else if (vdj == 'bcr'){
    req(myReactives$bcr_path)
    dt <- read.csv(myReactives$bcr_path)
  }
  dt <- dt %>% mutate(sample = str_remove(barcode, "^.+-"))
  split_column <- input$group_by
  data_list <- split(dt, dt[[split_column]])
  for (name in names(data_list)) {
    write.csv(data_list[[name]], file = paste0("immunarch/", name, ".csv"), row.names = FALSE, quote = FALSE)
  }
  imm <- repLoad("immunarch/")
  unlink('immunarch')
  return(imm)
}
