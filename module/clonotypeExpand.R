
clonotypeExpandUI <- function(id){
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns("heat_or_bar"), "Heatmap or barplot?", choices = list("heatmap" = "heatmap", "barplot" = "barplot"), selected = "barplot"),
      radioButtons(ns("count_or_percent"), "Use cell count or percent?", choices = list("count" = "count", "percent" = "percent"), selected = "count"),
      selectInput(ns("group_by"), "Group by", choices = list("sample" = "sample")),
      selectInput(ns("focus_group"), "Show top 20 in this group", choices = ""),
      selectInput(ns("palette"), "Color palette", choices = ""),
      checkboxInput(ns("legend"), "Show legend", value = TRUE),
      sliderInput(ns("plot_width"),  "Width",  min = 100, max = 2000, value = 500, step = 100),
      sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
    ),
    mainPanel(
      plotOutput(ns("clonotype_expand_plot")),
      downloadButton(ns("download_data"), "Download data (.csv)"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)")
    )
  )
}


clonotypeExpandServer <- function(id, data, group_cols) {
  moduleServer(id, function(input, output, session) {
    
    plot_width <- reactive(input$plot_width)
    plot_height <- reactive(input$plot_height)
    observe(updateSelectInput(session, "group_by", choices = group_cols))
    observe(updateSelectInput(session, "focus_group", choices = unique(data[input$group_by])))
    observeEvent(input$heat_or_bar, {
      if (input$heat_or_bar == "barplot") {
        updateSelectInput(session, "palette", choices = palette_list)
      } else if (input$heat_or_bar == "heatmap") {
        updateSelectInput(session, "palette", choices = list("Greys"="grey40", "Blues"="#2171B5", "Oranges"="#F16913", "Purples"="#807DBA"))
      }
    })
    
    # tally
    clonotypeTally <- reactive({
      my_tally(
        data,
        drop_na = TRUE,
        arrange = "desc",
        group_by = input$group_by,
        x = "raw_clonotype_id"
      )
    })
    
    # Extract top20
    top <- 20
    clonotypeTallyTop <- reactive({
      top_clonotype_ids <- clonotypeTally() %>%
        filter(if_all(all_of(input$group_by), ~ . == input$focus_group)) %>%
        top_n(top, count) %>% # this will result in more then 20 rows if there are ties
        head(top) %>%
        pull(raw_clonotype_id)
      clonotypeTally() %>%
        filter(raw_clonotype_id %in% top_clonotype_ids) %>%
        mutate(raw_clonotype_id = factor(raw_clonotype_id, levels = top_clonotype_ids)) # To order from top1 to 20 in the graph
    })
    
    # Get plot
    clonotypeExpandPlot <- reactive({
      
      # heatmap
      if(input$heat_or_bar == "heatmap") {
          if(input$count_or_percent == "count") {
            round_n <- 0
          } else {
            round_n <- 1
          }
        g <- my_heatmap(
          clonotypeTallyTop(),
          x = "raw_clonotype_id",
          y = input$group_by,
          fill = input$count_or_percent,
          round = round_n,
          legend_position = input$legend,
          high_color = input$palette,
          color_limit = NA
        )
      }
      
      # barplot
      if(input$heat_or_bar == "barplot") {
        g <- my_barplot(
          clonotypeTallyTop(),
          x = "raw_clonotype_id",
          y = input$count_or_percent,
          fill = input$group_by,
          position = "dodge",
          legend_position = input$legend,
          flip = FALSE,
          palette = input$palette
        )
      }
      
      return(g)
      
    })
    
    # output
    output$clonotype_expand_plot <- renderPlot(
      clonotypeExpandPlot(),
      width  = plot_width,
      height = plot_height
    )
    
    output$download_data <- downloadHandler(
      filename = function() {"clonotype_expand.csv"},
      content = function(file) {
        write_csv(clonotypeTally(), file)
      }
    )
    
    output$download_plot <- downloadHandler(
      filename = function() {"clonotype_expand.pdf"},
      content = function(file) {
        ggsave(file, plot = clonotypeExpandPlot(), width = plot_width(), height = plot_height(), unit = "px", dpi = "screen")
      }
    )
    
  })
}