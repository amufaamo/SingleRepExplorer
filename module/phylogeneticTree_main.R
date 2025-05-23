library(ggtree)
library(dowser)
library(alakazam)

phylogeneticTree_mainUI <- function(id) {
  ns <- NS(id)
  mainPanel(
    plotOutput(ns("tree")),
    downloadButton(ns("download_plot"), "Download plot (.pdf)")
  )
}


phylogeneticTree_mainServer <- function(id, bcr_IGH) {
  moduleServer(id, function(input, output, session) {
    
    plot_width <- reactive(input$plot_width)
    plot_height <- reactive(input$plot_height)
    
    phyloData <- eventReactive(input$run, {
      make_phylotree_data(bcr_IGH, clone = input$clone)
    })
    
    phyloTree <- reactive({
      clones <- dowser::formatClones(phyloData())
      trees  <- dowser::getTrees(clones)
      plots  <- dowser::plotTrees(trees)
      plots[[1]] +
        geom_tiplab() +
        # To avoid label cropping
        coord_cartesian(clip = 'off') + 
        theme(plot.margin = margin(6,200,6,6))
    })
    
    output$tree <- renderPlot(
      phyloTree(),
      width  = plot_width,
      height = plot_height
    )
    
    output$download_plot <- downloadHandler(
      filename = function() { paste0("phylogenetic_tree_", input$clone, ".pdf") },
      content = function(file) {
        ggsave(file, plot = phyloTree(), width = plot_width(), height = plot_height(), unit = "px", dpi = "screen")
      }
    )

    
  })
}