library(shiny)
library(Seurat)
library(ggplot2)
library(shinyjs)

cellCommunicationUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      useShinyjs(),
      h4("Cell-Cell Communication (CellChat)"),
      p("Infer signaling networks between cell clusters and compare networks across clonal expansion states."),
      selectInput(ns("group_by"), "Cell Annotation (Clusters)", choices = c("seurat_clusters", "sctype_celltype", "singler_celltype", "Celltype")),
      p(em("Hint: You can compare Highly Expanded vs Unexpanded clones if you have loaded Repertoire data.")),
      actionButton(ns("run_cellchat"), "Run CellChat Inference", class = "btn-danger", width = "100%", icon = icon("network-wired")),
      hr(),
      h4("Visualization"),
      selectInput(ns("plot_type"), "Plot Type", choices = c("Network Circle", "Bubble Plot")),
      conditionalPanel(
          condition = sprintf("input['%s'] == 'Network Circle'", ns("plot_type")),
          selectInput(ns("pathway_show"), "Signaling Pathway", choices = NULL)
      ),
      downloadButton(ns("download_plot"), "Download Plot (.pdf)")
    ),
    mainPanel(
      h3("Inferred Communication Networks"),
      p("Note: Running CellChat may take several minutes. Ensure your data has distinct clusters and valid cell counts. Requires starting from v2.0 Docker container."),
      plotOutput(ns("networkPlot"), height = "700px")
    )
  )
}

cellCommunicationServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    cellchat_results <- reactiveVal(NULL)
    
    observeEvent(myReactives$seurat_object, {
        req(myReactives$seurat_object)
        so <- myReactives$seurat_object
        is_categorical <- sapply(so@meta.data, function(x) is.character(x) || is.factor(x))
        choices <- names(so@meta.data)[is_categorical]
        choices <- choices[!grepl("^(RNA_snn|barcode|orig.ident)", choices)]
        updateSelectInput(session, "group_by", choices = choices, selected = "seurat_clusters")
    })
    
    observeEvent(input$run_cellchat, {
      req(myReactives$seurat_object)
      if (!requireNamespace("CellChat", quietly = TRUE)) {
          showNotification("CellChat package is not installed. Please build v2.0 docker image.", type="error")
          return()
      }
      so <- myReactives$seurat_object
      
      withProgress(message = "Running CellChat...", value = 0, {
          tryCatch({
              incProgress(0.1, detail="Preparing data...")
              data.input <- Seurat::GetAssayData(so, assay = "RNA", layer = "data")
              meta <- so@meta.data
              
              if(!input$group_by %in% colnames(meta)) stop("Selected annotation not found in metadata.")
              meta$labels <- as.character(meta[[input$group_by]])
              # CellChat fails on NAs
              meta$labels[is.na(meta$labels)] <- "Unknown"
              
              incProgress(0.2, detail="Creating CellChat object...")
              cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")
              
              incProgress(0.3, detail="Setting Database (Human)...")
              CellChatDB <- CellChat::CellChatDB.human 
              CellChat::CellChatDB.use <- CellChat::subsetDB(CellChatDB, search = "Secreted Signaling")
              cellchat@DB <- CellChat::CellChatDB.use
              
              incProgress(0.4, detail="Preprocessing data...")
              cellchat <- CellChat::subsetData(cellchat)
              cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
              cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
              
              incProgress(0.6, detail="Computing Communication Probabilities...")
              cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
              cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
              
              incProgress(0.8, detail="Aggregating Network...")
              cellchat <- CellChat::computeCommunProbPathway(cellchat)
              cellchat <- CellChat::aggregateNet(cellchat)
              
              cellchat_results(cellchat)
              
              # Update pathway choices
              updateSelectInput(session, "pathway_show", choices = cellchat@netP$pathways)
              
              showNotification("CellChat completed successfully!", type="message")
              incProgress(1)
          }, error = function(e){
              showNotification(paste("CellChat Error:", e$message), type="error")
          })
      })
    })
    
    plot_reactive <- reactive({
        req(cellchat_results())
        cc <- cellchat_results()
        if (input$plot_type == "Bubble Plot") {
            # ggplot object
            return(CellChat::netVisual_bubble(cc, sources.use = unique(cc@idents), targets.use = unique(cc@idents), remove.isolate = FALSE))
        }
        return(NULL)
    })
    
    output$networkPlot <- renderPlot({
        req(cellchat_results())
        cc <- cellchat_results()
        if (input$plot_type == "Network Circle") {
             req(input$pathway_show)
             if(!input$pathway_show %in% cc@netP$pathways) {
                 plot.new()
                 text(0.5, 0.5, "Pathway not valid")
                 return()
             }
             CellChat::netVisual_aggregate(cc, signaling = input$pathway_show, layout = "circle")
        } else {
             plot_reactive()
        }
    })
    
    output$download_plot <- downloadHandler(
        filename = function() { "CellChat_Network.pdf" },
        content = function(file) {
           pdf(file, width=10, height=10)
           if(input$plot_type == "Network Circle"){
               req(input$pathway_show, cellchat_results())
               CellChat::netVisual_aggregate(cellchat_results(), signaling = input$pathway_show, layout = "circle")
           } else {
               print(plot_reactive())
           }
           dev.off()
        }
    )
  })
}
