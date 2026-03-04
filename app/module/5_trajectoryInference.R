library(shiny)
library(ggplot2)
library(dplyr)
library(shinyjs)
library(Seurat)

trajectoryInferenceUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("Trajectory Inference Settings"),
      p("Build pseudotime lineages using slingshot based on your Dimensionality Reduction embeddings."),
      selectInput(ns("dim_red"), "Dimensionality Reduction", choices = c("umap", "pca", "tsne"), selected = "umap"),
      selectInput(ns("clusterLabels"), "Cluster Annotations", choices = NULL),
      p(em("Hint: You can select annotative columns like seurat_clusters, Celltype, or singler_celltype.")),
      selectInput(ns("startCluster"), "Starting Cluster / Root (Optional)", choices = c("Auto" = ""), selected = ""),
      
      actionButton(ns("run_slingshot"), "Run Slingshot", class = "btn-primary", width = "100%", icon = icon("project-diagram")),
      hr(),
      
      h4("Visualization Settings"),
      selectInput(ns("colorBy"), "Color Cells By:", choices = c("Clusters", "Pseudotime", "Clonotype Expansion", "Custom Gene")),
      
      # For Custom Gene
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Custom Gene'", ns("colorBy")),
        textInput(ns("custom_gene"), "Gene Name (e.g., CD8A)", value = "")
      ),
      
      # For Clonotype Expansion
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Clonotype Expansion'", ns("colorBy")),
        selectInput(ns("cloneType"), "Repertoire Type", choices = c("TCR", "BCR"), selected = "TCR")
      ),
      
      numericInput(ns("pointSize"), "Point Size", value = 1.0, min = 0.1, max = 5, step = 0.1),
      numericInput(ns("curveWidth"), "Curve Line Width", value = 1.0, min = 0.1, max = 5, step = 0.1)
    ),
    
    mainPanel(
      h3("Pseudotime Trajectory Map"),
      p("Displays the branching lineages and continuous developmental trajectories. The principal curves trace the centers of dense cell paths."),
      shinyjs::useShinyjs(),
      downloadButton(ns("download_plot"), "Download Plot (PDF)"),
      plotOutput(ns("trajectoryPlot"), height = "600px")
    )
  )
}

trajectoryInferenceServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update dimension reduction choices dynamically
    observeEvent(myReactives$seurat_object, {
       req(myReactives$seurat_object)
       so <- myReactives$seurat_object
       reds <- names(so@reductions)
       if(length(reds) > 0) {
           updateSelectInput(session, "dim_red", choices = reds, selected = "umap")
       }
    })
    
    # Update cluster label choices dynamically
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      # Find categorical columns
      is_cat <- sapply(so@meta.data, function(x) is.character(x) || is.factor(x) || is.integer(x))
      choices <- names(so@meta.data)[is_cat]
      choices <- choices[!grepl("^(RNA_snn|barcode|orig.ident)", choices)]
      updateSelectInput(session, "clusterLabels", choices = choices, selected = "seurat_clusters")
    })
    
    # Update Start Cluster choices based on selected Cluster Labels
    observeEvent(input$clusterLabels, {
      req(myReactives$seurat_object, input$clusterLabels)
      so <- myReactives$seurat_object
      if (input$clusterLabels %in% names(so@meta.data)) {
         clusters <- unique(na.omit(as.character(so@meta.data[[input$clusterLabels]])))
         updateSelectInput(session, "startCluster", choices = c("Auto" = "", sort(clusters)))
      }
    })
    
    # Run Slingshot
    slingshot_data <- eventReactive(input$run_slingshot, {
      req(myReactives$seurat_object, input$clusterLabels)
      
      if (!requireNamespace("slingshot", quietly = TRUE)) {
         showNotification("Please install the 'slingshot' package (e.g., wait for v2.0 Docker to finish building).", type = "error")
         return(NULL)
      }
      
      so <- myReactives$seurat_object
      
      # Validate reduction exists
      if (!(input$dim_red %in% names(so@reductions))) {
         showNotification(paste("Reduction", input$dim_red, "not found."), type = "error")
         return(NULL)
      }
      
      embeddings <- Embeddings(so, input$dim_red)
      clusters <- so@meta.data[[input$clusterLabels]]
      
      start_clus <- if(input$startCluster != "") input$startCluster else NULL
      
      withProgress(message = "Running Slingshot inference...", value = 0.5, {
        tryCatch({
          sds <- slingshot::slingshot(data = embeddings, clusterLabels = as.character(clusters), start.clus = start_clus)
          return(sds)
        }, error = function(e) {
          showNotification(paste("Slingshot Error:", e$message), type = "error")
          return(NULL)
        })
      })
    })
    
    # Plotting
    plot_reactive <- reactive({
      # require slingshot data
      sds <- slingshot_data()
      req(sds, myReactives$seurat_object)
      
      so <- myReactives$seurat_object
      embeddings <- Embeddings(so, input$dim_red) %>% as.data.frame()
      colnames(embeddings)[1:2] <- c("Dim1", "Dim2")
      embeddings$Cluster <- so@meta.data[[input$clusterLabels]]
      
      # Handle coloring
      if (input$colorBy == "Clusters") {
        embeddings$ColorValue <- factor(embeddings$Cluster)
        color_label <- input$clusterLabels
      } else if (input$colorBy == "Pseudotime") {
        # get average pseudotime across lineages
        pt <- slingshot::slingPseudotime(sds)
        embeddings$ColorValue <- rowMeans(pt, na.rm = TRUE)
        color_label <- "Avg Pseudotime"
      } else if (input$colorBy == "Clonotype Expansion") {
        col_name <- paste0(input$cloneType, "_cloneSize")
        if (col_name %in% names(so@meta.data)) {
          embeddings$ColorValue <- factor(so@meta.data[[col_name]])
        } else {
          embeddings$ColorValue <- "Not Available (Load TCR/BCR)"
        }
        color_label <- paste(input$cloneType, "Expansion")
      } else if (input$colorBy == "Custom Gene") {
        req(input$custom_gene)
        gene <- input$custom_gene
        if (gene %in% rownames(so)) {
          expr <- GetAssayData(so, layer = "data")[gene, ]
          embeddings$ColorValue <- expr
          color_label <- paste("Expression:", gene)
        } else {
          embeddings$ColorValue <- NA
          color_label <- "Gene Not Found"
        }
      }
      
      # Base plot
      p <- ggplot(embeddings, aes(x = Dim1, y = Dim2))
      
      if (is.numeric(embeddings$ColorValue)) {
         p <- p + geom_point(aes(color = ColorValue), size = input$pointSize, alpha=0.8) +
              scale_color_viridis_c(name = color_label, option="C")
      } else {
         # try to use d3/ggsci colors if discrete
         p <- p + geom_point(aes(color = ColorValue), size = input$pointSize, alpha=0.8) +
              scale_color_viridis_d(name = color_label, option="turbo", na.value="grey") +
              guides(color = guide_legend(override.aes = list(size = 3)))
      }
      
      p <- p + theme_classic() + labs(title = paste("Slingshot Trajectory on", toupper(input$dim_red)), x=paste0(toupper(input$dim_red), "_1"), y=paste0(toupper(input$dim_red), "_2"))
      
      # Extract Principal Curves
      curves <- slingshot::slingCurves(sds)
      for (i in seq_along(curves)) {
         curve_coords <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
         colnames(curve_coords)[1:2] <- c("Dim1", "Dim2")
         p <- p + geom_path(data = curve_coords, aes(x = Dim1, y = Dim2), linewidth = input$curveWidth, color = "black", inherit.aes = FALSE)
      }
      
      return(p)
    })
    
    output$trajectoryPlot <- renderPlot({
      plot_reactive()
    })
    
    output$download_plot <- downloadHandler(
      filename = function() { "Slingshot_Trajectory.pdf" },
      content = function(file) {
        ggsave(file, plot = plot_reactive(), device = "pdf", width = 10, height = 8, units = "in")
      }
    )
  })
}
