# Subset and Interactive Cell Selection Module
# Allows users to select cells on UMAP/TSNE via Lasso and create new subsets

library(shiny)
library(Seurat)
library(plotly)
library(dplyr)
library(DT)
library(bslib)

# --- UI Definition ---
subsetUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("Interactive Subsetting"),
      p("Lasso-select cells on the plot to define a new subset."),
      selectInput(ns("reduction"), "Reduction:", choices = c("umap", "tsne", "pca"), selected = "umap"),
      selectInput(ns("color_by"), "Color By:", choices = c("sample", "seurat_clusters", "sctype_celltype"), selected = "seurat_clusters"),
      hr(),
      h5("Selection Actions"),
      verbatimTextOutput(ns("selection_summary")),
      actionButton(ns("create_subset"), "Create Subset from Selected", icon = icon("cut"), class = "btn-secondary", width = "100%"),
      hr(),
      h5("Rerun Pipeline"),
      p("Rerunning the pipeline will normalize, find variable genes, and run PCA/UMAP on the subsetted data."),
      actionButton(ns("rerun_pipeline"), "Rerun Entire Pipeline", icon = icon("rocket"), class = "btn-primary", width = "100%"),
      hr(),
      downloadButton(ns("download_subset"), "Download Subset (.rds)")
    ),
    mainPanel(
      card(
        card_header("Interactive Dimensionality Reduction Plot"),
        plotlyOutput(ns("plotly_umap"), height = "600px")
      ),
      card(
        card_header("Selected Cells Info"),
        DT::DTOutput(ns("selected_cells_table"))
      )
    )
  )
}

# --- Server Logic ---
subsetServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 1. Update choices based on Seurat object
    observe({
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      
      # Update reduction choices
      reductions <- names(so@reductions)
      updateSelectInput(session, "reduction", choices = reductions, selected = if ("umap" %in% reductions) "umap" else reductions[1])
      
      # Update color_by choices (categorical metadata only)
      meta <- so@meta.data
      categorical_cols <- names(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "color_by", choices = categorical_cols, selected = if ("seurat_clusters" %in% categorical_cols) "seurat_clusters" else categorical_cols[1])
    })

    # 2. Render Plotly UMAP
    output$plotly_umap <- renderPlotly({
      req(myReactives$seurat_object, input$reduction, input$color_by)
      so <- myReactives$seurat_object
      
      # Extract coordinates and metadata
      embeds <- Embeddings(so, reduction = input$reduction)
      meta <- so@meta.data
      
      plot_df <- cbind(as.data.frame(embeds), meta[rownames(embeds), , drop = FALSE])
      plot_df$barcode <- rownames(plot_df)
      
      colnames(plot_df)[1:2] <- c("dim1", "dim2")
      
      p <- plot_ly(plot_df, 
                   x = ~dim1, y = ~dim2, 
                   color = ~.data[[input$color_by]], 
                   text = ~paste("Barcode:", barcode, "<br>", input$color_by, ":", .data[[input$color_by]]),
                   type = 'scatter', mode = 'markers', 
                   marker = list(size = 4, opacity = 0.7),
                   key = ~barcode,
                   source = "subset_plot") %>%
           layout(dragmode = "lasso", 
                  xaxis = list(title = colnames(embeds)[1]),
                  yaxis = list(title = colnames(embeds)[2]))
      
      p
    })

    # 3. Handle Selection
    selected_barcodes <- reactive({
      event <- event_data("plotly_selected", source = "subset_plot")
      if (is.null(event)) return(NULL)
      # plotly returns customdata or key if provided. Here we used 'key = ~barcode'
      return(event$key)
    })

    output$selection_summary <- renderText({
      bc <- selected_barcodes()
      if (is.null(bc)) "No cells selected yet." else paste("Cells selected:", length(bc))
    })

    output$selected_cells_table <- DT::renderDT({
      req(selected_barcodes())
      req(myReactives$seurat_object)
      
      meta <- myReactives$seurat_object@meta.data[selected_barcodes(), , drop = FALSE]
      DT::datatable(meta, options = list(pageLength = 5, scrollX = TRUE))
    })

    # 4. Create Subset
    observeEvent(input$create_subset, {
      req(selected_barcodes())
      withProgress(message = "Creating subset...", {
        so_subset <- subset(myReactives$seurat_object, cells = selected_barcodes())
        myReactives$seurat_object <- so_subset
        myReactives$meta.data <- so_subset@meta.data
        showNotification(paste("Subset created with", length(selected_barcodes()), "cells."), type = "message")
      })
    })

    # 5. Rerun Pipeline (Optional)
    observeEvent(input$rerun_pipeline, {
      req(myReactives$seurat_object)
      withProgress(message = "Rerunning analysis pipeline...", value = 0, {
        so <- myReactives$seurat_object
        
        incProgress(0.2, detail = "Normalizing...")
        so <- NormalizeData(so)
        
        incProgress(0.2, detail = "Variable Features...")
        so <- FindVariableFeatures(so, nfeatures = 2000)
        
        incProgress(0.2, detail = "Scaling & PCA...")
        so <- ScaleData(so)
        so <- RunPCA(so, npcs = 30)
        
        incProgress(0.2, detail = "Neighbors & Clusters...")
        so <- FindNeighbors(so, dims = 1:20)
        so <- FindClusters(so, resolution = 0.5)
        
        incProgress(0.2, detail = "UMAP...")
        so <- RunUMAP(so, dims = 1:20)
        
        myReactives$seurat_object <- so
        myReactives$meta.data <- so@meta.data
        showNotification("Pipeline rerun complete.", type = "message")
      })
    })

    # 6. Download
    output$download_subset <- downloadHandler(
        filename = function() { paste0("subset_", Sys.Date(), ".rds") },
        content = function(file) {
            saveRDS(myReactives$seurat_object, file)
        }
    )
  })
}
