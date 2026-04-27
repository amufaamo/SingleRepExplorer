trajectoryInferenceUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("Trajectory Inference Settings"),
      p("Build pseudotime lineages using slingshot based on your Dimensionality Reduction embeddings."),
      selectInput(ns("dim_red"), "Dimensionality Reduction", choices = c("umap", "pca", "tsne"), selected = "umap"),
      selectInput(ns("clusterLabels"), "Cluster Annotations", choices = NULL),
      p(em("Hint: You can select annotative columns like seurat_clusters, Celltype, or singler_celltype.")),
      
      # Starting cluster
      selectInput(ns("startCluster"), "Starting Cluster / Root (Optional)", choices = c("Auto" = ""), selected = ""),
      
      # Ending cluster (new)
      selectInput(ns("endCluster"), "Ending Cluster / Tip (Optional)", choices = c("Auto" = ""), selected = ""),
      
      # Cluster subset selection (new)
      hr(),
      checkboxInput(ns("use_subset"), "Use Only Selected Clusters", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s'] == true", ns("use_subset")),
        checkboxGroupInput(ns("subset_clusters"), "Select Clusters to Include:", choices = NULL, inline = FALSE)
      ),
      hr(),
      
      actionButton(ns("run_slingshot"), "Run Slingshot", class = "btn-primary", width = "100%", icon = icon("project-diagram")),
      hr(),
      
      h4("Plot Options"),
      selectInput(ns("colorBy"), "Color Cells By:", choices = c("Clusters", "Pseudotime", "Clonotype Expansion", "Custom Gene")),
      
      # For Custom Gene
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Custom Gene'", ns("colorBy")),
        selectizeInput(ns("custom_gene"), "Gene Name", choices = NULL,
                       options = list(placeholder = "e.g., CD8A", maxOptions = 50))
      ),
      
      # For Clonotype Expansion
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Clonotype Expansion'", ns("colorBy")),
        selectInput(ns("cloneType"), "Repertoire Type", choices = c("TCR", "BCR"), selected = "TCR")
      ),
      
      numericInput(ns("pointSize"), "Point Size", value = 1.0, min = 0.1, max = 5, step = 0.1),
      numericInput(ns("curveWidth"), "Curve Line Width", value = 1.0, min = 0.1, max = 5, step = 0.1),
      selectInput(ns("color_palette"), "Color Palette",
                  choices = c("Magma (C)" = "C", "Viridis (D)" = "D", "Plasma (E)" = "E", "Inferno (B)" = "B", "Cividis" = "cividis", "Turbo" = "turbo"),
                  selected = "turbo"),
      hr(),
      commonPlotOptions(ns)
    ),

    mainPanel(
      h3("Pseudotime Trajectory Map"),
      p("Displays the branching lineages and continuous developmental trajectories. The principal curves trace the centers of dense cell paths."),
      shinyjs::useShinyjs(),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      downloadButton(ns("download_pseudotime"), "Download Pseudotime Table (.xlsx)"),
      plotOutput(ns("trajectoryPlot"))
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
      choices <- choices[!grepl("^(RNA_snn|barcode|orig.ident|TCR_CT|BCR_CT|TCR_TRA|TCR_TRB)", choices)]
      updateSelectInput(session, "clusterLabels", choices = choices, selected = "seurat_clusters")
    })

    # Update gene list for Custom Gene autocomplete
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      genes <- rownames(myReactives$seurat_object)
      updateSelectizeInput(session, "custom_gene", choices = genes, server = TRUE)
    })
    
    # Update Start / End / Subset Cluster choices based on selected Cluster Labels
    observeEvent(input$clusterLabels, {
      req(myReactives$seurat_object, input$clusterLabels)
      so <- myReactives$seurat_object
      if (input$clusterLabels %in% names(so@meta.data)) {
         clusters <- sort(unique(na.omit(as.character(so@meta.data[[input$clusterLabels]]))))
         updateSelectInput(session, "startCluster", choices = c("Auto" = "", clusters))
         updateSelectInput(session, "endCluster",   choices = c("Auto" = "", clusters))
         updateCheckboxGroupInput(session, "subset_clusters", choices = clusters, selected = clusters)
      }
    })
    
    # Run Slingshot
    slingshot_data <- eventReactive(input$run_slingshot, {
      req(myReactives$seurat_object, input$clusterLabels)
      
      if (!requireNamespace("slingshot", quietly = TRUE)) {
         showNotification("Please install the 'slingshot' package.", type = "error")
         return(NULL)
      }
      
      so <- myReactives$seurat_object
      
      # Validate reduction exists
      if (!(input$dim_red %in% names(so@reductions))) {
         showNotification(paste("Reduction", input$dim_red, "not found."), type = "error")
         return(NULL)
      }
      
      embeddings    <- Embeddings(so, input$dim_red)
      clusters_all  <- as.character(so@meta.data[[input$clusterLabels]])
      clusters      <- clusters_all

      # --- Cluster subsetting ---
      kept_mask <- rep(TRUE, ncol(so))
      if (isTRUE(input$use_subset) && length(input$subset_clusters) > 0) {
        kept_mask <- clusters_all %in% input$subset_clusters
        if (sum(kept_mask) < 10) {
          showNotification("Too few cells remain after subsetting clusters. Please select more clusters.", type = "warning")
          return(NULL)
        }
        embeddings <- embeddings[kept_mask, , drop = FALSE]
        clusters   <- clusters_all[kept_mask]
      }

      start_clus <- if (nzchar(input$startCluster)) input$startCluster else NULL
      end_clus   <- if (nzchar(input$endCluster))   input$endCluster   else NULL

      withProgress(message = "Running Slingshot inference...", value = 0.5, {
        tryCatch({
          sds <- slingshot::slingshot(
            data          = embeddings,
            clusterLabels = clusters,
            start.clus    = start_clus,
            end.clus      = end_clus
          )
          return(list(sds = sds, kept_mask = kept_mask))
        }, error = function(e) {
          showNotification(paste("Slingshot Error:", e$message), type = "error")
          return(NULL)
        })
      })
    })
    
    # Plotting
    plot_reactive <- reactive({
      res <- slingshot_data()
      req(res, myReactives$seurat_object)
      sds <- res$sds
      
      so         <- myReactives$seurat_object
      embeddings <- Embeddings(so, input$dim_red) %>% as.data.frame()
      colnames(embeddings)[1:2] <- c("Dim1", "Dim2")
      embeddings$Cluster <- so@meta.data[[input$clusterLabels]]
      
      # Apply cluster subset mask
      if (isTRUE(input$use_subset) && length(input$subset_clusters) > 0) {
        keep <- as.character(embeddings$Cluster) %in% input$subset_clusters
        embeddings <- embeddings[keep, , drop = FALSE]
      }
      
      # Handle coloring
      if (input$colorBy == "Clusters") {
        # Preserve Seurat's original cluster level order so colors match the UMAP
        orig_levels <- levels(so@meta.data[[input$clusterLabels]])
        if (is.null(orig_levels)) orig_levels <- sort(unique(as.character(so@meta.data[[input$clusterLabels]])))
        embeddings$ColorValue <- factor(embeddings$Cluster, levels = orig_levels)
        color_label <- input$clusterLabels
      } else if (input$colorBy == "Pseudotime") {
        pt <- slingshot::slingPseudotime(sds)
        # Handle matrix return from slingPseudotime, average across lineages or pick primary if possible
        if(is.matrix(pt)) {
           embeddings$ColorValue <- rowMeans(pt, na.rm = TRUE)
        } else {
           embeddings$ColorValue <- pt
        }
        color_label <- "Avg Pseudotime"
      } else if (input$colorBy == "Clonotype Expansion") {
        col_name <- paste0(input$cloneType, "_cloneSize")
        if (col_name %in% names(so@meta.data)) {
          meta_sub <- so@meta.data
          if (isTRUE(input$use_subset) && length(input$subset_clusters) > 0) {
            keep <- as.character(meta_sub[[input$clusterLabels]]) %in% input$subset_clusters
            meta_sub <- meta_sub[keep, , drop = FALSE]
          }
          embeddings$ColorValue <- factor(meta_sub[[col_name]])
        } else {
          embeddings$ColorValue <- "Not Available (Load TCR/BCR)"
        }
        color_label <- paste(input$cloneType, "Expansion")
      } else if (input$colorBy == "Custom Gene") {
        req(input$custom_gene)
        gene <- input$custom_gene
        if (gene %in% rownames(so)) {
          expr <- GetAssayData(so, layer = "data")[gene, ]
          if (isTRUE(input$use_subset) && length(input$subset_clusters) > 0) {
            keep <- as.character(so@meta.data[[input$clusterLabels]]) %in% input$subset_clusters
            expr <- expr[keep]
          }
          embeddings$ColorValue <- expr
        } else {
          embeddings$ColorValue <- NA
        }
        color_label <- paste("Expression:", input$custom_gene)
      }
      
      # Base plot
      p <- ggplot(embeddings, aes(x = Dim1, y = Dim2))
      
      palette  <- input$color_palette %||% "C"
      leg_pos  <- input$legend        %||% "right"

      if (is.numeric(embeddings$ColorValue)) {
         p <- p + geom_point(aes(color = ColorValue), size = input$pointSize, alpha = 0.8) +
              scale_color_viridis_c(name = color_label, option = palette)
      } else if (input$colorBy == "Clusters") {
         p <- p + geom_point(aes(color = ColorValue), size = input$pointSize, alpha = 0.8) +
              scale_color_hue(name = color_label) +
              guides(color = guide_legend(override.aes = list(size = 3)))
      } else {
         p <- p + geom_point(aes(color = ColorValue), size = input$pointSize, alpha = 0.8) +
              scale_color_viridis_d(name = color_label, option = palette, na.value = "grey") +
              guides(color = guide_legend(override.aes = list(size = 3)))
      }

      p <- p + theme_classic() +
        theme(legend.position = leg_pos) +
        labs(
          title = paste("Slingshot Trajectory on", toupper(input$dim_red)),
          subtitle = if (isTRUE(input$use_subset)) paste("Clusters:", paste(sort(input$subset_clusters), collapse = ", ")) else NULL,
          x = paste0(toupper(input$dim_red), "_1"),
          y = paste0(toupper(input$dim_red), "_2")
        )
      
      # Extract Principal Curves
      curves <- slingshot::slingCurves(sds)
      for (i in seq_along(curves)) {
         curve_coords <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
         colnames(curve_coords)[1:2] <- c("Dim1", "Dim2")
         p <- p + geom_path(data = curve_coords, aes(x = Dim1, y = Dim2),
                            linewidth = input$curveWidth, color = "black", inherit.aes = FALSE)
      }
      
      return(p)
    })
    
    output$trajectoryPlot <- renderPlot({
      plot_reactive()
    }, width = reactive(input$plot_width %||% 800), height = reactive(input$plot_height %||% 600))

    output$download_plot <- downloadHandler(
      filename = function() { "Slingshot_Trajectory.pptx" },
      content = function(file) {
        save_plot_as_pptx(file, plot_reactive(), input$plot_width, input$plot_height)
      }
    )

    output$download_pseudotime <- downloadHandler(
      filename = function() { "pseudotime_table.xlsx" },
      content = function(file) {
        res <- slingshot_data()
        req(res, myReactives$seurat_object)
        so  <- myReactives$seurat_object
        pt  <- slingshot::slingPseudotime(res$sds)
        barcodes <- if (isTRUE(input$use_subset) && length(input$subset_clusters) > 0) {
          rownames(so@meta.data)[res$kept_mask]
        } else {
          rownames(so@meta.data)
        }
        df_pt <- as.data.frame(pt)
        df_pt$barcode <- barcodes
        df_pt <- dplyr::select(df_pt, barcode, dplyr::everything())
        openxlsx::write.xlsx(df_pt, file)
      }
    )
  })
}
