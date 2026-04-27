# Sequence Similarity and Network Analysis Module

# --- UI Definition ---
sequenceSimilarityUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("Repertoire Similarity"),
      p("Analyze sequence similarity between clonotypes using edit distances."),
      vdjType(ns),
      numericInput(ns("top_n"), "Top N Clonotypes to Analyze:", value = 50, min = 10, max = 200),
      selectInput(ns("dist_method"), "Distance Method:", choices = c("Levenshtein" = "lv", "Hamming" = "hamming")),
      numericInput(ns("dist_threshold"), "Distance Threshold for Network:", min = 0, max = 100, value = 8, step = 1),
      hr(),
      h5("Network Display"),
      selectInput(ns("color_by"), "Color Nodes By:", choices = c("sample", "seurat_clusters", "TCR_pair_v_gene", "BCR_pair_v_identity"), selected = "sample"),
      actionButton(ns("run_analysis"), "Run Network Analysis", icon = icon("project-diagram"), class = "btn-primary", width = "100%"),
      hr(),
      h5("Plot Options"),
      commonPlotOptions(ns, legend_selected = "right", width_value = 800, height_value = 700),
      div(style = "display:flex; gap:10px;",
        numericInput(ns("node_size_min"), "Min Node Size", value = 2,  min = 1, max = 20, step = 1),
        numericInput(ns("node_size_max"), "Max Node Size", value = 10, min = 2, max = 30, step = 1)
      ),
      numericInput(ns("node_label_size"), "Node Label Size",  value = 3,   min = 1, max = 10, step = 0.5),
      numericInput(ns("edge_alpha"),      "Edge Transparency", value = 0.5, min = 0.1, max = 1, step = 0.1),
      numericInput(ns("heatmap_text_size"), "Heatmap Text Size", value = 6, min = 4, max = 14, step = 1)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sequence Network",
                 br(),
                 p(strong("Note:"), " Nodes are linked if their CDR3 amino acid sequences are within the edit distance threshold."),
                 downloadButton(ns("download_network"), "Download Plot (.pptx)"),
                 br(), br(),
                 plotOutput(ns("network_plot"), height = "auto")
        ),
        tabPanel("Distance Heatmap",
                 br(),
                 downloadButton(ns("download_heatmap"), "Download Plot (.pptx)"),
                 br(), br(),
                 plotOutput(ns("dist_heatmap"), height = "auto")
        )
      )
    )
  )
}

# --- Server Logic ---
sequenceSimilarityServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    network_results <- eventReactive(input$run_analysis, {
      req(input$vdj_type)
      if (!requireNamespace("stringdist", quietly = TRUE)) {
        showNotification("Package 'stringdist' is required for sequence similarity.", type = "error")
        return(NULL)
      }
      if (!requireNamespace("ggnetwork", quietly = TRUE)) {
        showNotification("Package 'ggnetwork' is required for network plotting.", type = "error")
        return(NULL)
      }
      if (!requireNamespace("igraph", quietly = TRUE)) {
        showNotification("Package 'igraph' is required for network plotting.", type = "error")
        return(NULL)
      }
      df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(df)

      seq_col <- if (input$vdj_type == "tcr") "TCR_pair_CTaa" else "BCR_pair_CTaa"
      if (!seq_col %in% names(df)) {
        showNotification(paste("Sequence column", seq_col, "not found."), type = "error")
        return(NULL)
      }
      
      # Ensure cluster and sample are present safely
      cluster_col <- "seurat_clusters"
      if (!cluster_col %in% names(df)) {
        if (!is.null(myReactives$seurat_object) && cluster_col %in% names(myReactives$seurat_object@meta.data)) {
           meta <- myReactives$seurat_object@meta.data[, c(cluster_col), drop=FALSE]
           meta$barcode <- rownames(meta)
           df <- df %>% left_join(meta, by="barcode")
        } else {
           df[[cluster_col]] <- "Unknown"
        }
      }
      if (!"sample" %in% names(df)) df$sample <- "Unknown"

      top_clones <- df %>%
        filter(!is.na(.data[[seq_col]]), .data[[seq_col]] != "", !is.na(raw_clonotype_id)) %>%
        group_by(raw_clonotype_id) %>%
        summarise(
          sequence = gsub("^.*_", "", dplyr::first(.data[[seq_col]])), # Focus on TRB/IGH chain for convergence
          count    = n(),
          sample   = dplyr::first(.data[["sample"]]),
          cluster  = dplyr::first(.data[[cluster_col]]),
          .groups  = "drop"
        ) %>%
        arrange(desc(count)) %>%
        head(input$top_n)

      if (nrow(top_clones) < 5) {
        showNotification("Too few clones for network analysis.", type = "warning")
        return(NULL)
      }

      dist_mat <- tryCatch({
        if (input$dist_method == "hamming") {
          seq_lengths <- unique(nchar(top_clones$sequence))
          if (length(seq_lengths) > 1) {
            stop("Hamming distance requires equal-length sequences.")
          }
        }
        stringdist::stringdistmatrix(top_clones$sequence, top_clones$sequence, method = input$dist_method)
      }, error = function(e) {
        showNotification(paste("Distance computation error:", e$message, "Falling back to Levenshtein."), type = "warning")
        tryCatch(
          stringdist::stringdistmatrix(top_clones$sequence, top_clones$sequence, method = "lv"),
          error = function(e2) {
            showNotification(paste("Distance computation failed:", e2$message), type = "error")
            NULL
          }
        )
      })
      if (is.null(dist_mat)) return(NULL)
      colnames(dist_mat) <- rownames(dist_mat) <- top_clones$raw_clonotype_id

      adj_mat <- dist_mat <= input$dist_threshold
      diag(adj_mat) <- 0

      g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")

      V(g)$count    <- top_clones$count
      V(g)$sample   <- top_clones$sample
      V(g)$cluster  <- as.character(top_clones$cluster)
      V(g)$clonotype <- top_clones$raw_clonotype_id

      return(list(graph = g, dist_mat = dist_mat))
    })

    build_network_plot <- function() {
      res <- network_results()
      if (is.null(res)) {
        return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Run network analysis to display the plot.") + theme_void())
      }
      g <- res$graph

      l_matrix <- igraph::layout_with_fr(g)
      n <- ggnetwork::ggnetwork(g, layout = l_matrix)
      color_var    <- if (input$color_by %in% names(n)) input$color_by else "sample"
      nd_min       <- input$node_size_min  %||% 2
      nd_max       <- input$node_size_max  %||% 10
      lbl_size     <- input$node_label_size %||% 3
      edge_al      <- input$edge_alpha %||% 0.5
      leg_pos      <- input$legend %||% "right"

      ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(color = "grey80", alpha = edge_al) +
        geom_nodes(aes(color = .data[[color_var]], size = count)) +
        geom_nodetext(aes(label = clonotype), repel = TRUE, size = lbl_size) +
        scale_size_continuous(range = c(nd_min, nd_max)) +
        theme_void() +
        theme(legend.position = leg_pos) +
        labs(title    = paste(toupper(input$vdj_type), "Clonotype Sequence Network"),
             subtitle = paste("Threshold:", input$dist_threshold, "Method:", input$dist_method),
             color    = input$color_by, size = "Cell Count")
    }

    build_heatmap_plot <- function() {
      res <- network_results()
      if (is.null(res)) {
        return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Run network analysis to display the heatmap.") + theme_void())
      }

      dist_mat <- res$dist_mat
      df_melt  <- as.data.frame(as.table(dist_mat))
      colnames(df_melt) <- c("Clone1", "Clone2", "Distance")

      txt_sz  <- input$heatmap_text_size %||% 6
      leg_pos <- input$legend %||% "right"

      ggplot(df_melt, aes(x = Clone1, y = Clone2, fill = Distance)) +
        geom_tile() +
        scale_fill_viridis_c(direction = -1) +
        theme_minimal() +
        theme(
          axis.text.x     = element_text(angle = 90, hjust = 1, size = txt_sz),
          axis.text.y     = element_text(size = txt_sz),
          legend.position = leg_pos
        ) +
        labs(title = "Pairwise Edit Distances")
    }

    output$network_plot <- renderPlot({
      build_network_plot()
    }, width = reactive(input$plot_width %||% 800), height = reactive(input$plot_height %||% 700))

    output$dist_heatmap <- renderPlot({
      build_heatmap_plot()
    }, width = reactive(input$plot_width %||% 800), height = reactive(input$plot_height %||% 700))

    output$download_network <- downloadHandler(
        filename = function() { paste0("clonotype_network_", Sys.Date(), ".pptx") },
        content = function(file) {
           p <- build_network_plot()
           save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
        }
    )

    output$download_heatmap <- downloadHandler(
        filename = function() { paste0("clonotype_heatmap_", Sys.Date(), ".pptx") },
        content = function(file) {
           p <- build_heatmap_plot()
           save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
        }
    )
  })
}
