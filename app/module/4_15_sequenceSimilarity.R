# Sequence Similarity and Network Analysis Module
# Calculates CDR3 sequence distances and visualizes clonotype networks

library(shiny)
library(dplyr)
library(ggplot2)
library(igraph)
library(ggnetwork)
library(stringdist) # For Levenshtein distance

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
      sliderInput(ns("dist_threshold"), "Distance Threshold for Network:", min = 0, max = 10, value = 2, step = 1),
      hr(),
      h5("Network Display"),
      selectInput(ns("color_by"), "Color Nodes By:", choices = c("sample", "seurat_clusters", "TCR_pair_v_gene", "BCR_pair_v_identity"), selected = "sample"),
      actionButton(ns("run_analysis"), "Run Network Analysis", icon = icon("project-diagram"), class = "btn-primary", width = "100%"),
      hr(),
      downloadButton(ns("download_network"), "Download Network Plot")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sequence Network",
                 br(),
                 plotOutput(ns("network_plot"), height = "700px"),
                 p("Nodes are linked if their CDR3 sequences are within the edit distance threshold.")
        ),
        tabPanel("Distance Heatmap",
                 br(),
                 plotOutput(ns("dist_heatmap"), height = "600px")
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
      df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(df)
      
      # Determine sequence and ID columns
      seq_col <- if (input$vdj_type == "tcr") "TCR_pair_CTnt" else "BCR_pair_CTnt"
      if (!seq_col %in% names(df)) {
        showNotification(paste("Sequence column", seq_col, "not found."), type = "error")
        return(NULL)
      }

      # Prepare data: top N clonotypes
      top_clones <- df %>%
        filter(!is.na(.data[[seq_col]]) & .data[[seq_col]] != "") %>%
        group_by(raw_clonotype_id) %>%
        summarise(
          sequence = first(.data[[seq_col]]),
          count = n(),
          sample = first(sample), # Taking the first sample for simplicity in node coloring
          cluster = first(seurat_clusters),
          .groups = "drop"
        ) %>%
        arrange(desc(count)) %>%
        head(input$top_n)

      if (nrow(top_clones) < 5) {
        showNotification("Too few clones for network analysis.", type = "warning")
        return(NULL)
      }

      # Calculate distance matrix
      dist_mat <- stringdist::stringdistmatrix(top_clones$sequence, top_clones$sequence, method = input$dist_method)
      colnames(dist_mat) <- rownames(dist_mat) <- top_clones$raw_clonotype_id
      
      # Build Adjacency
      adj_mat <- dist_mat <= input$dist_threshold
      diag(adj_mat) <- 0 # Remove self-loops
      
      # Create Graph
      g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
      
      # Add attributes
      V(g)$count <- top_clones$count
      V(g)$sample <- top_clones$sample
      V(g)$cluster <- as.character(top_clones$cluster)
      V(g)$clonotype <- top_clones$raw_clonotype_id
      
      return(list(graph = g, dist_mat = dist_mat))
    })

    output$network_plot <- renderPlot({
      res <- network_results()
      req(res)
      g <- res$graph
      
      n <- ggnetwork(g, layout = "fruchtermanreingold")
      
      # Determine color variable
      color_var <- if (input$color_by %in% names(n)) input$color_by else "sample"

      # Base plot
      ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(color = "grey80", alpha = 0.5) +
        geom_nodes(aes(color = .data[[color_var]], size = count)) +
        geom_nodetext(aes(label = clonotype), repel = TRUE, size = 3) +
        scale_size_continuous(range = c(2, 10)) +
        theme_blank() +
        labs(title = paste(toupper(input$vdj_type), "Clonotype Sequence Network"), 
             subtitle = paste("Threshold:", input$dist_threshold, "Method:", input$dist_method),
             color = input$color_by, size = "Cell Count")
    })

    output$dist_heatmap <- renderPlot({
      res <- network_results()
      req(res)
      
      dist_mat <- res$dist_mat
      
      # Melt for ggplot
      df_melt <- as.data.frame(as.table(dist_mat))
      colnames(df_melt) <- c("Clone1", "Clone2", "Distance")
      
      ggplot(df_melt, aes(x = Clone1, y = Clone2, fill = Distance)) +
        geom_tile() +
        scale_fill_viridis_c(direction = -1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              axis.text.y = element_text(size = 6)) +
        labs(title = "Pairwise Edit Distances")
    })
    
    output$download_network <- downloadHandler(
        filename = function() { "clonotype_network.pdf" },
        content = function(file) {
           # Re-render plot for output
           ggsave(file, width = 10, height = 8)
        }
    )
  })
}
