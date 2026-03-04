# CoNGA-inspired Multimodal Integration Module
# Analyzes the overlap between Gene Expression Graph and Repertoire Sequence Graph

library(shiny)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringdist)

# --- UI Definition ---
congaIntegrationUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("CoNGA-style Integration"),
      p("Identify cells with coordinated phenotypic and clonal profiles by comparing Transcriptome KNN and Repertoire KNN graphs."),
      selectInput(ns("vdj_type"), "Repertoire Type:", choices = c("TCR", "BCR"), selected = "TCR"),
      selectInput(ns("chain_type"), "Sequence to use:", choices = c("Paired CDR3 (CTaa)", "Beta/Heavy chain", "Alpha/Light chain"), selected = "Paired CDR3 (CTaa)"),
      numericInput(ns("knn_k"), "K for Nearest Neighbors:", value = 15, min = 5, max = 50),
      actionButton(ns("run_conga"), "Run Integration", icon = icon("project-diagram"), class = "btn-primary", width = "100%"),
      hr(),
      p("The CoNGA Score represents the number of shared nearest neighbors between the gene expression and sequence similarity graphs.")
    ),
    mainPanel(
      h3("Multimodal Convergence Map"),
      p("Cells with high scores have both similar gene expression and similar VDJ sequences compared to their neighbors."),
      plotOutput(ns("conga_umap"), height = "600px"),
      hr(),
      h4("CoNGA Score Distribution by Cluster"),
      plotOutput(ns("conga_violin"), height = "400px")
    )
  )
}

# --- Server Logic ---
congaIntegrationServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    conga_results <- reactiveVal(NULL)
    
    observeEvent(input$run_conga, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      
      # Determine sequence column
      seq_col <- ""
      if (input$vdj_type == "TCR") {
        if (input$chain_type == "Paired CDR3 (CTaa)") seq_col <- "TCR_pair_CTaa"
        else if (input$chain_type == "Beta/Heavy chain") seq_col <- "TCR_TRB_cdr3"
        else seq_col <- "TCR_TRA_cdr3"
      } else {
        if (input$chain_type == "Paired CDR3 (CTaa)") seq_col <- "BCR_pair_CTaa"
        else if (input$chain_type == "Beta/Heavy chain") seq_col <- "BCR_IGH_cdr3_aa"
        else seq_col <- "BCR_IGK_cdr3_aa" # Simplification for heavy/light
      }
      
      if (!seq_col %in% names(so@meta.data)) {
        df_rep <- if(input$vdj_type == "TCR") myReactives$tcr_df else myReactives$bcr_df
        if (!is.null(df_rep) && seq_col %in% names(df_rep) && "barcode" %in% names(df_rep)) {
           # Fetch from repertoire df
           temp_df <- df_rep[, c("barcode", seq_col)]
           # Remove duplicates if any
           temp_df <- temp_df[!duplicated(temp_df$barcode), ]
           rownames(temp_df) <- temp_df$barcode
           common_cells <- intersect(colnames(so), temp_df$barcode)
           if (length(common_cells) > 0) {
               so <- AddMetaData(so, metadata = temp_df[common_cells, seq_col], col.name = seq_col)
           } else {
               showNotification(paste("Sequence data", seq_col, "found but no common cells mapped."), type = "error")
               return()
           }
        } else {
           showNotification(paste("Sequence data", seq_col, "not found in metadata or Repertoire dataframe. Please upload Repertoire data."), type = "error")
           return()
        }
      }
      
      withProgress(message = "Running Multimodal Integration...", value = 0, {
        
        # 1. Filter cells with valid sequences
        incProgress(0.1, detail = "Extracting valid cells...")
        valid_cells <- rownames(so@meta.data)[!is.na(so@meta.data[[seq_col]]) & so@meta.data[[seq_col]] != "NA_NA" & so@meta.data[[seq_col]] != ""]
        if (length(valid_cells) < 10) {
          showNotification("Not enough cells with valid sequences.", type = "error")
          return()
        }
        
        valid_seqs <- so@meta.data[valid_cells, seq_col]
        
        # 2. Extract Transcriptome KNN
        incProgress(0.3, detail = "Computing Transcriptome Graph...")
        pca_dims <- min(20, ncol(Embeddings(so, "pca")))
        if(pca_dims < 2) {
           showNotification("PCA not computed on Seurat Object.", type="error")
           return()
        }
        pca_emb <- Embeddings(so, "pca")[valid_cells, 1:pca_dims]
        rna_dist <- as.matrix(dist(pca_emb))
        k <- input$knn_k
        
        # Build RNA adjacency
        rna_adj <- matrix(0, nrow = length(valid_cells), ncol = length(valid_cells))
        for(i in 1:nrow(rna_dist)) {
           neighbors <- order(rna_dist[i,])[2:(k+1)]
           rna_adj[i, neighbors] <- 1
        }
        
        # 3. Compute Sequence KNN
        incProgress(0.6, detail = "Computing Sequence Graph (Levenshtein)...")
        if (length(valid_cells) > 5000) {
            showNotification("Optimization: Calculating distances for many cells. May take some time.", type = "warning", duration=5)
        }
        
        tryCatch({
           seq_dist <- stringdist::stringdistmatrix(valid_seqs, valid_seqs, method = "lv")
           
           seq_adj <- matrix(0, nrow = length(valid_cells), ncol = length(valid_cells))
           for(i in 1:nrow(seq_dist)) {
              neighbors <- order(seq_dist[i,])[2:(k+1)]
              seq_adj[i, neighbors] <- 1
           }
           
           incProgress(0.8, detail = "Intersecting Graphs...")
           # Intersection of adjacency = element-wise multiplication
           overlap_adj <- rna_adj * seq_adj
           
           # Score = degree in overlap graph (symmetric sum)
           overlap_symmetric <- pmax(overlap_adj, t(overlap_adj))
           conga_scores <- rowSums(overlap_symmetric)
           
           # Save back to Seurat object
           new_col <- paste0("CoNGA_Score_", input$vdj_type)
           
           current_scores <- rep(NA, ncol(so))
           names(current_scores) <- colnames(so)
           current_scores[valid_cells] <- conga_scores
           
           so <- AddMetaData(so, current_scores, col.name = new_col)
           myReactives$seurat_object <- so
           
           conga_results(list(col = new_col, k = k))
           showNotification("CoNGA analysis complete!", type = "message")
           
        }, error = function(e){
           showNotification(paste("Error in sequence computation:", e$message), type="error")
        })
        
        incProgress(1)
      })
    })
    
    output$conga_umap <- renderPlot({
      res <- conga_results()
      so <- myReactives$seurat_object
      
      if (is.null(res)) {
         return(DimPlot(so, reduction = "umap") + ggtitle("Run Integration to visualize scores"))
      }
      
      score_col <- res$col
      k <- res$k
      
      FeaturePlot(so, features = score_col, reduction = "umap", order = TRUE, pt.size = 1.2) +
        scale_color_viridis_c(option = "magma", direction = -1, name = "CoNGA Score\n(Shared Neighbors)", na.value="grey90") +
        labs(title = "Multimodal CoNGA Score Projection", subtitle = paste("Based on intersecting KNN graphs (k =", k, ")")) +
        theme(plot.title = element_text(face="bold"))
    })
    
    output$conga_violin <- renderPlot({
      res <- conga_results()
      req(res, myReactives$seurat_object)
      so <- myReactives$seurat_object
      
      score_col <- res$col
      
      # Determine clustering column to group by
      meta <- so@meta.data
      categorical_cols <- names(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      group_col <- if ("seurat_clusters" %in% categorical_cols) "seurat_clusters" else categorical_cols[1]
      
      VlnPlot(so, features = score_col, group.by = group_col, pt.size = 0) +
         geom_boxplot(width = 0.2, fill = "white") +
         theme(legend.position = "none") +
         labs(x = "Cluster", y = "CoNGA Score", title = paste("Scores by", group_col))
    })
  })
}
