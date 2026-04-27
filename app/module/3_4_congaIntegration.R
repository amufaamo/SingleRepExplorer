# CoNGA-inspired Multimodal Integration Module
# Analyzes the overlap between Gene Expression Graph and Repertoire Sequence Graph

# --- Helper: build dual-panel patchwork ---
.make_conga_dual <- function(so, score_col, k, threshold, pt_size, legend_pos, base_size) {

  # Panel A: CoNGA score gradient
  p_left <- FeaturePlot(so, features = score_col, reduction = "umap",
                        order = TRUE, pt.size = pt_size) +
    scale_color_viridis_c(option = "magma", direction = -1,
                          name = "CoNGA Score\n(Shared Neighbors)", na.value = "grey90") +
    labs(title    = "CoNGA Score Projection",
         subtitle = paste("Intersecting KNN graphs (k =", k, ")")) +
    theme_classic(base_size = base_size) +
    theme(plot.title    = element_text(face = "bold"),
          legend.position = legend_pos)

  # Panel B: Receptor-driven cells highlighted by cluster
  umap_df <- as.data.frame(Embeddings(so, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$score   <- so@meta.data[rownames(umap_df), score_col]

  # Use seurat_clusters if present; otherwise first categorical column
  meta <- so@meta.data
  cat_cols  <- names(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
  group_col <- if ("seurat_clusters" %in% cat_cols) "seurat_clusters" else cat_cols[1]
  umap_df$cluster <- as.factor(meta[rownames(umap_df), group_col])

  umap_df$is_high <- !is.na(umap_df$score) & umap_df$score >= threshold

  df_bg <- umap_df[!umap_df$is_high, ]
  df_hi <- umap_df[umap_df$is_high, ]
  n_hi  <- sum(umap_df$is_high, na.rm = TRUE)
  n_tot <- nrow(umap_df)

  p_right <- ggplot() +
    geom_point(data = df_bg, aes(UMAP_1, UMAP_2),
               color = "grey85", size = pt_size, alpha = 0.5) +
    geom_point(data = df_hi, aes(UMAP_1, UMAP_2, color = cluster),
               size = pt_size * 1.8) +
    scale_color_hue(name = group_col) +
    labs(title    = paste0("Receptor-Driven Cells (Score \u2265 ", threshold, ")"),
         subtitle = paste0(n_hi, " / ", n_tot, " cells"),
         x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = base_size) +
    theme(plot.title     = element_text(face = "bold"),
          legend.position = legend_pos)

  p_left | p_right
}

# --- UI Definition ---
congaIntegrationUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("CoNGA-style Integration"),
      p("Identify cells with coordinated phenotypic and clonal profiles by comparing Transcriptome KNN and Repertoire KNN graphs."),
      selectInput(ns("vdj_type"), "Repertoire Type:", choices = c("TCR", "BCR"), selected = "TCR"),
      selectInput(ns("chain_type"), "Sequence to use:",
                  choices = c("Paired CDR3 (CTaa)", "Beta/Heavy chain", "Alpha/Light chain"),
                  selected = "Paired CDR3 (CTaa)"),
      numericInput(ns("knn_k"), "K for Nearest Neighbors:", value = 15, min = 5, max = 50),
      actionButton(ns("run_conga"), "Run Integration",
                   icon = icon("project-diagram"), class = "btn-primary", width = "100%"),
      hr(),
      p("The CoNGA Score = number of shared nearest neighbors between the gene expression and sequence similarity graphs."),
      hr(),
      h5("Plot Options"),
      numericInput(ns("highlight_threshold"), "Highlight Threshold (Score \u2265):",
                   value = 1, min = 0, step = 1),
      commonPlotOptions(ns, legend_selected = "right", width_value = 1200, height_value = 500),
      pointSizeInput(ns, value = 1.2, min = 0.1, max = 5, step = 0.1)
    ),
    mainPanel(
      h3("Multimodal Convergence Map"),
      p("Left: CoNGA score gradient. Right: Cells above threshold colored by cluster (others grey)."),
      downloadButton(ns("download_conga_umap"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns("conga_umap"), height = "500px"),
      hr(),
      h4("CoNGA Score Distribution by Cluster"),
      downloadButton(ns("download_conga_violin"), "Download Plot (.pptx)"),
      br(), br(),
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
        else seq_col <- "BCR_IGK_cdr3_aa"
      }

      if (!seq_col %in% names(so@meta.data)) {
        df_rep <- if (input$vdj_type == "TCR") myReactives$tcr_df else myReactives$bcr_df
        if (!is.null(df_rep) && seq_col %in% names(df_rep) && "barcode" %in% names(df_rep)) {
          temp_df <- df_rep[, c("barcode", seq_col)]
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
          showNotification(paste("Sequence data", seq_col, "not found. Please upload Repertoire data."), type = "error")
          return()
        }
      }

      withProgress(message = "Running Multimodal Integration...", value = 0, {

        # 1. Filter cells with valid sequences
        incProgress(0.1, detail = "Extracting valid cells...")
        valid_cells <- rownames(so@meta.data)[
          !is.na(so@meta.data[[seq_col]]) &
          so@meta.data[[seq_col]] != "NA_NA" &
          so@meta.data[[seq_col]] != ""
        ]
        if (length(valid_cells) < 10) {
          showNotification("Not enough cells with valid sequences.", type = "error")
          return()
        }
        n_cells    <- length(valid_cells)
        valid_seqs <- so@meta.data[valid_cells, seq_col]
        k          <- input$knn_k

        # 2. RNA KNN using RANN (O(N*k) memory — no N×N distance matrix)
        incProgress(0.25, detail = "Computing Transcriptome KNN (RANN)...")
        pca_dims <- min(20, ncol(Embeddings(so, "pca")))
        if (pca_dims < 2) {
          showNotification("PCA not computed on Seurat Object.", type = "error")
          return()
        }
        pca_emb    <- Embeddings(so, "pca")[valid_cells, 1:pca_dims, drop = FALSE]
        rna_nn_idx <- RANN::nn2(pca_emb, k = k + 1)$nn.idx[, -1, drop = FALSE]  # n_cells × k
        rm(pca_emb); gc()

        # 3. Sequence KNN — operate on unique sequences only
        incProgress(0.45, detail = "Building Sequence Index...")
        unique_seqs  <- unique(valid_seqs)
        n_unique     <- length(unique_seqs)
        cell_seq_idx <- match(valid_seqs, unique_seqs)

        seq_nn_idx <- NULL

        tryCatch({
          if (n_unique <= 5000) {
            incProgress(0.55, detail = paste0("Levenshtein on ", n_unique, " unique sequences..."))
            seq_dist_mat <- stringdist::stringdistmatrix(unique_seqs, unique_seqs, method = "lv")
            k_seq        <- min(k, n_unique - 1L)
            seq_nn_idx   <- t(apply(seq_dist_mat, 1L, function(d) head(order(d)[-1L], k_seq)))
            rm(seq_dist_mat); gc()
          } else {
            showNotification(
              paste0(n_unique, " unique sequences (>5000). Using exact-match clonal CoNGA. ",
                     "Subset data for full Levenshtein analysis."),
              type = "warning", duration = 12
            )
          }
        }, error = function(e) {
          showNotification(paste("Sequence distance error:", e$message), type = "error")
        })

        # 4. Intersect graphs (index-based — no dense adjacency matrices)
        incProgress(0.75, detail = "Intersecting Graphs...")
        conga_scores <- vapply(seq_len(n_cells), function(i) {
          rna_nbrs <- rna_nn_idx[i, ]
          if (!is.null(seq_nn_idx)) {
            nn_unique <- unique(c(cell_seq_idx[i], seq_nn_idx[cell_seq_idx[i], ]))
            nn_unique <- nn_unique[!is.na(nn_unique)]
            seq_nbrs  <- which(cell_seq_idx %in% nn_unique)
            seq_nbrs  <- seq_nbrs[seq_nbrs != i]
          } else {
            seq_nbrs <- which(cell_seq_idx == cell_seq_idx[i])
            seq_nbrs <- seq_nbrs[seq_nbrs != i]
          }
          length(intersect(rna_nbrs, seq_nbrs))
        }, FUN.VALUE = 0L)

        rm(rna_nn_idx); gc()

        # 5. Save to Seurat object
        new_col <- paste0("CoNGA_Score_", input$vdj_type)
        current_scores <- rep(NA_real_, ncol(so))
        names(current_scores) <- colnames(so)
        current_scores[valid_cells] <- conga_scores
        so <- AddMetaData(so, current_scores, col.name = new_col)
        myReactives$seurat_object <- so

        # Update threshold slider max based on actual max score
        max_score <- max(conga_scores, na.rm = TRUE)
        updateNumericInput(session, "highlight_threshold", value = 1L)

        conga_results(list(col = new_col, k = k))
        incProgress(1)
        showNotification("CoNGA analysis complete!", type = "message")
      })
    })

    # --- Dual-panel UMAP ---
    output$conga_umap <- renderPlot({
      res <- conga_results()
      so  <- myReactives$seurat_object
      req(so)

      validate(need(!is.null(res), "Run Integration to visualize CoNGA scores."))

      .make_conga_dual(
        so        = so,
        score_col = res$col,
        k         = res$k,
        threshold = input$highlight_threshold %||% 1,
        pt_size   = input$point_size %||% 1.2,
        legend_pos = input$legend %||% "right",
        base_size  = input$base_font_size %||% 14
      )
    }, width  = reactive(input$plot_width  %||% 1200),
       height = reactive(input$plot_height %||% 500))

    # --- Violin by cluster ---
    output$conga_violin <- renderPlot({
      res <- conga_results()
      req(res, myReactives$seurat_object)
      so <- myReactives$seurat_object

      meta      <- so@meta.data
      cat_cols  <- names(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      group_col <- if ("seurat_clusters" %in% cat_cols) "seurat_clusters" else cat_cols[1]

      VlnPlot(so, features = res$col, group.by = group_col, pt.size = 0) +
        geom_boxplot(width = 0.2, fill = "white") +
        theme(legend.position = "none") +
        labs(x = "Cluster", y = "CoNGA Score", title = paste("Scores by", group_col))
    }, width  = reactive(input$plot_width  %||% 1200),
       height = reactive(input$plot_height %||% 400))

    # --- Downloads ---
    output$download_conga_umap <- downloadHandler(
      filename = function() { "conga_umap.pptx" },
      content  = function(file) {
        res <- conga_results(); req(res, myReactives$seurat_object)
        p <- .make_conga_dual(
          so        = myReactives$seurat_object,
          score_col = res$col,
          k         = res$k,
          threshold = input$highlight_threshold %||% 1,
          pt_size   = input$point_size %||% 1.2,
          legend_pos = input$legend %||% "right",
          base_size  = input$base_font_size %||% 14
        )
        save_plot_as_pptx(file, p, input$plot_width %||% 1200, input$plot_height %||% 500)
      }
    )

    output$download_conga_violin <- downloadHandler(
      filename = function() { "conga_violin.pptx" },
      content  = function(file) {
        res <- conga_results(); req(res, myReactives$seurat_object)
        so        <- myReactives$seurat_object
        meta      <- so@meta.data
        cat_cols  <- names(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
        group_col <- if ("seurat_clusters" %in% cat_cols) "seurat_clusters" else cat_cols[1]
        p <- VlnPlot(so, features = res$col, group.by = group_col, pt.size = 0) +
          geom_boxplot(width = 0.2, fill = "white") +
          theme(legend.position = "none") +
          labs(x = "Cluster", y = "CoNGA Score", title = paste("Scores by", group_col))
        save_plot_as_pptx(file, p, input$plot_width %||% 1200, input$plot_height %||% 400)
      }
    )
  })
}
