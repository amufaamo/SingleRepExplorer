# Phylogenetic Tree Analysis Module

# Optional heavy packages — check availability at runtime
# NOTE: dowser dependency removed (2026-04-08); tree built directly via msa + ape + phangorn
.phylo_pkgs_ok <- requireNamespace("ggtree", quietly = TRUE) &&
                  requireNamespace("msa", quietly = TRUE) &&
                  requireNamespace("ape", quietly = TRUE) &&
                  requireNamespace("phangorn", quietly = TRUE)

# --- UI Definition ---
phylogeneticTreeUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("BCR Lineage Tree"),
      p("Visualize B-cell evolution within a clonotype."),
      selectInput(ns("selected_clone"), "Select Clone:", choices = NULL),
      actionButton(ns("run_phylotree"), "Build Tree", icon = icon("tree"), class = "btn-primary", width = "100%"),
      hr(),
      h5("Plot Options"),
      selectInput(ns("color_by"), "Color Tips By:", choices = c("None", "sample", "seurat_clusters", "BCR_pair_c_gene", "exact_subclonotype_id"), selected = "sample"),
      selectInput(ns("label_by"), "Label Tips By:", choices = c("barcode", "exact_subclonotype_id", "BCR_pair_CTaa"), selected = "barcode"),
      sliderInput(ns("label_size"), "Label Size:", min = 0, max = 10, value = 3, step = 0.5),
      sliderInput(ns("label_offset"), "Label Offset:", min = 0, max = 0.5, value = 0.05, step = 0.01),
      numericInput(ns("tip_size"), "Tip Point Size", value = 3, min = 0.5, max = 10, step = 0.5),
      numericInput(ns("tree_line_width"), "Tree Line Width", value = 0.8, min = 0.1, max = 5, step = 0.1),
      hr(),
      commonPlotOptions(ns, legend_selected = "right", width_value = 700, height_value = 600),
      hr(),
      # downloadButton moved to mainPanel
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Phylogenetic Tree", 
                 br(),
                 textOutput(ns("plot_info")),
                 br(),
                 downloadButton(ns("download_plot"), "Download Tree (.pptx)"),
                 br(), br(),
                 plotOutput(ns("plot"), width = "100%", height = "auto")
        ),
        tabPanel("Clone Data", 
                 br(),
                 tags$h4("Sequences in Selected Clone"),
                 DT::DTOutput(ns('table'))
        )
      )
    )
  )
}

# --- Server Logic ---
phylogeneticTreeServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Early warning if packages are not installed
    if (!.phylo_pkgs_ok) {
      output$plot_info <- renderText({
        "BCR Lineage Tree requires 'ggtree', 'msa', 'ape', and 'phangorn' packages. Please install them and restart the app."
      })
      output$plot <- renderPlot({
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "Required packages (ggtree/msa/ape/phangorn)\nnot installed. Contact administrator.", size = 6, hjust = 0.5) +
          theme_void()
      })
      return()
    }

    # 1. Update clone choices
    observe({
      req(myReactives$bcr_df)
      df <- myReactives$bcr_df

      # Allow selection of all clones; enforce sequence count at run time
      clones_summary <- tryCatch({
        df %>%
          dplyr::filter(!is.na(raw_clonotype_id) & raw_clonotype_id != "None") %>%
          dplyr::count(raw_clonotype_id, name = "n") %>%
          dplyr::filter(n >= 3) %>%
          dplyr::arrange(raw_clonotype_id)
      }, error = function(e) {
        return(tibble::tibble(raw_clonotype_id = character(0), n = integer(0)))
      })

      if (nrow(clones_summary) > 0) {
        choice_labels <- paste0(clones_summary$raw_clonotype_id, " (n=", clones_summary$n, ")")
        choices <- stats::setNames(clones_summary$raw_clonotype_id, choice_labels)
        current <- isolate(input$selected_clone)
        selected <- if (!is.null(current) && current %in% clones_summary$raw_clonotype_id) {
          current
        } else {
          clones_summary$raw_clonotype_id[1]
        }
        updateSelectInput(session, "selected_clone", choices = choices, selected = selected)
      } else {
        updateSelectInput(session, "selected_clone", choices = c("No suitable clones found" = ""), selected = "")
      }
    })

    # 2. Extract selected clone data
    selected_clone_data <- eventReactive(input$run_phylotree, {
      req(myReactives$bcr_df, input$selected_clone, input$selected_clone != "")
      selected_clonotype_id <- input$selected_clone
      df <- myReactives$bcr_df
      
      # Extract BCR Heavy chain nucleotide sequence
      nt_col <- "BCR_IGH_full_length_nt"
      if (!nt_col %in% names(df)) {
        seq_col_fallback <- "BCR_pair_CTnt"
        if (seq_col_fallback %in% names(df)) {
             df[[nt_col]] <- sapply(strsplit(as.character(df[[seq_col_fallback]]), "_"), function(x) if(length(x) >= 1) x[1] else NA)
        } else {
             showNotification("Sequencing nucleotide column not found in database.", type = "error")
             return(NULL)
        }
      }

      clone_subset <- df %>%
        dplyr::filter(
          raw_clonotype_id == selected_clonotype_id,
          !is.na(.data[[nt_col]]) & .data[[nt_col]] != "",
          !is.na(barcode) & barcode != ""
        ) %>%
        dplyr::distinct(barcode, .keep_all = TRUE)

      if (nrow(clone_subset) < 2) {
        showNotification("Insufficient sequences (< 2) for this clone.", type = "warning")
        return(NULL)
      }
      return(clone_subset)
    })

    # 3. Build phylogenetic tree directly from raw 10x sequences (no dowser/Change-O needed)
    phylogenetic_tree <- reactive({
      req(selected_clone_data())
      clone_subset <- selected_clone_data()
      nt_col <- "BCR_IGH_full_length_nt"

      withProgress(message = "Building tree...", value = 0, {
        tryCatch({
          # Extract sequences and labels
          seqs <- as.character(clone_subset[[nt_col]])
          labels <- as.character(clone_subset$barcode)
          names(seqs) <- labels

          incProgress(0.2, detail = "Aligning sequences...")

          # Multiple sequence alignment via msa (Biostrings DNAStringSet)
          # ClustalOmega writes temp files to cwd - switch to writable dir
          dna_ss <- Biostrings::DNAStringSet(seqs)
          names(dna_ss) <- labels
          old_wd <- getwd()
          setwd(tempdir())
          msa_result <- tryCatch(
            msa::msa(dna_ss, method = "ClustalOmega"),
            error = function(e) {
              message("[PhyloTree] ClustalOmega failed, trying Muscle: ", e$message)
              msa::msa(dna_ss, method = "Muscle")
            }
          )
          setwd(old_wd)

          # Convert to ape DNAbin
          aligned_seqs <- as.character(msa_result)  # named character matrix
          aligned_mat <- do.call(rbind, strsplit(as.character(Biostrings::unmasked(
            Biostrings::DNAMultipleAlignment(msa_result))), ""))
          rownames(aligned_mat) <- labels
          dnabin <- ape::as.DNAbin(aligned_mat)

          incProgress(0.5, detail = "Computing distances...")

          # Distance matrix and NJ tree
          d <- ape::dist.dna(dnabin, model = "N", pairwise.deletion = TRUE)
          # Handle zero-distance: add tiny jitter to prevent NJ collapse
          d[d == 0] <- 1e-8

          incProgress(0.7, detail = "Building NJ tree...")
          tree <- ape::nj(d)

          # Root at midpoint for nicer display
          tree <- tryCatch(phangorn::midpoint(tree), error = function(e) tree)

          # Fix negative edge lengths (NJ + midpoint rooting can produce them)
          if (any(tree$edge.length < 0)) {
            tree$edge.length[tree$edge.length < 0] <- 0
          }

          # Prepare metadata for ggtree %<+% join
          # ggtree uses column named 'label' as join key (must match tree$tip.label)
          meta_cols <- intersect(
            names(clone_subset),
            c("barcode", "sample", "seurat_clusters", "c_gene",
              "BCR_pair_c_gene", "exact_subclonotype_id", "BCR_pair_CTaa",
              "origin", "donor")
          )
          meta_df <- clone_subset[, meta_cols, drop = FALSE]
          # Set 'label' = barcode (= tip labels), then drop 'barcode' to avoid conflicts
          meta_df$label <- as.character(meta_df$barcode)
          meta_df$barcode <- NULL
          # Remove duplicates (ggtree %<+% needs unique labels)
          meta_df <- meta_df[!duplicated(meta_df$label), ]

          incProgress(1, detail = "Done")

          # Return tree (phylo) + metadata as a list
          list(tree = tree, meta = meta_df)
        }, error = function(e) {
          message("[PhyloTree] ERROR: ", e$message)
          showNotification(paste("Tree building error:", e$message), type = "error")
          NULL
        })
      })
    })

    # 5. Build plot
    tree_plot_object <- reactive({
      req(phylogenetic_tree())
      result <- phylogenetic_tree()
      tree_data <- result$tree
      meta_df   <- result$meta
      
      color_col <- input$color_by
      label_col <- input$label_by
      
      tip_sz   <- input$tip_size       %||% 3
      lw       <- input$tree_line_width %||% 0.8
      leg_pos  <- input$legend          %||% "right"

      # Build basic tree plot (without %<+% to avoid label duplication bugs)
      p <- ggtree::ggtree(tree_data) +
           ggtree::geom_tree(linewidth = lw) +
           ggtree::geom_treescale(x = 0, y = -1, fontsize = 3)

      # Manually merge metadata into plot data via left_join on 'label'
      # meta_df has 'label' column = barcode = tip.label
      plot_data <- p$data
      # Remove any meta columns that already exist in plot_data (except 'label')
      meta_cols_to_add <- setdiff(names(meta_df), names(plot_data))
      if (length(meta_cols_to_add) > 0) {
        join_df <- meta_df[, c("label", meta_cols_to_add), drop = FALSE]
        plot_data <- dplyr::left_join(plot_data, join_df, by = "label")
      }
      p$data <- plot_data

      if (color_col != "None" && color_col %in% names(p$data)) {
        p <- p + ggtree::geom_tippoint(aes(color = .data[[color_col]]), size = tip_sz)
      } else {
        p <- p + ggtree::geom_tippoint(size = tip_sz, color = "steelblue")
      }

      # 'barcode' was renamed to 'label' — map for display
      display_label_col <- if (label_col == "barcode") "label" else label_col
      if (display_label_col %in% names(p$data)) {
        p <- p + ggtree::geom_tiplab(aes(label = .data[[display_label_col]]),
                                     size   = input$label_size,
                                     offset = input$label_offset,
                                     align  = TRUE)
      }

      p <- p + theme(legend.position = leg_pos) +
           labs(title = paste("Clonotype:", input$selected_clone)) +
           ggtree::xlim_tree(max(p$data$x, na.rm = TRUE) * 1.5)
      
      return(p)
    })

    # 6. Outputs
    output$plot_info <- renderText({
      req(phylogenetic_tree())
      result <- phylogenetic_tree()
      n_tips <- ape::Ntip(result$tree)
      paste("Tree built for clone:", input$selected_clone, "| Sequences:", n_tips)
    })

    output$plot <- renderPlot({
      req(tree_plot_object())
      tree_plot_object()
    }, width  = reactive(input$plot_width  %||% 700),
       height = reactive({
         if (!is.null(input$plot_height) && input$plot_height != 600) return(input$plot_height)
         result <- phylogenetic_tree()
         n <- if (!is.null(result)) ape::Ntip(result$tree) else 30
         max(600, n * 20)
       }))

    output$table <- DT::renderDT({
      req(selected_clone_data())
      DT::datatable(selected_clone_data(), options = list(pageLength = 10, scrollX = TRUE))
    })

    output$download_plot <- downloadHandler(
      filename = function() { paste0("tree_", input$selected_clone, ".pptx") },
      content = function(file) {
        save_plot_as_pptx(file, tree_plot_object(), input$plot_width, input$plot_height)
      }
    )
  })
}
