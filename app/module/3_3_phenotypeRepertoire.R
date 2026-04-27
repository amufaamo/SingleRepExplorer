phenotypeRepertoireUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      shinyjs::useShinyjs(),
      h4("Phenotype-Repertoire Correlation Engine"),
      p("Identify genes whose expression directly correlates with clonotype expansion sizes."),
      selectInput(ns("repertoire_type"), "Repertoire Type",
                  choices = c("TCR", "BCR"), selected = "TCR"),
      selectInput(ns("cor_method"), "Correlation Method",
                  choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
                  selected = "spearman"),
      numericInput(ns("min_cells"), "Minimum Cells per Gene", value = 10, min = 1),
      numericInput(ns("p_val_adj"), "Adjusted P-value threshold (FDR)", value = 0.05,
                   min = 0, max = 1, step = 0.01),
      numericInput(ns("top_n"), "Show Top N significant genes in plot", value = 6,
                   min = 1, max = 20),
      actionButton(ns("run_cor"), "Run Correlation Engine",
                   class = "btn-primary", width = "100%", icon = icon("calculator")),
      hr(),
      h5("Plot Options"),
      div(style = "display: flex; gap: 10px;",
        numericInput(ns("plot_width"), "Plot Width (px)", min = 100, max = 2000, value = 700, step = 100),
        numericInput(ns("plot_height"), "Plot Height (px)", min = 100, max = 2000, value = 600, step = 100)
      ),
      pointSizeInput(ns, value = 0.7, min = 0.1, max = 5, step = 0.1),
      hr(),
      h5("Custom Gene List"),
      checkboxInput(ns("use_custom_genes"), "Use custom gene list", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s'] == true", ns("use_custom_genes")),
        textAreaInput(ns("custom_genes_text"),
                      label = "Gene names (one per line or comma-separated):",
                      placeholder = "e.g.\nCD8A\nGNLY\nTOX",
                      rows = 5)
      )
    ),
    mainPanel(
      h3("Gene Correlation with Clone Size"),
      p("Genes with positive correlation represent markers of clonal expansion (e.g., exhaustion, effector state)."),
      shinyjs::hidden(
        div(id = ns("download_panel"),
            downloadButton(ns("download_table"), "Download Full Results (.xlsx)"),
            br(), br()
        )
      ),
      DTOutput(ns("cor_table")),
      hr(),
      h3("Expression vs Clone Size (Top Significant Genes)"),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns("scatter_plots"))
    )
  )
}

phenotypeRepertoireServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cor_results <- reactiveVal(NULL)

    observeEvent(input$run_cor, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object

      col_name <- paste0(input$repertoire_type, "_cloneSize")
      if (!col_name %in% names(so@meta.data)) {
        showNotification(
          paste0(col_name, " not found in Seurat metadata.\n",
                 "Please upload ", input$repertoire_type, " data first (Data Upload tab)."),
          type = "error", duration = 8
        )
        return()
      }

      clone_sizes <- as.numeric(so@meta.data[[col_name]])
      valid_cells <- which(!is.na(clone_sizes) & clone_sizes > 0)

      if (length(valid_cells) < 10) {
        showNotification(
          paste0("Only ", length(valid_cells), " cells have valid clone size > 0. ",
                 "Minimum 10 required for correlation analysis."),
          type = "warning", duration = 8
        )
        return()
      }

      cor_results(NULL)
      shinyjs::hide("download_panel")

      withProgress(message = paste("Analysing", input$repertoire_type, "correlations..."), value = 0, {

        incProgress(0.03, detail = "Subsetting cells...")
        so_sub <- subset(so, cells = colnames(so)[valid_cells])
        clone_sizes_sub <- as.numeric(so_sub@meta.data[[col_name]])

        counts <- Seurat::GetAssayData(so_sub, assay = "RNA", layer = "data")

        incProgress(0.05, detail = "Filtering genes...")
        cells_per_gene <- Matrix::rowSums(counts > 0)
        genes_to_keep  <- names(which(cells_per_gene >= input$min_cells))
        counts         <- counts[genes_to_keep, , drop = FALSE]
        num_genes      <- nrow(counts)

        if (num_genes == 0) {
          showNotification(
            paste0("No genes passed the minimum cell filter (>= ", input$min_cells, " cells). ",
                   "Try lowering the threshold."),
            type = "warning", duration = 8
          )
          return()
        }

        method_str <- input$cor_method
        chunk_size <- 500
        n_chunks   <- ceiling(num_genes / chunk_size)
        res_chunks <- vector("list", n_chunks)

        incProgress(0.02, detail = paste0("Computing correlations for ", num_genes, " genes..."))

        for (ci in seq_len(n_chunks)) {
          idx   <- seq((ci - 1) * chunk_size + 1, min(ci * chunk_size, num_genes))
          chunk <- as.matrix(counts[idx, , drop = FALSE])

          res_chunks[[ci]] <- apply(chunk, 1, function(x) {
            if (stats::sd(x) == 0) return(c(NA_real_, NA_real_))
            ct <- suppressWarnings(
              stats::cor.test(x, clone_sizes_sub, method = method_str, exact = FALSE)
            )
            c(ct$estimate[[1]], ct$p.value)
          })

          incProgress(
            0.85 / n_chunks,
            detail = paste0("Chunk ", ci, "/", n_chunks,
                            " (", min(ci * chunk_size, num_genes), "/", num_genes, " genes)")
          )
        }

        res <- do.call(cbind, lapply(res_chunks, function(r) {
          if (is.matrix(r)) r else matrix(r, nrow = 2)
        }))

        cor_df <- data.frame(
          Gene        = rownames(counts),
          Correlation = res[1, ],
          P_value     = res[2, ]
        )
        cor_df        <- stats::na.omit(cor_df)
        cor_df$P_adj  <- stats::p.adjust(cor_df$P_value, method = "fdr")
        cor_df        <- dplyr::arrange(cor_df, dplyr::desc(Correlation))

        cor_results(cor_df)
        shinyjs::show("download_panel")

        n_sig <- sum(cor_df$P_adj <= input$p_val_adj, na.rm = TRUE)
        incProgress(0.05, detail = "Done.")
        showNotification(
          paste0("Completed! ", n_sig, " significant genes found (FDR < ", input$p_val_adj, ")."),
          type = "message", duration = 6
        )
      })
    })

    output$cor_table <- renderDT({
      req(cor_results())
      sig_df <- dplyr::filter(cor_results(), P_adj <= input$p_val_adj)
      datatable(
        sig_df,
        options  = list(pageLength = 15, scrollX = TRUE, searching = TRUE),
        rownames = FALSE
      ) |>
        DT::formatSignif(columns = c("Correlation", "P_value", "P_adj"), digits = 4)
    })

    output$download_table <- downloadHandler(
      filename = function() {
        paste0(input$repertoire_type, "_Phenotype_Correlation_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        openxlsx::write.xlsx(cor_results(), file)
      }
    )

    make_scatter_plot <- reactive({
      req(myReactives$seurat_object)

      use_custom <- isTRUE(input$use_custom_genes)
      if (use_custom) {
        raw_text <- input$custom_genes_text %||% ""
        top_genes <- unique(trimws(unlist(strsplit(raw_text, "[,\n\r]+"))))
        top_genes <- top_genes[nzchar(top_genes)]
      } else {
        req(cor_results())
        top_genes <- cor_results() |>
          dplyr::filter(P_adj <= input$p_val_adj) |>
          dplyr::slice_head(n = input$top_n) |>
          dplyr::pull(Gene)
      }

      if (length(top_genes) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = "No significant genes found.\nTry relaxing the P-value threshold.",
                     size = 5, hjust = 0.5) +
            theme_void()
        )
      }

      so       <- myReactives$seurat_object
      col_name <- paste0(input$repertoire_type, "_cloneSize")

      if (!col_name %in% names(so@meta.data)) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste0(col_name, " not found.\nPlease upload ", input$repertoire_type, " data first."),
                     size = 5, hjust = 0.5) +
            theme_void()
        )
      }

      valid_cells <- rownames(so@meta.data)[
        !is.na(so@meta.data[[col_name]]) & so@meta.data[[col_name]] > 0
      ]
      so_sub <- subset(so, cells = valid_cells)

      # Filter to genes available in the assay
      available_genes <- rownames(Seurat::GetAssayData(so_sub, assay = "RNA", layer = "data"))
      top_genes <- intersect(top_genes, available_genes)
      if (length(top_genes) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = "None of the specified genes found in the dataset.",
                     size = 5, hjust = 0.5) +
            theme_void()
        )
      }

      df_list <- lapply(top_genes, function(g) {
        expr <- as.numeric(Seurat::GetAssayData(so_sub, assay = "RNA", layer = "data")[g, ])
        data.frame(
          Expression = expr,
          CloneSize  = as.numeric(so_sub@meta.data[[col_name]]),
          Gene       = g
        )
      })
      plot_df <- dplyr::bind_rows(df_list)

      # Correlation labels (only available when not using custom genes)
      cor_data <- if (!use_custom && !is.null(cor_results())) {
        cor_results() |>
          dplyr::filter(Gene %in% top_genes) |>
          dplyr::mutate(
            label = paste0("r = ", round(Correlation, 3),
                           "\nFDR = ", formatC(P_adj, format = "e", digits = 2))
          )
      } else {
        data.frame(Gene = character(0), label = character(0))
      }
      plot_df <- dplyr::left_join(plot_df, cor_data[, c("Gene", "label")], by = "Gene")
      if (!"label" %in% names(plot_df)) plot_df$label <- NA_character_

      pt_size   <- input$point_size %||% 0.7
      base_font <- input$base_font_size %||% 13
      leg_pos   <- input$legend %||% "none"

      label_df <- dplyr::distinct(plot_df, Gene, label) |> dplyr::filter(!is.na(label))

      p <- ggplot(plot_df, aes(x = CloneSize, y = Expression)) +
        geom_point(alpha = 0.4, size = pt_size, color = "steelblue") +
        geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 0.8)

      if (nrow(label_df) > 0) {
        p <- p + geom_text(
          data = label_df,
          aes(label = label),
          x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
          size = base_font * 0.25, color = "black", inherit.aes = FALSE
        )
      }

      p <- p +
        facet_wrap(~ Gene, scales = "free_y") +
        theme_minimal(base_size = base_font) +
        labs(
          x     = paste(input$repertoire_type, "Clone Size (# cells with same clonotype)"),
          y     = "Normalized Expression",
          title = paste(input$repertoire_type,
                        "– Gene Expression vs Clone Size (Top Significant Genes)")
        ) +
        theme(
          strip.text      = element_text(face = "bold", size = base_font),
          legend.position = leg_pos
        )
    })

    output$scatter_plots <- renderPlot({
      make_scatter_plot()
    }, res = 96,
       width  = reactive(input$plot_width  %||% 700),
       height = reactive(input$plot_height %||% 600))

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(input$repertoire_type, "_Phenotype_Scatter_", Sys.Date(), ".pptx")
      },
      content = function(file) {
        p <- make_scatter_plot()
        save_plot_as_pptx(file, p,
               input$plot_width  %||% 700,
               input$plot_height %||% 600)
      }
    )
  })
}
