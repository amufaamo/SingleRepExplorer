library(shiny)
library(Seurat)
library(DT)
library(dplyr)
library(plotly)
library(ggplot2)
library(ggplot2)
library(shinyjs)
library(enrichR)



differentialGeneExpressionUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      selectInput(
        ns("analysis_type"),
        label = "Comparison Type",
        choices = c(
          "To entire dataset" = "all",
          "Between selected group(s)" = "two",
          "Multi-group (LRT: ANOVA-like)" = "multi"
        ),
        selected = "all"
      ),
      selectInput(ns("group_by"), "Group by", choices = c("sample", "seurat_clusters"), selected = "sample"),
      conditionalPanel(
        condition = "input.analysis_type == 'two'",
        ns = ns,
        selectInput(ns("target_cluster"), label = "Target cluster(s)", choices = NULL, selected = NULL, multiple = TRUE),
        selectInput(ns("reference_cluster"), label = "Reference cluster(s)", choices = NULL, selected = NULL, multiple = TRUE),
      ),
      conditionalPanel(
        condition = "input.analysis_type == 'multi'",
        ns = ns,
        selectInput(ns("multi_clusters"), label = "Target clusters (requires 3+)", choices = NULL, selected = NULL, multiple = TRUE)
      ),
      numericInput(ns("logfc"), "Threshold of Log fold change (default: 0.1)", min = 0, max = 1, value = 0.1, step = 0.01),
      numericInput(ns("minpct"), "Minimum percentage of cells expressing a gene (%) (default: 0.01)", min = 0, max = 1, value = 0.01, step = 0.01),
      conditionalPanel(
        condition = "input.analysis_type == 'all'",
        ns = ns,
        selectInput(
          ns("choice"),
          "Show Features",
          choices = list(
            "Only Up-regulated" = "positive",
            "Only Down-regulated" = "negative",
            "Both" = "both"
          ),
          selected = "positive"
        ),
      ),
      actionButton(ns("run"), "Run"),
      conditionalPanel(
        condition = "input.analysis_type == 'all' && output.table_markers",
        ns = ns,
        numericInput(ns("num"), label = "Number of markers", value = 10),
        numericInput(ns("p_val_adj_marker"), "Threshold of adjusted p-value (default: 0.05)", min = 0, max = 1, value = 0.05, step = 0.01),
      ),
      conditionalPanel(
        condition = "input.analysis_type == 'two' && output.table_two",
        ns = ns,
        textInput(ns("gene_input"), "Highlight Genes (comma separated)", value = ""),
        numericInput(ns("p_val_adj"), "Threshold of adjusted p-value (default: 0.05)", min = 0, max = 1, value = 0.05, step = 0.01),
        numericInput(ns("highlight_size"), "Highlight Size", min = 0.001, max = 10, value = 2),
        numericInput(ns("other_size"), "Other Points Size", min = 0.001, max = 10, value = 1),
        numericInput(ns("point_alpha"), "Point Alpha", min = 0.1, max = 1, value = 0.7, step = 0.1),
        numericInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 400, step = 100),
        numericInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 400, step = 100),
        hr(),
        hr(),
        radioButtons(ns("enrich_method"), "Enrichment Method",
                     choices = c("EnrichR (ORA)" = "enrichr", "GSEA (fgsea)" = "gsea"),
                     selected = "enrichr"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'enrichr'", ns("enrich_method")),
          selectInput(ns("enrich_db"), "Pathway Database (Enrichr)", 
                      choices = c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "KEGG_2021_Human", "Reactome_2022"),
                      selected = "GO_Biological_Process_2023")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'gsea'", ns("enrich_method")),
          selectInput(ns("gsea_category"), "GSEA Category (MSigDB)",
                      choices = c("H (Hallmark)" = "H", "C2 (Curated/Pathways)" = "C2", "C5 (Gene Ontology)" = "C5"),
                      selected = "H")
        )
      ),
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.analysis_type == 'all'",
        ns = ns,
        h3("Most Differentially Expressed Features (markers) for Each Group"),
        downloadButton(ns("download_marker_table"), "Download marker table (.csv)"),
        DTOutput(ns("table_markers")),
        hr(),
        h3("Top Markers Expression (Dot Plot)"),
        p("Displays the expression of the top markers across all clusters intuitively."),
        downloadButton(ns("download_dotplot"), "Download Dot Plot (.pdf)"),
        plotOutput(ns("dotPlot"), height = "500px"),
        hr(),
        h3("Differentially Expressed Features across Groups"),
        downloadButton(ns("download_table_all"), "Download all table (.csv)"),
        DTOutput(ns("table_all")),
      ),
      conditionalPanel(
        condition = "input.analysis_type == 'two'",
        ns = ns,
        h3("Volcano Plot"),
        downloadButton(ns("download_plot"), "Download Plot (PDF)"),
        plotlyOutput(ns("volcanoPlot")),
        hr(),
        h3("Pathway Enrichment Analysis"),
        p("Automated gene pathway enrichment analysis. Choose EnrichR for ORA on top markers, or GSEA to evaluate the entire ranked gene list."),
        downloadButton(ns("download_enrichment_plot"), "Download Enrichment Plot (.pdf)"),
        plotOutput(ns("enrichmentPlot"), height = "400px"),
        hr(),
        h3("Differentially Expressed Features Between Selected Group(s)"),
        downloadButton(ns("download_table_two"), "Download all table (.csv)"), # Keep using all
        DTOutput(ns("table_two")),
      ),
      conditionalPanel(
        condition = "input.analysis_type == 'multi'",
        ns = ns,
        h3("Multi-group Variance (Kruskal-Wallis)"),
        p("Genes that significantly vary across 3 or more clusters (ANOVA-like)."),
        downloadButton(ns("download_table_multi"), "Download multi table (.csv)"),
        DTOutput(ns("table_multi"))
      )
    )
  )
}



differentialGeneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    # Initialize disabled state
    shinyjs::disable("download_marker_table")
    shinyjs::disable("download_table_all")
    shinyjs::disable("download_table_two")
    shinyjs::disable("download_table_multi")


    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_for_marker(session, input, myReactives)
    })


    # twoのときに、クラスターの種類
    observe({
      req(myReactives$seurat_object, input$group_by, input$analysis_type == "two")
      all_clusters <- unique(myReactives$seurat_object@meta.data[[input$group_by]])

      updated_target_choices <- setdiff(all_clusters, input$reference_cluster)
      numeric_choices_target <- as.numeric(updated_target_choices)
      sorted_choices_target <- if (all(!is.na(numeric_choices_target))) {
        sort(numeric_choices_target)
      } else {
        sort(updated_target_choices)
      }
      updateSelectInput(session, "target_cluster", choices = sorted_choices_target, selected = input$target_cluster)

      updated_reference_choices <- setdiff(all_clusters, input$target_cluster)
      numeric_choices_reference <- as.numeric(updated_reference_choices)
      sorted_choices_reference <- if (all(!is.na(numeric_choices_reference))) {
        sort(numeric_choices_reference)
      } else {
        sort(updated_reference_choices)
      }
      updateSelectInput(session, "reference_cluster", choices = sorted_choices_reference, selected = input$reference_cluster)
    })
    
    # multiのときに、クラスターの種類
    observe({
      req(myReactives$seurat_object, input$group_by, input$analysis_type == "multi")
      all_clusters <- unique(myReactives$seurat_object@meta.data[[input$group_by]])
      numeric_choices <- as.numeric(all_clusters)
      sorted_choices <- if (all(!is.na(numeric_choices))) sort(numeric_choices) else sort(all_clusters)
      updateSelectInput(session, "multi_clusters", choices = sorted_choices, selected = input$multi_clusters)
    })


    observeEvent(input$run, {
      req(myReactives$seurat_object, input$group_by)
      so <- myReactives$seurat_object
      Idents(so) <- input$group_by

      if (input$analysis_type == "all") {
        withProgress(message = "Calculating...", value = 0, {
          clusters <- unique(Idents(so))
          n_clusters <- length(clusters)
          all_data_list <- list()

          for (i in seq_along(clusters)) {
            cluster <- clusters[i]
            incProgress(1 / n_clusters, detail = paste("Processing cluster", cluster))

            # Use FindMarkers with only.pos based on input$choice
            markers <- FindMarkers(
              so,
              ident.1 = cluster,
              logfc.threshold = input$logfc,
              min.pct = input$minpct,
              only.pos = (input$choice == "positive"), # Correct only.pos handling
              verbose = FALSE
            )
            if (input$choice == "negative") {
              markers <- markers[markers$avg_log2FC < 0, ]
            } else if (input$choice == "both") {
              # No filter
            } else {
              markers <- markers[markers$avg_log2FC > 0, ]
            }
            markers$cluster <- cluster
            markers$gene <- rownames(markers)
            all_data_list[[i]] <- markers
          }
          myReactives$all_data_all <- do.call(rbind, all_data_list)
        })
      } else if (input$analysis_type == "two") {
        req(input$target_cluster, input$reference_cluster)
        withProgress(message = "Calculating...", value = 0, {
          myReactives$all_data_two <- FindMarkers(
            so,
            ident.1 = input$target_cluster,
            ident.2 = input$reference_cluster,
            logfc.threshold = input$logfc,
            min.pct = input$minpct
          )
          myReactives$all_data_two$gene <- rownames(myReactives$all_data_two) # Add gene names

          incProgress(1) # Full progress for single comparison
        })
      } else if (input$analysis_type == "multi") {
        req(length(input$multi_clusters) >= 3)
        withProgress(message = "Calculating Kruskal-Wallis...", value = 0, {
          so_sub <- subset(so, idents = input$multi_clusters)
          counts <- Seurat::GetAssayData(so_sub, assay = "RNA", layer = "data")
          groups <- as.character(Idents(so_sub))
          
          incProgress(0.1, detail = "Filtering genes...")
          cells_per_gene <- Matrix::rowSums(counts > 0)
          min_cells <- max(length(groups) * input$minpct, 3)
          genes_to_keep <- names(which(cells_per_gene >= min_cells))
          counts_matrix <- as.matrix(counts[genes_to_keep, , drop = FALSE])
          
          incProgress(0.4, detail = "Running statistical tests (may take 1-2 minutes)...")
          if (nrow(counts_matrix) > 0) {
              res <- suppressWarnings(apply(counts_matrix, 1, function(x) {
                  if(stats::sd(x) == 0) return(NA)
                  tryCatch({ stats::kruskal.test(x ~ groups)$p.value }, error = function(...) NA)
              }))
              pval_df <- data.frame(gene = names(res), P_value = as.numeric(res))
              pval_df <- na.omit(pval_df)
              pval_df$p_val_adj <- stats::p.adjust(pval_df$P_value, method = "BH")
              pval_df <- pval_df %>% dplyr::arrange(P_value)
              myReactives$all_data_multi <- pval_df
          } else {
              myReactives$all_data_multi <- data.frame(gene=character(0), P_value=numeric(0), p_val_adj=numeric(0))
          }
          incProgress(1)
        })
      }
    })


    # Marker Table (analysis_type == "all")
    observe({
      req(input$analysis_type == "all")
      req(myReactives$all_data_all)

      myReactives$marker_table <- myReactives$all_data_all %>%
        dplyr::filter(p_val_adj <= input$p_val_adj_marker) %>%
        group_by(cluster) %>%
        arrange(desc(abs(avg_log2FC))) %>%
        slice_head(n = input$num) %>%
        summarise(genes = paste(gene, collapse = ", "), .groups = "drop")
    })

    output$table_markers <- renderDT({
      req(myReactives$marker_table)
      datatable(
        myReactives$marker_table,
        options = list(paging = FALSE)
      )
    })

    # All Data Table (both analysis types)
    output$table_all <- renderDT({
      req(myReactives$all_data_all) # Use all_data
      datatable(myReactives$all_data_all,
        filter = "top",
        options = list(pageLength = 10)
      )
    })

    # All Data Table (both analysis types)
    output$table_two <- renderDT({
      req(myReactives$all_data_two) # Use all_data
      datatable(myReactives$all_data_two,
        filter = "top",
        options = list(pageLength = 10)
      )
    })
    
    # Enable/disable download buttons based on table content
    observe({
      if (input$analysis_type == "all") {
        if (is.null(myReactives$marker_table) || nrow(myReactives$marker_table) == 0) {
          shinyjs::disable("download_marker_table")
        } else {
          shinyjs::enable("download_marker_table")
        }
        if (is.null(myReactives$all_data_all) || nrow(myReactives$all_data_all) == 0) {
          shinyjs::disable("download_table_all")
        } else {
          shinyjs::enable("download_table_all")
        }
      } else if (input$analysis_type == "two") {
        if (is.null(myReactives$all_data_two) || nrow(myReactives$all_data_two) == 0) {
          shinyjs::disable("download_table_two")
        } else {
          shinyjs::enable("download_table_two")
        }
      } else if (input$analysis_type == "multi") {
         if (is.null(myReactives$all_data_multi) || nrow(myReactives$all_data_multi) == 0) {
             shinyjs::disable("download_table_multi")
         } else {
             shinyjs::enable("download_table_multi")
         }
      }
    })

    # All Data Table (multi-group)
    output$table_multi <- renderDT({
      req(myReactives$all_data_multi)
      datatable(myReactives$all_data_multi,
        filter = "top",
        options = list(pageLength = 10)
      )
    })
    
    output$download_table_multi <- downloadHandler(
      filename = function() { "multi_varied_genes.csv" },
      content = function(file) { write.csv(myReactives$all_data_multi, file, row.names = FALSE) }
    )

    volcano_plot_data <- reactive({
      req(myReactives$all_data_two)

      filtered_data <- myReactives$all_data_two %>%
        dplyr::filter(p_val_adj <= input$p_val_adj)

      volcano_data <- filtered_data %>%
        mutate(log_p_value = -log10(p_val_adj))

      highlighted <- if (!is.null(input$gene_input) && nzchar(input$gene_input)) {
        genes <- strsplit(input$gene_input, ",\\s*")[[1]]
        tolower(volcano_data$gene) %in% tolower(genes) # Use the gene column
      } else {
        rep(FALSE, nrow(volcano_data))
      }

      num_highlighted <- sum(highlighted)

      p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log_p_value, text = gene)) + # Use the gene column
        geom_point(
          aes(
            color = ifelse(highlighted, "Highlighted", "Other"),
            size = ifelse(highlighted, input$highlight_size, input$other_size)
          ),
          alpha = 0.7
        ) +
        scale_color_manual(
          values = c("Highlighted" = "red", "Other" = "gray"),
          labels = c("Highlighted" = paste0("Highlighted (", num_highlighted, ")"), "Other" = "Other")
        ) +
        scale_size_identity() +
        geom_text(
          data = subset(volcano_data, highlighted),
          aes(label = gene), # Use the gene column
          size = 3,
          vjust = -0.5,
          hjust = -0.1,
          check_overlap = TRUE
        ) +
        labs(
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted p-value"
        ) +
        theme_classic() +
        theme(legend.position = "none")

      return(p)
    })

    output$volcanoPlot <- renderPlotly({
      volcano_plot_data() %>%
        ggplotly(tooltip = "text", width = input$plot_width, height = input$plot_height) %>%
        layout(
          plot_bgcolor = "white",
          panel_grid = list(visible = FALSE),
          xaxis = list(showline = TRUE, linecolor = "black"),
          yaxis = list(showline = TRUE, linecolor = "black"),
          showlegend = FALSE
        )
    })
    
    # Dot Plot implementation for 'all' mode
    dotplot_reactive <- reactive({
      req(input$analysis_type == "all", myReactives$marker_table, myReactives$all_data_all)
      so <- myReactives$seurat_object
      Idents(so) <- input$group_by
      
      # Extract list of top genes dynamically
      top_n <- min(input$num, 5) # Cap at 5 per cluster to prevent overcrowding the plot
      top_markers <- myReactives$all_data_all %>%
        dplyr::filter(p_val_adj <= input$p_val_adj_marker) %>%
        group_by(cluster) %>%
        arrange(desc(abs(avg_log2FC))) %>%
        slice_head(n = top_n)
        
      genes_to_plot <- unique(top_markers$gene)
      req(length(genes_to_plot) > 0)
      
      p <- Seurat::DotPlot(so, features = genes_to_plot) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10),
              axis.text.y = element_text(size=12)) +
        labs(x = "Top Marker Genes", y = "Groups", title = paste("Top Marker Expression Across", tools::toTitleCase(input$group_by)))
      return(p)
    })
    
    output$dotPlot <- renderPlot({
      dotplot_reactive()
    })
    
    # Automated Pathway Enrichment for 'two' mode
    enrichment_reactive <- reactive({
      req(input$analysis_type == "two", myReactives$all_data_two)
      
      if (input$enrich_method == "enrichr") {
          # --- EnrichR (ORA) Logic ---
          # Filter to highly significant and upregulated relative genes
          sig_genes <- myReactives$all_data_two %>%
            dplyr::filter(p_val_adj <= input$p_val_adj, avg_log2FC > 0) %>%
            pull(gene)
            
          if (length(sig_genes) == 0) {
            return(ggplot() + annotate("text", x = 1, y = 1, label = "No significant Up-regulated genes found to enrich.") + theme_void())
          }
          
          dbs <- input$enrich_db
          
          tryCatch({
            enriched <- enrichR::enrichr(sig_genes, dbs)
            res <- enriched[[dbs]]
            if (nrow(res) == 0) { return(ggplot() + annotate("text", x = 1, y = 1, label = "No enriched pathways found.") + theme_void()) }
            
            top_res <- res %>% arrange(P.value) %>% head(10) %>%
              mutate(Term = stringr::str_wrap(stringr::str_remove(Term, " \\(GO:.*\\)"), width = 50)) %>% 
              mutate(Term = factor(Term, levels = rev(Term)))
              
            p <- ggplot(top_res, aes(x = -log10(P.value), y = Term)) +
              geom_col(fill = "steelblue") + theme_minimal(base_size = 14) +
              labs(x = "-Log10(P-value)", y = "", title = paste("Enriched Pathways (", dbs, ")", sep=""),
                   subtitle = paste("Based on", length(sig_genes), "significantly up-regulated genes")) +
              theme(axis.text.y = element_text(size = 11, color="black"), plot.title = element_text(face="bold"))
            return(p)
          }, error = function(e) { return(ggplot() + annotate("text", x = 1, y = 1, label = paste("EnrichR API Error:", e$message)) + theme_void()) })
          
      } else if (input$enrich_method == "gsea") {
          # --- GSEA Logic using fgsea and msigdbr ---
          if (!requireNamespace("fgsea", quietly = TRUE) || !requireNamespace("msigdbr", quietly = TRUE)) {
             return(ggplot() + annotate("text", x=1, y=1, label="Error: 'fgsea' or 'msigdbr' packages not installed.\nPlease run update/build docker image.") + theme_void())
          }
          tryCatch({
              # Prepare ranks
              ranked_df <- myReactives$all_data_two %>% filter(!is.na(avg_log2FC)) %>% arrange(desc(avg_log2FC))
              ranks <- setNames(ranked_df$avg_log2FC, ranked_df$gene)
              
              if(length(ranks) == 0) return(ggplot() + annotate("text", x=1,y=1,label="No genes available for ranking.") + theme_void())
              
              # Prepare pathways
              md <- msigdbr::msigdbr(species = "Homo sapiens", category = input$gsea_category)
              pathways <- split(x = md$gene_symbol, f = md$gs_name)
              
              # Run fgsea
              # suppress warnings about ties
              suppressWarnings({
                  fgseaRes <- fgsea::fgsea(pathways = pathways, stats = ranks, minSize=15, maxSize=500, nproc=1)
              })
              
              if(nrow(fgseaRes) == 0) return(ggplot() + annotate("text", x=1,y=1, label="No pathways enriched.") + theme_void())
              
              top_pathways <- fgseaRes %>% filter(padj < 0.25) %>% arrange(padj) %>% head(15) %>%
                  mutate(Pathway = stringr::str_replace_all(pathway, "_", " ")) %>%
                  mutate(Pathway = stringr::str_remove(Pathway, "^(HALLMARK |GO |KEGG |REACTOME )")) %>%
                  mutate(Pathway = stringr::str_wrap(Pathway, width=45)) %>%
                  mutate(Pathway = factor(Pathway, levels = rev(Pathway))) %>%
                  mutate(Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"))
              
              if(nrow(top_pathways) == 0) return(ggplot() + annotate("text", x=1,y=1, label="No significant pathways found (FDR < 0.25).") + theme_void())
              
              p <- ggplot(top_pathways, aes(x = NES, y = Pathway, fill = padj)) +
                  geom_col() +
                  scale_fill_viridis_c(direction = -1, name="FDR p-adj") +
                  theme_minimal(base_size = 14) +
                  labs(x = "Normalized Enrichment Score (NES)", y = "", 
                       title = paste("GSEA Top Pathways (", input$gsea_category, ")"),
                       subtitle = "Ranked by avg_log2FC (All Genes)") +
                  theme(axis.text.y = element_text(size = 10, color="black"), plot.title = element_text(face="bold"))
                  
              return(p)
          }, error = function(e) { return(ggplot() + annotate("text", x=1, y=1, label=paste("GSEA Error:", e$message)) + theme_void()) })
      }
    })

    output$enrichmentPlot <- renderPlot({
      enrichment_reactive()
    })

    output$download_marker_table <- downloadHandler(
      filename = function() {
        "marker_table.csv"
      },
      content = function(file) {
        write.csv(myReactives$marker_table, file, row.names = FALSE)
      }
    )

    output$download_table_all <- downloadHandler(
      filename = function() {
        "all_table.csv"
      },
      content = function(file) {
        write.csv(myReactives$all_data_all, file, row.names = TRUE) # Use all_data
      }
    )

    output$download_table_two <- downloadHandler(
      filename = function() {
        "all_table.csv"
      },
      content = function(file) {
        write.csv(myReactives$all_data_two, file, row.names = TRUE) # Use all_data
      }
    )

   # PDF Download Handler
    output$download_plot <- downloadHandler(
      filename = function() {
        "Volcano_plot.pdf"
      },
      content = function(file) {
        ggsave(file, plot = volcano_plot_data(), device = "pdf", width = input$plot_width / 72, height = input$plot_height / 72, units = "in", dpi = 300)
      }
    )
    
    # Download Handlers for new plots
    output$download_dotplot <- downloadHandler(
      filename = function() { "Marker_DotPlot.pdf" },
      content = function(file) {
        ggsave(file, plot = dotplot_reactive(), device = "pdf", width = 14, height = 7, units = "in")
      }
    )
    
    output$download_enrichment_plot <- downloadHandler(
      filename = function() { "Pathway_Enrichment_Plot.pdf" },
      content = function(file) {
        ggsave(file, plot = enrichment_reactive(), device = "pdf", width = 10, height = 6, units = "in")
      }
    )


    # Enable/disable download buttons based on table content
    # Move this to top observation area, already done above!
  })
}
