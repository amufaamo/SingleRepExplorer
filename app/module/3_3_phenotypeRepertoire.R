library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(shinyjs)

phenotypeRepertoireUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      useShinyjs(),
      h4("Phenotype-Repertoire Correlation Engine"),
      p("Identify genes whose expression directly correlates with clonotype expansion sizes."),
      selectInput(ns("repertoire_type"), "Repertoire Type", choices = c("TCR", "BCR"), selected = "TCR"),
      selectInput(ns("cor_method"), "Correlation Method", choices = c("spearman", "pearson"), selected = "spearman"),
      numericInput(ns("min_cells"), "Minimum Cells per Gene", value = 10, min = 1),
      numericInput(ns("p_val_adj"), "Adjusted P-val threshold", value = 0.05, min = 0, step = 0.01),
      numericInput(ns("top_n"), "Plots for Top N genes", value = 6, min = 1, max = 20),
      actionButton(ns("run_cor"), "Run Correlation Engine", class = "btn-primary", width = "100%", icon = icon("calculator"))
    ),
    mainPanel(
      h3("Gene Correlation with Clone Size"),
      p("Genes positively correlated with Clone Size represent markers of clonal expansion (e.g., exhaustion, effector state)."),
      downloadButton(ns("download_table"), "Download Results (.csv)"),
      DTOutput(ns("cor_table")),
      hr(),
      h3("Expression vs Clone Size (Top Genes)"),
      plotOutput(ns("scatter_plots"), height = "600px")
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
        showNotification(paste(col_name, "not found. Ensure", input$repertoire_type, "data was uploaded properly."), type = "error")
        return()
      }
      
      clone_sizes <- as.numeric(so@meta.data[[col_name]])
      valid_cells <- which(!is.na(clone_sizes) & clone_sizes > 0)
      
      if (length(valid_cells) < 10) {
        showNotification("Not enough cells with valid repertoire data for correlation.", type="warning")
        return()
      }
      
      withProgress(message = paste("Analysing", input$repertoire_type, "correlations..."), value = 0, {
        # Subset to valid cells
        so_sub <- subset(so, cells = colnames(so)[valid_cells])
        clone_sizes_sub <- as.numeric(so_sub@meta.data[[col_name]])
        
        counts <- Seurat::GetAssayData(so_sub, assay = "RNA", layer = "data")
        
        # Filter genes
        incProgress(0.1, detail = "Filtering genes...")
        cells_per_gene <- Matrix::rowSums(counts > 0)
        genes_to_keep <- names(which(cells_per_gene >= input$min_cells))
        counts <- counts[genes_to_keep, , drop = FALSE]
        
        # Calculate correlation for each gene
        incProgress(0.3, detail = "Computing correlations (this might take a moment)...")
        counts_matrix <- as.matrix(counts)
        num_genes <- nrow(counts_matrix)
        
        if (num_genes == 0) {
             showNotification("No genes passed the minimum cell expression filter.", type="warning")
             return()
        }
        
        method_str <- input$cor_method
        
        # We can optimize via fast correlator but using apply is acceptable for a few thousand genes
        res <- suppressWarnings(apply(counts_matrix, 1, function(x) {
            if (stats::sd(x) == 0) return(c(NA, NA))
            c <- stats::cor.test(x, clone_sizes_sub, method = method_str, exact = FALSE)
            c(c$estimate, c$p.value)
        }))
        
        cor_df <- data.frame(
            Gene = rownames(counts_matrix),
            Correlation = res[1,],
            P_value = res[2,]
        )
        cor_df <- na.omit(cor_df)
        cor_df$P_adj <- stats::p.adjust(cor_df$P_value, method = "fdr")
        cor_df <- cor_df %>% dplyr::arrange(dplyr::desc(Correlation))
        
        cor_results(cor_df)
        incProgress(1)
      })
    })
    
    output$cor_table <- renderDT({
      req(cor_results())
      datatable(cor_results() %>% filter(P_adj <= input$p_val_adj), options = list(pageLength = 10))
    })
    
    output$download_table <- downloadHandler(
      filename = function() { paste0(input$repertoire_type, "_Phenotype_Correlation.csv") },
      content = function(file) { write.csv(cor_results(), file, row.names = FALSE) }
    )
    
    output$scatter_plots <- renderPlot({
      req(cor_results(), myReactives$seurat_object)
      top_genes <- cor_results() %>% filter(P_adj <= input$p_val_adj) %>% head(input$top_n) %>% pull(Gene)
      if (length(top_genes) == 0) return(ggplot() + annotate("text", x=1,y=1,label="No significant genes found.") + theme_void())
      
      so <- myReactives$seurat_object
      col_name <- paste0(input$repertoire_type, "_cloneSize")
      
      # For visualization limit to cells with repertoire info
      valid_cells <- rownames(so@meta.data)[!is.na(so@meta.data[[col_name]]) & so@meta.data[[col_name]] > 0]
      so_sub <- subset(so, cells = valid_cells)
      
      df_list <- lapply(top_genes, function(g) {
         expr <- Seurat::GetAssayData(so_sub, assay="RNA", layer="data")[g, ]
         data.frame(Expression = expr, CloneSize = factor(so_sub@meta.data[[col_name]]), Gene = g)
      })
      plot_df <- dplyr::bind_rows(df_list)
      
      ggplot(plot_df, aes(x = CloneSize, y = Expression, fill = CloneSize)) +
        geom_violin(scale = "width", alpha = 0.5) +
        geom_jitter(width=0.2, alpha=0.5, size=0.5) +
        geom_smooth(aes(group=1), method="lm", color="red", se = FALSE) +
        facet_wrap(~ Gene, scales = "free_y") +
        theme_minimal(base_size = 14) +
        labs(x = "Clone Size", y = "Normalized Expression") +
        theme(legend.position = "none", strip.text = element_text(face="bold", size=12))
    })
  })
}
