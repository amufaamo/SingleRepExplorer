# Phylogenetic Tree Analysis Module for BCR Data
# Uses 'dowser' and 'ggtree' for lineage tree construction and visualization

library(shiny)
library(ggtree)
library(dowser)
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)

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
      hr(),
      downloadButton(ns("download_plot"), "Download Tree (PDF)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Phylogenetic Tree", 
                 br(),
                 textOutput(ns("plot_info")),
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

    # 1. Update clone choices
    observe({
      req(myReactives$bcr_df)
      df <- myReactives$bcr_df

      # Find clones with multiple subclotypes (suitable for tree)
      clones_summary <- tryCatch({
        df %>%
          dplyr::filter(!is.na(raw_clonotype_id) & raw_clonotype_id != "None") %>%
          dplyr::group_by(raw_clonotype_id) %>%
          dplyr::summarise(n_distinct_subclones = dplyr::n_distinct(barcode), .groups = 'drop') %>%
          dplyr::filter(n_distinct_subclones >= 2) %>%
          dplyr::pull(raw_clonotype_id) %>%
          sort()
      }, error = function(e) {
        return(character(0))
      })

      if (length(clones_summary) > 0) {
        updateSelectInput(session, "selected_clone", choices = clones_summary, selected = clones_summary[1])
      } else {
        updateSelectInput(session, "selected_clone", choices = c("No suitable clones found" = ""), selected = "")
      }
    })

    # 2. Extract selected clone data
    selected_clone_data <- eventReactive(input$run_phylotree, {
      req(myReactives$bcr_df, input$selected_clone, input$selected_clone != "")
      selected_clonotype_id <- input$selected_clone
      df <- myReactives$bcr_df

      clone_subset <- df %>%
        dplyr::filter(
          raw_clonotype_id == selected_clonotype_id,
          !is.na(BCR_IGH_full_length_nt) & BCR_IGH_full_length_nt != "",
          !is.na(barcode) & barcode != ""
        ) %>%
        dplyr::distinct(barcode, .keep_all = TRUE)

      if (nrow(clone_subset) < 2) {
        showNotification("Insufficient sequences (< 2) for this clone.", type = "warning")
        return(NULL)
      }
      return(clone_subset)
    })

    # 3. Format for dowser
    phylotree_input_data <- reactive({
      req(selected_clone_data())
      clone_subset <- selected_clone_data()

      # Filter out columns with all NAs to prevent formatting errors
      cols_to_keep <- names(clone_subset)[colSums(!is.na(clone_subset)) > 0]
      clone_subset <- clone_subset[, cols_to_keep]

      formatted_data <- tryCatch({
        dowser::formatClones(
          clone_subset,
          clone = "raw_clonotype_id",
          seq = "BCR_IGH_full_length_nt",
          id = 'barcode'
        )
      }, error = function(e) {
        showNotification(paste("Formatting error:", e$message), type = "error")
        NULL
      })
      return(formatted_data)
    })

    # 4. Build tree
    phylogenetic_tree <- reactive({
      req(phylotree_input_data())
      formatted_clones <- phylotree_input_data()

      withProgress(message = "Building tree...", {
        tree_list <- tryCatch({
          dowser::getTrees(formatted_clones, build = "nj")
        }, error = function(e) {
          showNotification(paste("Tree building error:", e$message), type = "error")
          NULL
        })
      })
      
      req(tree_list, length(tree_list) > 0)
      return(tree_list[[1]]) # Return the first treedata object
    })

    # 5. Build plot
    tree_plot_object <- reactive({
      req(phylogenetic_tree())
      tree_data <- phylogenetic_tree()
      
      color_col <- input$color_by
      label_col <- input$label_by
      
      # Tip data check
      tip_data <- as_tibble(tree_data) %>% dplyr::filter(isTip)
      
      p <- ggtree(tree_data) +
           geom_tree(linewidth = 0.8) +
           geom_treescale(x = 0, y = -1, fontsize = 3)

      if (color_col != "None" && color_col %in% names(tip_data)) {
        p <- p + geom_tippoint(aes(color = .data[[color_col]]), size = 3)
      } else {
        p <- p + geom_tippoint(size = 3, color = "steelblue")
      }

      if (label_col %in% names(tip_data)) {
        p <- p + geom_tiplab(aes(label = .data[[label_col]]), 
                            size = input$label_size, 
                            offset = input$label_offset, 
                            align = TRUE)
      }
      
      p <- p + theme(legend.position = "right") +
           labs(title = paste("Clonotype:", input$selected_clone)) +
           xlim_tree(max(p$data$x, na.rm=TRUE) * 1.5)
      
      return(p)
    })

    # 6. Outputs
    output$plot_info <- renderText({
      req(phylogenetic_tree())
      tree_obj <- phylogenetic_tree()
      paste("Tree built for clone:", input$selected_clone, "| Sequences:", Ntip(tree_obj))
    })

    output$plot <- renderPlot({
      req(tree_plot_object())
      tree_plot_object()
    }, height = function() {
      tree_obj <- phylogenetic_tree()
      if (!is.null(tree_obj)) max(500, Ntip(tree_obj) * 20) else 600
    })

    output$table <- DT::renderDT({
      req(selected_clone_data())
      DT::datatable(selected_clone_data(), options = list(pageLength = 10, scrollX = TRUE))
    })

    output$download_plot <- downloadHandler(
      filename = function() { paste0("tree_", input$selected_clone, ".pdf") },
      content = function(file) {
        tree_obj <- phylogenetic_tree()
        h <- if (!is.null(tree_obj)) max(6, Ntip(tree_obj) * 0.3) else 8
        ggsave(file, plot = tree_plot_object(), width = 10, height = h)
      }
    )
  })
}