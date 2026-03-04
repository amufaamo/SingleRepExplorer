# BCR SHM and Isotype Analysis Module
# Analyzes Somatic Hypermutation (SHM) rates and Isotype (C-gene) distribution

library(shiny)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DT)

# --- UI Definition ---
shmIsotypeUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("BCR SHM & Isotype"),
      p("Analyze B-cell somatic hypermutation and class switching."),
      selectInput(ns("group_by"), "Group By:", choices = c("sample", "seurat_clusters", "BCR_pair_c_gene"), selected = "sample"),
      hr(),
      h5("SHM Options"),
      radioButtons(ns("shm_plot_type"), "SHM Plot Type:", choices = c("Violin" = "violin", "Boxplot" = "boxplot", "Histogram" = "hist")),
      hr(),
      h5("Isotype Options"),
      checkboxInput(ns("show_percentage"), "Show Percentage for Isotypes", value = TRUE),
      hr(),
      downloadButton(ns("download_shm_plot"), "Download SHM Plot"),
      downloadButton(ns("download_isotype_plot"), "Download Isotype Plot")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("SHM Rate",
                 br(),
                 plotOutput(ns("shm_plot"), height = "500px"),
                 br(),
                 tags$h4("SHM Summary Statistics"),
                 tableOutput(ns("shm_stats"))
        ),
        tabPanel("Isotype Distribution",
                 br(),
                 plotOutput(ns("isotype_plot"), height = "500px"),
                 br(),
                 tags$h4("Isotype Counts"),
                 DT::DTOutput(ns("isotype_table"))
        )
      )
    )
  )
}

# --- Server Logic ---
shmIsotypeServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive data for BCR analysis
    bcr_data <- reactive({
      req(myReactives$bcr_df)
      df <- myReactives$bcr_df
      
      # Ensure v_identity is numeric
      v_id_col <- "BCR_pair_v_identity"
      if (v_id_col %in% names(df)) {
        df[[v_id_col]] <- as.numeric(df[[v_id_col]])
        # Calculate SHM frequency (100 - %identity)
        df$shm_freq <- 100 - df[[v_id_col]]
      } else {
        # Fallback if v_identity is missing
        df$shm_freq <- NA_real_
      }
      return(df)
    })

    # --- SHM Plotting ---
    output$shm_plot <- renderPlot({
      df <- bcr_data()
      req(df, "shm_freq" %in% names(df))
      
      # Remove NAs for plotting
      df_plot <- df %>% filter(!is.na(shm_freq))
      if (nrow(df_plot) == 0) return(ggplot() + annotate("text", x=1, y=1, label="No SHM data available.") + theme_void())

      p <- ggplot(df_plot, aes(x = .data[[input$group_by]], y = shm_freq, fill = .data[[input$group_by]]))
      
      if (input$shm_plot_type == "violin") {
        p <- p + geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill = "white")
      } else if (input$shm_plot_type == "boxplot") {
        p <- p + geom_boxplot()
      } else {
        p <- ggplot(df_plot, aes(x = shm_freq, fill = .data[[input$group_by]])) + 
             geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
             facet_wrap(~ .data[[input$group_by]])
      }

      p + theme_pubr() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(y = "SHM Frequency (100 - %V-Identity)", x = input$group_by, title = "Somatic Hypermutation Rate")
    })

    output$shm_stats <- renderTable({
      df <- bcr_data()
      req(df, "shm_freq" %in% names(df))
      
      df %>% 
        filter(!is.na(shm_freq)) %>%
        group_by(.data[[input$group_by]]) %>%
        summarise(
          Mean = mean(shm_freq),
          Median = median(shm_freq),
          SD = sd(shm_freq),
          Count = n(),
          .groups = "drop"
        )
    })

    # --- Isotype Plotting ---
    output$isotype_plot <- renderPlot({
      df <- bcr_data()
      isotype_col <- "BCR_pair_c_gene"
      req(df, isotype_col %in% names(df))
      
      df_isotype <- df %>% 
        filter(!is.na(.data[[isotype_col]]) & .data[[isotype_col]] != "")
      
      if (nrow(df_isotype) == 0) return(ggplot() + annotate("text", x=1, y=1, label="No isotype data available.") + theme_void())

      p <- ggplot(df_isotype, aes(x = .data[[input$group_by]], fill = .data[[isotype_col]]))
      
      if (input$show_percentage) {
        p <- p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent)
      } else {
        p <- p + geom_bar(position = "stack")
      }

      p + theme_pubr() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = input$group_by, y = "Proportion", fill = "Isotype (C-gene)", title = "Isotype Distribution")
    })

    output$isotype_table <- DT::renderDT({
        df <- bcr_data()
        isotype_col <- "BCR_pair_c_gene"
        req(df, isotype_col %in% names(df))
        
        counts <- df %>%
            group_by(.data[[input$group_by]], .data[[isotype_col]]) %>%
            summarise(Count = n(), .groups = "drop")
            
        DT::datatable(counts, options = list(pageLength = 10))
    })

    # --- Downloads ---
    output$download_shm_plot <- downloadHandler(
        filename = function() { "shm_plot.pdf" },
        content = function(file) {
            # Recalculate p for download (simplified)
            df <- bcr_data()
            p <- ggplot(df, aes(x = .data[[input$group_by]], y = shm_freq, fill = .data[[input$group_by]])) + 
                 geom_boxplot() + theme_pubr()
            ggsave(file, plot = p, width = 8, height = 6)
        }
    )
  })
}
