# BCR SHM and Isotype Analysis Module

# --- UI Definition ---
shmIsotypeUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("BCR SHM & Isotype"),
      p("Analyze B-cell somatic hypermutation and class switching."),
      groupByInput(ns),
      hr(),
      h5("SHM Options"),
      radioButtons(ns("shm_plot_type"), "SHM Plot Type:", choices = c("Violin" = "violin", "Boxplot" = "boxplot", "Histogram" = "hist")),
      jitterPointsInput(ns),
      hr(),
      h5("Isotype Options"),
      checkboxInput(ns("show_percentage"), "Show Percentage for Isotypes", value = TRUE),
      hr(),
      h5("Plot Options"),
      commonPlotOptions(ns, legend_selected = "right", width_value = 700, height_value = 500),
      hr(),
      # downloadButtons moved to mainPanel
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("SHM Rate",
                 br(),
                 downloadButton(ns("download_shm_plot"), "Download Plot (.pptx)"),
                 br(), br(),
                 plotOutput(ns("shm_plot")),
                 br(),
                 tags$h4("SHM Summary Statistics"),
                 tableOutput(ns("shm_stats"))
        ),
        tabPanel("Isotype Distribution",
                 br(),
                 downloadButton(ns("download_isotype_plot"), "Download Plot (.pptx)"),
                 br(), br(),
                 plotOutput(ns("isotype_plot")),
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

    observeEvent(myReactives$seurat_object, {
       update_group_by_select_input(session, myReactives)
    })

    bcr_data <- reactive({
      req(myReactives$bcr_df)
      df <- myReactives$bcr_df

      # Strip Rle / S4 / list columns to plain vectors to avoid ggplot2/dplyr crashes
      df <- as.data.frame(lapply(df, function(col) {
        if (inherits(col, "Rle") || inherits(col, "List")) {
          vec <- tryCatch(as.vector(col), error = function(e) as.character(col))
          if (is.list(vec))
            vapply(vec, function(x) if (length(x) == 0) NA_character_ else as.character(x[[1]]), NA_character_)
          else vec
        } else if (is.list(col)) {
          vapply(col, function(x) if (length(x) == 0) NA_character_ else as.character(x[[1]]), NA_character_)
        } else {
          col
        }
      }), stringsAsFactors = FALSE)

      # Ensure group_by column is in df
      group_col <- input$group_by %||% "sample"
      if (!group_col %in% names(df) && !is.null(myReactives$seurat_object)) {
         if (group_col %in% names(myReactives$seurat_object@meta.data)) {
            if ("barcode" %in% names(df)) {
              meta <- myReactives$seurat_object@meta.data[, c(group_col), drop=FALSE]
              meta$barcode <- rownames(meta)
              df <- df %>% left_join(meta, by="barcode")
            } else {
              df[[group_col]] <- "Unknown"
            }
         } else {
            df[[group_col]] <- "Unknown"
         }
      }

      # Pick a v_identity column. Preferred: BCR_pair_v_identity (from combineBCR + restore).
      # Fallback: BCR_IGH_v_identity (always populated directly from the raw CSV).
      # Last resort: any column ending with "v_identity" (handles .x/.y suffix cases from
      # older qs2 sessions saved before the bcr_csv_to_dataframe join-conflict fix).
      v_id_candidates <- c("BCR_pair_v_identity", "BCR_IGH_v_identity")
      v_id_col <- v_id_candidates[v_id_candidates %in% names(df)][1]
      if (is.na(v_id_col)) {
        suffixed <- grep("v_identity(\\.x)?$", names(df), value = TRUE)
        if (length(suffixed) > 0) v_id_col <- suffixed[1]
      }

      if (!is.na(v_id_col) && !is.null(v_id_col)) {
        vals <- suppressWarnings(as.numeric(as.character(df[[v_id_col]])))
        df$shm_freq <- 100 - vals
        message(sprintf("[shmIsotype] Using '%s' for SHM calculation (%d non-NA values).",
                        v_id_col, sum(!is.na(vals))))
      } else {
        df$shm_freq <- NA_real_
        message(sprintf("[shmIsotype] No v_identity column found. Available BCR columns: %s",
                        paste(grep("^BCR_", names(df), value = TRUE), collapse = ", ")))
      }
      return(df)
    })

    build_shm_plot <- function() {
      df <- bcr_data()
      validate(need(!is.null(input$group_by) && nzchar(input$group_by), "Select a grouping variable."))
      validate(need("shm_freq" %in% names(df), "SHM data not found in the repertoire table."))

      df_plot <- df %>% filter(!is.na(shm_freq))
      validate(need(nrow(df_plot) > 0, "No SHM data available."))

      angle     <- as.numeric(input$x_axis_angle %||% "45")
      hjust_val <- if (angle == 0) 0.5 else 1
      base_font <- input$base_font_size %||% 12
      leg_pos   <- input$legend %||% "right"
      show_jit  <- isTRUE(input$show_jitter)
      jit_size  <- input$jitter_size %||% 0.5

      p <- ggplot(df_plot, aes(x = .data[[input$group_by]], y = shm_freq, fill = .data[[input$group_by]]))

      if (input$shm_plot_type == "violin") {
        p <- p + geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill = "white")
        if (show_jit) p <- p + geom_jitter(size = jit_size, alpha = 0.4, width = 0.1)
      } else if (input$shm_plot_type == "boxplot") {
        p <- p + geom_boxplot()
        if (show_jit) p <- p + geom_jitter(size = jit_size, alpha = 0.4, width = 0.1)
      } else {
        p <- ggplot(df_plot, aes(x = shm_freq, fill = .data[[input$group_by]])) +
             geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
             facet_wrap(as.formula(paste0("~", input$group_by)))
      }

      p + theme_classic(base_size = base_font) +
          theme(
            axis.text.x     = element_text(angle = angle, hjust = hjust_val),
            legend.position = leg_pos
          ) +
          labs(y = "SHM Frequency (100 - %V-Identity)", x = input$group_by,
               title = "Somatic Hypermutation Rate", fill = input$group_by)
    }

    pick_isotype_col <- function(df) {
      candidates <- c("BCR_pair_c_gene", "BCR_IGH_c_gene")
      found <- candidates[candidates %in% names(df)][1]
      if (is.na(found)) {
        suffixed <- grep("c_gene(\\.x)?$", names(df), value = TRUE)
        if (length(suffixed) > 0) found <- suffixed[1]
      }
      if (is.na(found)) return(NULL) else return(found)
    }

    build_isotype_plot <- function() {
      df <- bcr_data()
      isotype_col <- pick_isotype_col(df)
      validate(need(!is.null(input$group_by) && nzchar(input$group_by), "Select a grouping variable."))
      validate(need(!is.null(isotype_col),
                    paste0("Isotype column not found. Available BCR columns: ",
                           paste(grep("^BCR_", names(df), value = TRUE), collapse = ", "))))

      df_isotype <- df %>%
        filter(!is.na(.data[[isotype_col]]) & .data[[isotype_col]] != "")
      validate(need(nrow(df_isotype) > 0, "No isotype data available."))

      angle     <- as.numeric(input$x_axis_angle %||% "45")
      hjust_val <- if (angle == 0) 0.5 else 1
      base_font <- input$base_font_size %||% 12
      leg_pos   <- input$legend %||% "right"

      p <- ggplot(df_isotype, aes(x = .data[[input$group_by]], fill = .data[[isotype_col]]))

      if (input$show_percentage) {
        p <- p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent)
        y_lab <- "Proportion"
      } else {
        p <- p + geom_bar(position = "stack")
        y_lab <- "Count"
      }

      p + theme_classic(base_size = base_font) +
          theme(
            axis.text.x     = element_text(angle = angle, hjust = hjust_val),
            legend.position = leg_pos
          ) +
          labs(x = input$group_by, y = y_lab,
               fill = "Isotype (C-gene)", title = "Isotype Distribution")
    }

    output$shm_plot <- renderPlot({
      build_shm_plot()
    }, width = reactive(input$plot_width %||% 700), height = reactive(input$plot_height %||% 500))

    output$shm_stats <- renderTable({
      df <- bcr_data()
      validate(need(!is.null(input$group_by) && nzchar(input$group_by), "Select a grouping variable."))
      validate(need("shm_freq" %in% names(df), "SHM data not found in the repertoire table."))

      df %>%
        filter(!is.na(shm_freq)) %>%
        group_by(.data[[input$group_by]]) %>%
        summarise(
          Mean   = mean(shm_freq),
          Median = median(shm_freq),
          SD     = sd(shm_freq),
          Count  = n(),
          .groups = "drop"
        )
    })

    output$isotype_plot <- renderPlot({
      build_isotype_plot()
    }, width = reactive(input$plot_width %||% 700), height = reactive(input$plot_height %||% 500))

    output$isotype_table <- DT::renderDT({
        df <- bcr_data()
        isotype_col <- pick_isotype_col(df)
        validate(need(!is.null(input$group_by) && nzchar(input$group_by), "Select a grouping variable."))
        validate(need(!is.null(isotype_col), "Isotype column not found in the repertoire table."))

        counts <- df %>%
            group_by(.data[[input$group_by]], .data[[isotype_col]]) %>%
            summarise(Count = n(), .groups = "drop")

        DT::datatable(counts, options = list(pageLength = 10))
    })

    output$download_shm_plot <- downloadHandler(
        filename = function() { paste0("shm_plot_", input$group_by, "_", Sys.Date(), ".pptx") },
        content = function(file) {
            p <- build_shm_plot()
            save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
        }
    )

    output$download_isotype_plot <- downloadHandler(
        filename = function() { paste0("isotype_plot_", input$group_by, "_", Sys.Date(), ".pptx") },
        content = function(file) {
            p <- build_isotype_plot()
            save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
        }
    )
  })
}
