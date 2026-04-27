# CDR3 Physicochemical Properties Module

# Kyte-Doolittle Hydrophobicity Scale
kd_scale <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, 
              H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, 
              T=-0.7, W=-0.9, Y=-1.3, V=4.2)

# Charge at pH 7
charge_scale <- c(R=1, K=1, H=0.1, D=-1, E=-1)

# Helper function to compute properties
compute_aa_properties <- function(sequences) {
  # Initialize vectors
  hydro <- numeric(length(sequences))
  charge <- numeric(length(sequences))
  length_aa <- numeric(length(sequences))
  
  for(i in seq_along(sequences)) {
    seq <- sequences[i]
    if(is.na(seq) || nchar(seq) == 0 || grepl("NA", seq) || seq == "") {
      hydro[i] <- NA; charge[i] <- NA; length_aa[i] <- NA
      next
    }
    
    # Split to chars
    chars <- strsplit(seq, "")[[1]]
    valid_chars <- chars[chars %in% names(kd_scale)]
    
    if(length(valid_chars) == 0) {
      hydro[i] <- NA; charge[i] <- NA; length_aa[i] <- NA
      next
    }
    
    length_aa[i] <- length(valid_chars)
    hydro[i] <- mean(kd_scale[valid_chars], na.rm = TRUE)
    
    # Sum of charges / length (normalized charge Density)
    seq_charges <- charge_scale[valid_chars]
    seq_charges[is.na(seq_charges)] <- 0
    charge[i] <- sum(seq_charges)
  }
  
  return(data.frame(Length = length_aa, Hydrophobicity = hydro, NetCharge = charge))
}

# --- UI Definition ---
cdr3PropertiesUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("CDR3 Physicochemical Properties"),
      p("Analyze the biochemical properties of CDR3 amino acid sequences (e.g., Hydrophobicity, Net Charge, Length)."),
      vdjType(ns),
      selectInput(ns("chain_type"), "Target Chain:", choices = c("Alpha/Heavy", "Beta/Light"), selected = "Beta/Light"),
      groupByInput(ns),
      selectInput(ns("property"), "Property to Visualize:", choices = c("Hydrophobicity (Kyte-Doolittle)", "Net Charge", "Sequence Length")),
      hr(),
      actionButton(ns("run_calc"), "Compute Properties", class = "btn-primary", width = "100%", icon = icon("flask")),
      hr(),
      h5("Plot Options"),
      commonPlotOptions(ns),
      jitterPointsInput(ns),
      hr(),
      # downloadButton moved to mainPanel
    ),
    mainPanel(
      h3("Biochemical Profile Distribution"),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns("prop_plot"), height = "500px"),
      hr(),
      h4("Summary Statistics"),
      DT::DTOutput(ns("prop_table"))
    )
  )
}

# --- Server Logic ---
cdr3PropertiesServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(myReactives$seurat_object, {
       update_group_by_select_input(session, myReactives)
    })
    
    prop_data <- eventReactive(input$run_calc, {
      req(input$vdj_type)
      df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(df)
      
      # Determine sequence column
      seq_col <- ""
      if (input$vdj_type == "tcr") {
        seq_col <- if (input$chain_type == "Alpha/Heavy") "TCR_TRA_cdr3" else "TCR_TRB_cdr3"
      } else {
         # BCR
        seq_col <- if (input$chain_type == "Alpha/Heavy") "BCR_IGH_cdr3_aa" else "BCR_IGK_cdr3_aa"
      }
      
      if (!seq_col %in% names(df)) {
        seq_col_fallback <- paste0(toupper(input$vdj_type), "_pair_CTaa")
        if (seq_col_fallback %in% names(df)) {
             # naive split for pairs, take the requested part
             split_idx <- if (input$chain_type == "Alpha/Heavy") 1 else 2
             df[[seq_col]] <- sapply(strsplit(as.character(df[[seq_col_fallback]]), "_"), function(x) if(length(x) >= split_idx) x[split_idx] else NA)
        } else {
             showNotification(paste("Sequence column", seq_col, "not found."), type="error")
             return(NULL)
        }
      }
      
      group_col <- input$group_by
      if (!"barcode" %in% names(df)) {
        showNotification("Barcode column not found in repertoire data.", type = "error")
        return(NULL)
      }
      if (!group_col %in% names(df)) {
        so <- myReactives$seurat_object
        if (!is.null(so) && group_col %in% names(so@meta.data)) {
          meta_sub <- data.frame(barcode = rownames(so@meta.data),
                                 stringsAsFactors = FALSE)
          meta_sub[[group_col]] <- so@meta.data[[group_col]]
          df <- dplyr::left_join(df, meta_sub, by = "barcode")
        }
      }
      if (!group_col %in% names(df)) {
        showNotification(paste("Group column", group_col, "not found in repertoire or metadata."), type = "error")
        return(NULL)
      }
      
      withProgress(message = "Computing structural properties...", {
        # --- Step 1: filter to valid rows using base R to avoid S4/Rle issues ---
        keep <- !is.na(df[[seq_col]]) & !is.na(df[[group_col]]) &
                df[[seq_col]] != "" & df[[seq_col]] != "NA"
        df_unique <- df[keep, , drop = FALSE]

        # --- Step 2: strip ALL Rle / S4 / list columns to plain base-R vectors ---
        df_unique <- as.data.frame(lapply(df_unique, function(col) {
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

        # --- Step 3: deduplicate using base R (avoids dplyr re-creating Rle internally) ---
        if ("raw_clonotype_id" %in% names(df_unique)) {
          dup_key <- paste(df_unique[[group_col]], df_unique[["raw_clonotype_id"]], sep = "\t")
          df_unique <- df_unique[!duplicated(dup_key), , drop = FALSE]
        } else {
          df_unique <- df_unique[!duplicated(df_unique[["barcode"]]), , drop = FALSE]
        }

        # --- Step 4: keep only the columns needed ---
        keep_cols <- intersect(c("barcode", group_col, seq_col), names(df_unique))
        df_unique <- df_unique[, keep_cols, drop = FALSE]
        
        props <- compute_aa_properties(df_unique[[seq_col]])
        df_unique <- bind_cols(df_unique, props)
        
        return(df_unique %>% filter(!is.na(Length)))
      })
    })
    
    prop_plot_obj <- reactive({
      df <- prop_data()
      validate(need(!is.null(df), "No data available. Click 'Compute Properties' after selecting valid inputs."))
      
      group_col <- input$group_by
      y_col <- if(grepl("Hydrophobicity", input$property)) {
         "Hydrophobicity" 
      } else if(grepl("Charge", input$property)) {
         "NetCharge"
      } else "Length"
      
      y_label <- if (y_col == "Hydrophobicity") "Mean Kyte-Doolittle Score" else if (y_col == "NetCharge") "Net Charge (pH 7.0)" else "Amino Acid Length"
      
      p_cdr3 <- ggplot(df, aes(x = .data[[group_col]], y = .data[[y_col]], fill = .data[[group_col]])) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
        theme_classic(base_size = 14) +
        labs(title = paste("CDR3", input$property, "Distribution"), x = tools::toTitleCase(gsub("_", " ", group_col)), y = y_label) +
        theme(
          axis.text.x = element_text(angle = as.numeric(input$x_axis_angle %||% "45"), hjust = 1, face = "bold"),
          legend.position = input$legend %||% "none"
        )
      if (isTRUE(input$show_jitter)) {
        p_cdr3 <- p_cdr3 + geom_jitter(width = 0.2, size = input$jitter_size %||% 0.5, alpha = 0.5, color = "black")
      }
      p_cdr3
    })
    
    output$prop_plot <- renderPlot({
      prop_plot_obj()
    }, res=96, width=reactive(input$plot_width %||% 500L), height=reactive(input$plot_height %||% 500L))
    
    output$prop_table <- DT::renderDT({
       df <- prop_data()
       validate(need(!is.null(df), "No data available. Click 'Compute Properties' after selecting valid inputs."))
       group_col <- input$group_by
       
       y_col <- if(grepl("Hydrophobicity", input$property)) {
         "Hydrophobicity" 
       } else if(grepl("Charge", input$property)) {
         "NetCharge"
       } else "Length"
       
       summary_df <- df %>%
         group_by(.data[[group_col]]) %>%
         summarise(
           Mean = mean(.data[[y_col]], na.rm=TRUE),
           Median = median(.data[[y_col]], na.rm=TRUE),
           SD = sd(.data[[y_col]], na.rm=TRUE),
           Unique_Clones = n()
         ) %>%
         mutate(across(where(is.numeric), ~round(., 3)))
         
       DT::datatable(summary_df, options = list(pageLength = 10, dom = "t"))
    })
    
    output$download_plot <- downloadHandler(
      filename = function() {
        y_col_name <- if (grepl("Hydrophobicity", input$property)) "Hydrophobicity"
                      else if (grepl("Charge", input$property)) "NetCharge"
                      else "Length"
        paste0(input$vdj_type, "_cdr3_", y_col_name, "_", Sys.Date(), ".pptx")
      },
      content = function(file) {
         save_plot_as_pptx(file, prop_plot_obj(), input$plot_width, input$plot_height)
      }
    )
  })
}
