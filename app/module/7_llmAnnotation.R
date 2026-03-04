# LLM-Based Automated Cell Type Annotation Module
# Uses OpenAI API to predict cell types from cluster marker genes

library(shiny)
library(Seurat)
library(dplyr)
library(httr)
library(jsonlite)
library(DT)

# --- UI Definition ---
llmAnnotationUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h4("LLM Cell Type Annotation"),
      p("Automatically annotate clusters using Large Language Models (OpenAI) based on top marker genes."),
      passwordInput(ns("api_key"), "OpenAI API Key:", value = ""),
      selectInput(ns("model_name"), "Model:", choices = c("gpt-4o-mini", "gpt-4o", "gpt-3.5-turbo"), selected = "gpt-4o-mini"),
      hr(),
      h5("Clustering Input"),
      selectInput(ns("cluster_col"), "Target Cluster Column:", choices = NULL),
      numericInput(ns("top_n_genes"), "Top Genes per Cluster for LLM:", value = 15, min = 5, max = 50),
      actionButton(ns("run_llm"), "Run Auto-Annotation", icon = icon("robot"), class = "btn-primary", width = "100%"),
      hr(),
      h5("Save Results"),
      p("Apply the generated labels back to the Seurat object metadata."),
      actionButton(ns("apply_annotation"), "Apply to Seurat Metadata", icon = icon("check"), class = "btn-success", width = "100%")
    ),
    mainPanel(
      card(
        card_header("Annotation Results"),
        DT::DTOutput(ns("annotation_table"))
      ),
      br(),
      card(
        card_header("System Log"),
        verbatimTextOutput(ns("log_output"))
      )
    )
  )
}

# --- Server Logic ---
llmAnnotationServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Observe Seurat object to update cluster choices
    observe({
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      meta <- so@meta.data
      
      categorical_cols <- names(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      default_col <- if ("seurat_clusters" %in% categorical_cols) "seurat_clusters" else categorical_cols[1]
      
      updateSelectInput(session, "cluster_col", choices = categorical_cols, selected = default_col)
    })
    
    # Reactive values
    annotation_results <- reactiveVal(NULL)
    log_messages <- reactiveVal("Awaiting user execution...")
    
    observeEvent(input$run_llm, {
      req(myReactives$seurat_object, input$cluster_col)
      
      if (trimws(input$api_key) == "") {
        showNotification("Please provide a valid OpenAI API Key.", type = "error")
        log_messages("Error: OpenAI API key is missing.")
        return()
      }
      
      withProgress(message = "Running LLM Auto-Annotation...", value = 0, {
        so <- myReactives$seurat_object
        cluster_col <- input$cluster_col
        top_n <- input$top_n_genes
        
        incProgress(0.1, detail = "Setting identities...")
        Idents(so) <- cluster_col
        
        log_messages("Finding marker genes for all clusters...")
        incProgress(0.2, detail = paste("Finding marker genes for:", cluster_col))
        
        markers <- tryCatch({
          # only positive markers to reduce noise and computation time
          FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
        }, error = function(e) {
          log_messages(paste("Error finding markers:", e$message))
          NULL
        })
        
        if (is.null(markers) || nrow(markers) == 0) {
          showNotification("FindAllMarkers returned no genes or failed.", type = "error")
          return()
        }
        
        incProgress(0.5, detail = "Formatting prompt for LLM...")
        log_messages("Processing markers and formatting prompt...")
        
        top_markers <- markers %>%
          group_by(cluster) %>%
          slice_max(n = top_n, order_by = avg_log2FC)
        
        prompt_data <- ""
        for (clust in unique(top_markers$cluster)) {
          genes <- top_markers %>% filter(cluster == clust) %>% pull(gene) %>% paste(collapse = ", ")
          prompt_data <- paste0(prompt_data, "- Cluster '", clust, "': ", genes, "\n")
        }
        
        system_prompt <- "You are an expert single-cell bioinformatician. I will provide you with a list of cell clusters and their top positive marker genes (based on log fold change). Identify the most likely biological cell type for each cluster. Respond strictly in JSON format. The JSON must contain a single key 'annotations' pointing to an array of objects. Each object must have three keys: 'Cluster' (string), 'CellType' (string), and 'Reasoning' (string explaining your choice briefly)."
        user_prompt <- paste0("Here are the marker genes per cluster:\n", prompt_data, "\nPlease provide the JSON output.")
        
        incProgress(0.7, detail = paste("Contacting OpenAI API (", input$model_name, ")..."))
        log_messages(paste("Sending request to OpenAI using model:", input$model_name))
        
        json_body <- list(
          model = input$model_name,
          messages = list(
            list(role = "system", content = system_prompt),
            list(role = "user", content = user_prompt)
          ),
          temperature = 0.2, 
          response_format = list(type = "json_object")
        )
        
        resp <- tryCatch({
          httr::POST(
            url = "https://api.openai.com/v1/chat/completions",
            httr::add_headers(Authorization = paste("Bearer", trimws(input$api_key))),
            httr::content_type_json(),
            body = jsonlite::toJSON(json_body, auto_unbox = TRUE),
            httr::timeout(60) # 60 seconds timeout
          )
        }, error = function(e) {
          log_messages(paste("Network error contacting OpenAI:", e$message))
          NULL
        })
        
        if (is.null(resp)) {
          showNotification("Failed to connect to OpenAI.", type = "error")
          return()
        }
        
        incProgress(0.9, detail = "Parsing LLM response...")
        
        if (httr::status_code(resp) == 200) {
          resp_content <- httr::content(resp, as = "text", encoding = "UTF-8")
          parsed <- jsonlite::fromJSON(resp_content)
          llm_text <- parsed$choices$message$content
          
          res_json <- tryCatch({
            jsonlite::fromJSON(llm_text)
          }, error = function(e) {
            log_messages(paste("JSON Parsing Error:", e$message, "\nRaw Text:", llm_text))
            NULL
          })
          
          if (!is.null(res_json) && "annotations" %in% names(res_json)) {
            res_df <- as.data.frame(res_json$annotations)
            annotation_results(res_df)
            log_messages("Successfully received and parsed annotations from LLM.")
            showNotification("LLM Annotation completed.", type = "message")
          } else {
            log_messages("Valid JSON returned, but 'annotations' key was missing.")
            showNotification("Invalid JSON format from LLM.", type = "error")
          }
        } else {
          err_msg <- httr::content(resp, as = "text", encoding = "UTF-8")
          log_messages(paste("API Error:", httr::status_code(resp), "\nDetails:", err_msg))
          showNotification(paste("API Error:", httr::status_code(resp)), type = "error")
        }
      })
    })
    
    output$annotation_table <- DT::renderDT({
      req(annotation_results())
      DT::datatable(annotation_results(), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$log_output <- renderText({
      log_messages()
    })
    
    observeEvent(input$apply_annotation, {
      req(myReactives$seurat_object, annotation_results(), input$cluster_col)
      
      ann_df <- annotation_results()
      if (!"Cluster" %in% names(ann_df) || !"CellType" %in% names(ann_df)) {
        showNotification("Cannot apply: Results map is incomplete.", type = "error")
        return()
      }
      
      so <- myReactives$seurat_object
      cluster_col <- input$cluster_col
      
      # Build dictionary
      mapping <- setNames(ann_df$CellType, ann_df$Cluster)
      
      # Vectorized replacement
      current_clusters <- as.character(so@meta.data[[cluster_col]])
      new_labels <- mapping[current_clusters]
      
      # Handle potential NAs if LLM missed a cluster
      missing_idx <- is.na(new_labels)
      new_labels[missing_idx] <- current_clusters[missing_idx]
      
      col_name <- paste0("llm_annotation_", Sys.Date())
      so@meta.data[[col_name]] <- as.factor(new_labels)
      
      # Overwrite Seurat Object
      myReactives$seurat_object <- so
      
      # Trigger Universal UI update
      myReactives$grouping_updated <- isolate(myReactives$grouping_updated) + 1
      
      log_messages(paste("Applied annotations as new column in metadata:", col_name))
      showNotification(paste("Annotations applied to metadata column:", col_name), type = "message", duration = 8)
    })
  })
}
