# 必要なライブラリがapp.Rで読み込まれていることを前提としています
# library(shiny)
# library(Seurat)
# library(ggplot2)
# library(dplyr)
# library(DT)
# library(shinyjs)
# library(shinyalert)

#=================================================
# UI部分
#=================================================
reductionPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      useShinyjs(),
      
      h4("Plot Display Options"),
      reductionInput(ns),
      div(style = "display: flex; align-items: center; gap: 10px;",
        div(style = "flex-grow: 1;", groupByInput(ns)),
        div(style = "flex-shrink: 0;", actionButton(ns("rename_btn"), "Rename", icon = icon("pen"), style = "height: 100%;"))
      ),
      # ★ 論文用の図作成に適したnumericInputに変更
      pointSizeInput(ns),
      labelSizeInput(ns),
      commonPlotOptions(ns),
      
      hr(),

      # --- 検索＆ハイライト機能 ---
      h4("Search & Highlight"),
      checkboxInput(ns("clonotype_search_toggle"), "Enable Highlight by Search", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s'] == true", ns("clonotype_search_toggle")),
        selectInput(ns("search_source"), "1. Select Data Source", choices = NULL),
        selectInput(ns("search_column"), "2. Select Column to Search", choices = NULL),
        selectizeInput(ns("search_values"), "3. Select Value(s) to Highlight", choices = NULL, multiple = TRUE, options = list(placeholder = 'Type or select values')),
        actionButton(ns("search_highlight_btn"), "Apply Highlight", icon = icon("search"))
      ),
      hr(),
      h4("Re-run Annotation"),
      selectInput(ns("annotation_method"), "Annotation Method", choices = c("scType (Marker-based)" = "sctype", "SingleR (Reference-based)" = "singler")),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'sctype'", ns("annotation_method")),
        selectInput(ns("sctype_tissue"), "Target Tissue (scType)", choices = c("Immune system", "Pancreas", "Heart", "Brain", "Placenta", "Kidney", "Liver", "Lung", "Muscle", "Intestine"))
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'singler'", ns("annotation_method")),
        selectInput(ns("singler_ref"), "Reference Dataset (celldex)",
                    choices = list("Human Primary Cell Atlas (Broad)" = "HumanPrimaryCellAtlasData", "Monaco Immune Data (Immune)" = "MonacoImmuneData", "Database of Immune Cell Expression (DICE)" = "DatabaseImmuneCellExpressionData", "Novershtern Hematopoietic Data" = "NovershternHematopoieticData", "Blueprint Encode Data" = "BlueprintEncodeData", "Mouse RNA-seq Data (Broad)" = "MouseRNAseqData", "ImmGen (Mouse Immune)" = "ImmGenData"))
      ),
      actionButton(ns("run_annotation_btn"), "Run Annotation", icon = icon("tags"), class = "btn-info")
    ),
    mainPanel(
      h3("Plot"),
      div(style = "display: flex; justify-content: space-between; align-items: center;",
          downloadButton(ns("download_plot"), "Download plot (.pdf)"),
          div(
            actionButton(ns("reset_subset_btn"), "Reset Subset", icon = icon("undo"), class = "btn-secondary"),
            actionButton(ns("subset_btn"), "Subset Brushed Cells", icon = icon("cut"), class = "btn-warning")
          )
      ),
      br(),
      plotOutput(ns("plot"), brush = brushOpts(id = ns("plot_brush"), resetOnNew = TRUE)),
      hr(),
      h3("Coordinates Table"),
      downloadButton(ns("download_table"), "Download Table (.csv)"),
      DT::dataTableOutput(ns("coordinates_table"))
    )
  )
}

#=================================================
# Server部分
#=================================================
reductionPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    highlighted_cells <- reactiveVal(NULL)

    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      
      # TCR/BCR status列の追加 (変更なし)
      if (!"TCR_BCR_status" %in% names(so@meta.data)) {
        clonotype_cols <- c("CTaa", "CTnt", "raw_clonotype_id", "clonotype_id", "clonotype") 
        if (any(clonotype_cols %in% names(so@meta.data))) {
          col_to_use <- clonotype_cols[clonotype_cols %in% names(so@meta.data)][1]
          so$TCR_BCR_status <- ifelse(is.na(so@meta.data[[col_to_use]]) | so@meta.data[[col_to_use]] == "None", "Absent", "Present")
          so$TCR_BCR_status <- factor(so$TCR_BCR_status, levels = c("Present", "Absent"))
          myReactives$seurat_object <- so
        }
      }
      
      update_reduction_choices(session, myReactives)
      # ✅ 修正済みの関数を呼び出し
      update_group_by_select_input(session, myReactives)
      
      search_sources <- c("Seurat Metadata" = "metadata")
      if (!is.null(myReactives$tcr_df)) search_sources <- c(search_sources, "TCR Clonotypes" = "tcr_df")
      if (!is.null(myReactives$bcr_df)) search_sources <- c(search_sources, "BCR Clonotypes" = "bcr_df")
      updateSelectInput(session, "search_source", choices = search_sources)
      
    }, ignoreNULL = TRUE)

    # --- 検索機能のロジック ---
    
    # ハイライト機能を無効にしたら、ハイライト細胞リストをクリア
    observeEvent(input$clonotype_search_toggle, {
      if (!isTRUE(input$clonotype_search_toggle)) {
        highlighted_cells(NULL)
      }
    })
    
    # 検索ソースに応じて、検索対象の列を更新
    observeEvent(input$search_source, {
      req(input$search_source)
      source_name <- input$search_source
      
      choices <- if (source_name == "metadata") {
        req(myReactives$seurat_object)
        meta <- myReactives$seurat_object@meta.data
        is_categorical <- sapply(names(meta), function(col) is.character(meta[[col]]) || is.factor(meta[[col]]) || is.logical(meta[[col]]))
        names(meta)[is_categorical]
      } else if (!is.null(myReactives[[source_name]])) {
        names(myReactives[[source_name]])
      } else {
        NULL
      }
      updateSelectInput(session, "search_column", choices = choices)
    }, ignoreInit = TRUE)
    
    # 検索列に応じて、検索候補の値を更新
    observeEvent(input$search_column, {
      req(input$search_source, input$search_column, nzchar(input$search_column))
      source_name <- input$search_source
      col_name <- input$search_column
      
      unique_vals <- if (source_name == "metadata") {
        req(myReactives$seurat_object)
        col_data <- myReactives$seurat_object@meta.data[[col_name]]
        if (is.logical(col_data)) c("TRUE", "FALSE") else unique(col_data)
      } else if (!is.null(myReactives[[source_name]])) {
        unique(myReactives[[source_name]][[col_name]])
      } else {
        NULL
      }
      updateSelectizeInput(session, "search_values", choices = sort(na.omit(unique_vals)), server = TRUE)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    # ハイライトボタンが押された時の処理
    observeEvent(input$search_highlight_btn, {
      req(input$search_source, input$search_column, input$search_values)
      source_name <- input$search_source
      col_name <- input$search_column
      vals_to_find <- input$search_values
      
      barcodes_to_highlight <- if (source_name == "metadata") {
        meta_df <- myReactives$seurat_object@meta.data
        if (is.logical(meta_df[[col_name]])) {
          vals_to_find <- as.logical(vals_to_find)
        }
        rownames(meta_df)[meta_df[[col_name]] %in% vals_to_find]
      } else {
        search_df <- myReactives[[source_name]]
        search_df %>%
          filter(.data[[col_name]] %in% vals_to_find) %>%
          pull(barcode)
      }
      
      valid_cells <- intersect(barcodes_to_highlight, colnames(myReactives$seurat_object))
      highlighted_cells(valid_cells)
    })

    # --- プロット描画ロジック (簡潔版) ---
    plot <- reactive({
      req(myReactives$seurat_object, input$reduction)
      so <- myReactives$seurat_object
      
      # 検索によるハイライトが有効な場合
      if (isTRUE(input$clonotype_search_toggle) && !is.null(highlighted_cells()) && length(highlighted_cells()) > 0) {
        p <- DimPlot(
          so,
          reduction = input$reduction,
          cells.highlight = highlighted_cells(),
          cols.highlight = 'darkblue',
          pt.size = input$point_size,
          label = FALSE
        ) + guides(color = "none") + theme(legend.position = input$legend)

      } else { # 通常のプロット
        req(input$group_by)
        p <- DimPlot(
          so,
          reduction = input$reduction,
          group.by = input$group_by,
          label = TRUE,
          pt.size = input$point_size,
          label.size = input$label_size
        ) + theme(legend.position = input$legend)
      }
      return(p)
    })

    output$plot <- renderPlot({ plot() }, width = reactive(input$plot_width), height = reactive(input$plot_height))
    
    # --- Rename機能 (変更なし) ---
    observeEvent(input$rename_btn, {
      req(myReactives$seurat_object, input$group_by)
      so <- myReactives$seurat_object
      group_by_col <- input$group_by
      current_levels <- if (is.factor(so[[group_by_col, drop = TRUE]])) { levels(so[[group_by_col, drop = TRUE]]) } else { sort(unique(so[[group_by_col, drop = TRUE]])) }
      showModal(modalDialog(
        title = "Rename/Regroup Clusters",
        textInput(ns("new_group_name"), "New Group Name", placeholder = "e.g., Celltype_manual"), hr(),
        lapply(current_levels, function(level) { textInput(ns(paste0("rename_", level)), label = paste("New name for cluster:", level), value = level) }),
        footer = tagList(modalButton("Cancel"), actionButton(ns("save_rename"), "Save"))
      ))
    })

    observeEvent(input$save_rename, {
      req(input$new_group_name)
      if (grepl("\\s|[^A-Za-z0-9_.]", input$new_group_name)) {
        if(require(shinyalert)) shinyalert("Error", "Group name cannot contain spaces or special characters.", type = "error")
        return()
      }
      so <- myReactives$seurat_object
      group_by_col <- input$group_by
      current_levels <- if (is.factor(so[[group_by_col, drop = TRUE]])) { levels(so[[group_by_col, drop = TRUE]]) } else { sort(unique(so[[group_by_col, drop = TRUE]])) }
      
      new_names_map <- setNames(sapply(current_levels, function(level) input[[paste0("rename_", level)]]), current_levels)
      original_vector <- so[[group_by_col, drop = TRUE]]
      new_group_vector <- factor(unname(new_names_map[as.character(original_vector)]))
      
      so@meta.data[[input$new_group_name]] <- new_group_vector
      myReactives$seurat_object <- so
      
      # ✅ 修正済みの関数を呼び出し
      update_group_by_select_input(session, myReactives, selected = input$new_group_name)
      removeModal()
    })

    # --- テーブルとダウンロード機能 (変更なし) ---
    coordinates_data <- reactive({
      req(myReactives$seurat_object, input$reduction, input$group_by)
      so <- myReactives$seurat_object
      reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
      table_data <- cbind(so@meta.data, reduction_data)
      table_data$Barcode <- rownames(table_data)
      reduction_prefix <- sub("_$", "", colnames(reduction_data)[1]) %>% sub("(?<=[A-Za-z])_1$", "", ., perl=TRUE)
      # For UMAP_1, regex sub("_1$", "") returns "UMAP". 
      # Then starts_with("UMAP_") matches UMAP_1, UMAP_2
      if (grepl("_1$", colnames(reduction_data)[1])) {
          reduction_prefix <- sub("_1$", "", colnames(reduction_data)[1])
      } else {
          reduction_prefix <- "UMAP" # fallback
      }
      table_data %>% select(Barcode, all_of(input$group_by), starts_with(paste0(reduction_prefix, "_")))
    })
    
    # --- Interactive Subsetting Logic ---
    observeEvent(input$subset_btn, {
      req(myReactives$seurat_object, input$plot_brush)
      so <- myReactives$seurat_object
      
      coords <- coordinates_data()
      # Find the exact column names for x and y
      # They typically look like UMAP_1 and UMAP_2
      coord_cols <- grep("_[12]$", colnames(coords), value = TRUE)
      
      if(length(coord_cols) < 2) {
         if(require(shinyalert)) shinyalert("Error", "Could not determine X/Y axes for selection.", type="error")
         return()
      }
      x_col <- coord_cols[1]
      y_col <- coord_cols[2]
      
      selected_rows <- brushedPoints(coords, input$plot_brush, xvar = x_col, yvar = y_col)
      
      if(nrow(selected_rows) == 0){
        if(require(shinyalert)) shinyalert("Warning", "No cells selected. Please draw a box on the plot first.", type="warning")
        return()
      }
      
      selected_barcodes <- selected_rows$Barcode
      
      showModal(modalDialog(
        title = "Confirm Subset",
        HTML(paste0("You have selected <b>", length(selected_barcodes), "</b> cells. Do you want to strictly keep only these cells?<br><br><b>Note:</b> This will permanently alter the current session. All other tabs will immediately update to reflect only these cells.")),
        br(), br(),
        checkboxInput(ns("rerun_umap_checkbox"), "Re-run Normalization, PCA, and UMAP after subsetting", value = FALSE),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_subset_btn"), "Yes, Subset Data", class = "btn-danger")
        )
      ))
      
      observeEvent(input$confirm_subset_btn, {
         removeModal()
         
         # Backup original seurat object if not done yet
         if (is.null(myReactives$full_seurat_object)) {
             myReactives$full_seurat_object <- so
         }

         # Perform Subsetting
         so_sub <- subset(so, cells = selected_barcodes)
         
         if (input$rerun_umap_checkbox) {
            withProgress(message = "Re-running standard pipeline...", value = 0, {
              incProgress(0.2, detail = "NormalizeData")
              so_sub <- NormalizeData(so_sub)
              incProgress(0.4, detail = "FindVariableFeatures")
              so_sub <- FindVariableFeatures(so_sub)
              incProgress(0.6, detail = "ScaleData & PCA")
              so_sub <- ScaleData(so_sub)
              so_sub <- RunPCA(so_sub)
              incProgress(0.8, detail = "Neighbors & Clusters")
              so_sub <- FindNeighbors(so_sub, dims = 1:10)
              so_sub <- FindClusters(so_sub)
              incProgress(0.9, detail = "UMAP")
              so_sub <- RunUMAP(so_sub, dims = 1:10)
            })
         }
         
         # Update global Reactives
         myReactives$seurat_object <- so_sub
         if(require(shinyalert)) shinyalert("Success", paste("Dataset successfully subsetted to", ncol(so_sub), "cells."), type="success")
      }, ignoreInit = TRUE, once = TRUE)
    })

    # --- Reset Subsetting Logic ---
    observeEvent(input$reset_subset_btn, {
        if (!is.null(myReactives$full_seurat_object)) {
            myReactives$seurat_object <- myReactives$full_seurat_object
            myReactives$full_seurat_object <- NULL
            if(require(shinyalert)) shinyalert("Dataset Reset", "The dataset has been returned to its original state before subsetting.", type="success")
        } else {
            if(require(shinyalert)) shinyalert("Info", "No subset has been performed yet.", type="info")
        }
    })

    # --- Re-run Annotation Logic ---
    observeEvent(input$run_annotation_btn, {
       req(myReactives$seurat_object)
       so <- myReactives$seurat_object
       # Find if the helper functions exist (they are defined in 1_uploadCellranger.R)
       if (!exists("run_sctype_and_update_seurat")) {
           if(require(shinyalert)) shinyalert("Error", "Required annotation functions not found.", type="error")
           return()
       }
       withProgress(message = "Running Annotation...", value = 0, {
           if (input$annotation_method == "sctype") {
               incProgress(0.2, detail = "Running scType cell typing...")
               so <- run_sctype_and_update_seurat(so, tissue = input$sctype_tissue)
           } else if (input$annotation_method == "singler") {
               incProgress(0.2, detail = "Running SingleR cell typing...")
               so <- run_singler_and_update_seurat(so, ref_name = input$singler_ref)
           }
           myReactives$seurat_object <- so
       })
       # Refresh group_by choices
       update_group_by_select_input(session, myReactives, selected = if(input$annotation_method == "sctype") "sctype_celltype" else "singler_celltype")
       if(require(shinyalert)) shinyalert("Success", "Annotation completed successfully.", type="success")
    })
    
    output$coordinates_table <- DT::renderDataTable({ coordinates_data() })
    
    output$download_plot <- downloadHandler(
      filename = function() { paste0("DimPlot_", input$reduction, "_", input$group_by, ".pdf") },
      content = function(file) {
        req(plot())
        ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300, device = "pdf")
      }
    )
    
    output$download_table <- downloadHandler(
      filename = function() { paste0("coordinates_", input$reduction, "_", Sys.Date(), ".csv") },
      content = function(file) { write.csv(coordinates_data(), file, row.names = TRUE) }
    )
  })
}


#=================================================
# ヘルパー関数群
#=================================================

update_reduction_choices <- function(session, myReactives) {
  req(myReactives$seurat_object)
  reduction_names <- names(myReactives$seurat_object@reductions)
  choices <- stats::setNames(reduction_names, toupper(reduction_names))
  default_selection <- if ("umap" %in% reduction_names) "umap" else reduction_names[1]
  updateSelectInput(session, "reduction", choices = choices, selected = default_selection)
}

# ✅ 以前の議論で完成した最終版の関数
update_group_by_select_input <- function(session, myReactives, selected = NULL) {
  req(myReactives$seurat_object)
  meta_data <- myReactives$seurat_object@meta.data
  
  is_categorical <- sapply(names(meta_data), function(col) {
    is.character(meta_data[[col]]) || is.factor(meta_data[[col]]) || is.logical(meta_data[[col]])
  })
  potential_choices <- names(meta_data)[is_categorical]

  minus_column <- c("orig.ident", "barcode")
  choices <- setdiff(potential_choices, minus_column)
  choices <- choices[!grepl("^(RNA_snn_res\\.|TCR|BCR|raw_clonotype|clonotype)", choices, ignore.case=TRUE)]

  shiny::validate(shiny::need(length(choices) > 0, "No suitable metadata columns found for grouping."))

  selected_value <- if (!is.null(selected) && selected %in% choices) {
    selected
  } else if ("sample" %in% choices) {
    "sample"
  } else {
    choices[1]
  }
  
  updateSelectInput(session, "group_by", choices = choices, selected = selected_value)
}

# --- UIコンポーネント定義関数 (numericInputへ変更) ---

reductionInput <- function(ns) { selectInput(ns("reduction"), "Reduction", choices = NULL) }
groupByInput <- function(ns) { selectInput(ns("group_by"), "Group by", choices = NULL) }
pointSizeInput <- function(ns) { numericInput(ns("point_size"), "Point Size", value = 0.1, min = 0, max = 10, step = 0.1) }
labelSizeInput <- function(ns) { numericInput(ns("label_size"), "Label Size", value = 4, min = 0, max = 20, step = 1) }

commonPlotOptions <- function(ns) {
  tagList(
    selectInput(ns("legend"), "Legend Position", choices = c("right", "left", "top", "bottom", "none"), selected = "right"),
    numericInput(ns("plot_width"), "Plot Width (px)", value = 600, min = 200, max = 2000, step = 50),
    numericInput(ns("plot_height"), "Plot Height (px)", value = 500, min = 200, max = 2000, step = 50)
  )
}