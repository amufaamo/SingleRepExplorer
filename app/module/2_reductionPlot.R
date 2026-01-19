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
      )
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
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
      reduction_prefix <- sub("_$", "", names(dimnames(reduction_data))[1]) # UMAP_ or PC_ -> UMAP or PC
      table_data %>% select(all_of(input$group_by), starts_with(paste0(reduction_prefix, "_")))
    })
    
    output$coordinates_table <- DT::renderDataTable({ coordinates_data() })
    
    output$download_plot <- downloadHandler(
      filename = function() { paste0("DimPlot_", input$reduction, "_", inputg$group_by, ".pdf") },
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
  choices <- choices[!grepl("^(RNA_snn_res\\.|TCR|BCR)", choices)]

  validate(need(length(choices) > 0, "No suitable metadata columns found for grouping."))

  selected_value <- if (!is.null(selected) && selected %in% choices) {
    selected
  } else if ("sample" %in% choices) {
    "sample"
  } else {
    choices[1]
  }
  
  updateSelectInput(session$ns("group_by"), choices = choices, selected = selected_value)
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