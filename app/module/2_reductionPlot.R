# Required libraries are assumed to be loaded in app.R

#=================================================
# UI部分
#=================================================
reductionPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      useShinyjs(),

      # ---- Plot Options ----
      h4("Plot Options"),
      reductionInput(ns),
      groupByInput(ns),
      pointSizeInput(ns),
      labelSizeInput(ns),
      commonPlotOptions(ns),

      hr(),

      # ==============================
      # Annotation Section
      # ==============================
      h4(icon("tags"), " Annotation", style = "font-weight: 600;"),
      tabsetPanel(
        id = ns("annot_tab"),

        # ---- Auto Annotation Tab ----
        tabPanel("Auto",
          br(),
          radioButtons(ns("annotation_method"), "Method:",
            choices = c(
              "scType (Marker-based)"     = "sctype",
              "SingleR (Reference-based)" = "singler"
            ),
            selected = "sctype"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'sctype'", ns("annotation_method")),
            selectInput(ns("sctype_tissue"), "Tissue Type",
              choices = c(
                "Immune system", "Pancreas", "Heart", "Brain",
                "Placenta", "Kidney", "Liver", "Lung", "Muscle", "Intestine"
              ),
              selected = "Immune system"
            )
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'singler'", ns("annotation_method")),
            selectInput(ns("singler_ref"), "Reference Dataset",
              choices = list(
                "Human Primary Cell Atlas (Broad)" = "HumanPrimaryCellAtlasData",
                "Monaco Immune Data (Immune)"      = "MonacoImmuneData",
                "DICE (Immune)"                    = "DatabaseImmuneCellExpressionData",
                "Novershtern Hematopoietic"        = "NovershternHematopoieticData",
                "Blueprint Encode"                 = "BlueprintEncodeData",
                "Mouse RNA-seq (Broad)"            = "MouseRNAseqData",
                "ImmGen (Mouse Immune)"            = "ImmGenData"
              )
            )
          ),
          br(),
          actionButton(ns("run_annotation_btn"), "Run Annotation",
            icon = icon("tags"), class = "btn-info",
            style = "width: 100%;"
          )
        ),

        # ---- Manual Annotation Tab ----
        tabPanel("Manual",
          br(),
          tabsetPanel(
            id = ns("manual_annot_tab"),

            # -- Sub-tab ①: Rename Labels (In-place) --
            tabPanel("\u2460 Rename Labels",
              br(),
              div(class = "alert alert-info",
                style = "padding: 8px; font-size: 0.82em; margin-bottom: 10px;",
                icon("info-circle"),
                " Rename labels within an existing metadata column (in-place)."
              ),
              selectInput(ns("rename_col"), "Column to Edit", choices = NULL),
              uiOutput(ns("rename_labels_ui")),
              br(),
              actionButton(ns("rename_inplace_btn"), "Apply Rename",
                icon = icon("check"), class = "btn-warning",
                style = "width: 100%;"
              )
            ),

            # -- Sub-tab ②: Create New Group --
            tabPanel("\u2461 Create New Group",
              br(),
              div(class = "alert alert-info",
                style = "padding: 8px; font-size: 0.82em; margin-bottom: 10px;",
                icon("copy"),
                " Copy a column with new label names into a new metadata column."
              ),
              selectInput(ns("newgroup_source_col"), "Source Column", choices = NULL),
              textInput(ns("newgroup_col_name"), "New Column Name",
                value = "", placeholder = "e.g., CellType_manual"
              ),
              uiOutput(ns("newgroup_remap_ui")),
              br(),
              actionButton(ns("newgroup_save_btn"), "Save as New Column",
                icon = icon("save"), class = "btn-success",
                style = "width: 100%;"
              )
            ),

            # -- Sub-tab ③: Lasso Select --
            tabPanel("\u2462 Lasso Select",
              br(),
              div(class = "alert alert-info",
                style = "padding: 8px; font-size: 0.82em; margin-bottom: 10px;",
                icon("hand-pointer"),
                HTML(" Select cells <b>directly on the main UMAP</b> using the lasso tool.")
              ),
              tags$ol(
                style = "font-size: 0.82em; color: #555; padding-left: 18px; margin-bottom: 10px;",
                tags$li("The main plot switches to interactive mode automatically"),
                tags$li("Click the ", tags$b("lasso icon (\u29be)"), " in the plot toolbar (top-right)"),
                tags$li("Draw a free-hand selection around target cells"),
                tags$li("Set target column and label below, then click ", tags$b("Assign"))
              ),
              radioButtons(ns("lasso_col_mode"), "Target Column:",
                choices = c(
                  "Create new column" = "new",
                  "Assign to existing column" = "existing"
                ),
                selected = "new"
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'new'", ns("lasso_col_mode")),
                textInput(ns("lasso_new_col_name"), "New Column Name",
                  value = "LassoGroup",
                  placeholder = "e.g., manual_annotation"
                )
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'existing'", ns("lasso_col_mode")),
                selectInput(ns("lasso_target_col"), "Existing Column", choices = NULL)
              ),
              textInput(ns("lasso_label"), "Label to Assign",
                value = "Group_1", placeholder = "e.g., CD8 T cell"
              ),
              actionButton(ns("lasso_assign_btn"), "Assign Selected Cells",
                icon = icon("check-circle"), class = "btn-success",
                style = "width: 100%;"
              )
            )
          ) # end nested tabsetPanel (manual_annot_tab)
        ) # end tabPanel "Manual"
      ), # end outer tabsetPanel (annot_tab)

      hr(),

      # ==============================
      # Subset & Re-cluster Section
      # ==============================
      h4(icon("scissors"), " Subset & Re-cluster", style = "font-weight: 600;"),

      selectInput(ns("subset_source_col"), "Column", choices = NULL),
      selectizeInput(ns("subset_selected_groups"), "Groups to keep:",
        choices = NULL, multiple = TRUE,
        options = list(
          placeholder = "Select groups to keep...",
          plugins = list("remove_button")
        )
      ),
      div(
        style = "margin-top: -8px; margin-bottom: 6px;",
        actionLink(ns("subset_select_all"),   "Select all",   style = "font-size: 0.82em;"),
        tags$span(" | ", style = "font-size: 0.82em; color: #aaa;"),
        actionLink(ns("subset_deselect_all"), "Deselect all", style = "font-size: 0.82em;")
      ),

      checkboxInput(ns("recluster_check"), "Re-run Clustering (Normalization, PCA, UMAP, t-SNE)", value = TRUE),

      fluidRow(
        column(6,
          actionButton(ns("do_subset_btn"), "Subset",
            icon = icon("cut"), class = "btn-warning btn-sm",
            style = "width: 100%; margin-top: 4px;"
          )
        ),
        column(6,
          actionButton(ns("reset_subset_btn2"), "Reset Subset",
            icon = icon("undo"), class = "btn-secondary btn-sm",
            style = "width: 100%; margin-top: 4px;"
          )
        )
      ),

      hr(),

      # ---- Search & Highlight ----
      h4("Search & Highlight"),
      checkboxInput(ns("clonotype_search_toggle"), "Enable Highlight by Search", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s'] == true", ns("clonotype_search_toggle")),
        selectInput(ns("search_source"), "1. Select Data Source", choices = NULL),
        selectInput(ns("search_column"), "2. Select Column to Search", choices = NULL),
        selectizeInput(ns("search_values"), "3. Select Value(s) to Highlight",
          choices = NULL, multiple = TRUE,
          options = list(placeholder = "Type or select values")
        ),
        actionButton(ns("search_highlight_btn"), "Apply Highlight", icon = icon("search"), style = "margin-bottom: 5px; width: 100%;"),
        actionButton(ns("subset_highlight_btn"), "Subset Highlighted Cells", icon = icon("cut"), class = "btn-warning", style = "width: 100%;")
      )
    ), # end sidebarPanel

    mainPanel(
      h3("Plot"),
      div(style = "display: flex; justify-content: space-between; align-items: center;",
        downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
        div(
          actionButton(ns("reset_subset_btn"), "Reset Subset",
            icon = icon("undo"), class = "btn-secondary"
          ),
          actionButton(ns("subset_btn"), "Subset Selected Cells",
            icon = icon("cut"), class = "btn-warning"
          )
        )
      ),
      br(),
      uiOutput(ns("dynamic_plot_ui"))
    ) # end mainPanel
  ) # end sidebarLayout
}


#=================================================
# Server部分
#=================================================
reductionPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    highlighted_cells        <- reactiveVal(NULL)
    subset_barcodes_pending  <- reactiveVal(NULL)

    # ---- Helper: categorical metadata columns ----
    categorical_meta_cols <- reactive({
      req(myReactives$seurat_object)
      meta <- myReactives$seurat_object@meta.data
      is_cat <- sapply(names(meta), function(col) {
        is.character(meta[[col]]) || is.factor(meta[[col]])
      })
      cols <- names(meta)[is_cat]
      cols[!grepl(
        "^(orig\\.ident|barcode|RNA_snn_res|TCR_CT|BCR_CT|TCR_clonotype_id|BCR_clonotype_id|raw_clonotype|TCR_TRA|TCR_TRB|BCR_IGH|BCR_IGL|BCR_IGK)",
        cols, ignore.case = TRUE
      )]
    })

    # ---- Helper: level names for a column ----
    get_col_levels <- function(so, col) {
      if (!col %in% colnames(so@meta.data)) return(character(0))
      vals <- so@meta.data[[col]]
      if (is.factor(vals)) levels(vals) else sort(unique(as.character(vals)))
    }

    # ---- Update all dropdowns when Seurat object loads ----
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object

      # TCR/BCR status column
      if (!"TCR_BCR_status" %in% names(so@meta.data)) {
        clonotype_cols <- c("CTaa", "CTnt", "raw_clonotype_id", "clonotype_id", "clonotype")
        if (any(clonotype_cols %in% names(so@meta.data))) {
          col_to_use <- clonotype_cols[clonotype_cols %in% names(so@meta.data)][1]
          so$TCR_BCR_status <- ifelse(
            is.na(so@meta.data[[col_to_use]]) | so@meta.data[[col_to_use]] == "None",
            "Absent", "Present"
          )
          so$TCR_BCR_status <- factor(so$TCR_BCR_status, levels = c("Present", "Absent"))
          myReactives$seurat_object <- so
        }
      }

      update_reduction_choices(session, myReactives)
      update_group_by_select_input(session, myReactives)

      # Search source
      search_sources <- c("Seurat Metadata" = "metadata")
      if (!is.null(myReactives$tcr_df)) search_sources <- c(search_sources, "TCR Clonotypes" = "tcr_df")
      if (!is.null(myReactives$bcr_df)) search_sources <- c(search_sources, "BCR Clonotypes" = "bcr_df")
      updateSelectInput(session, "search_source", choices = search_sources)

    }, ignoreNULL = TRUE)

    # Update annotation-related dropdowns whenever categorical_meta_cols changes
    observe({
      cols <- categorical_meta_cols()
      current_group <- isolate(input$group_by)
      default_col   <- if (!is.null(current_group) && current_group %in% cols) current_group else cols[1]

      updateSelectInput(session, "rename_col",          choices = cols, selected = default_col)
      updateSelectInput(session, "newgroup_source_col", choices = cols, selected = default_col)
      updateSelectInput(session, "lasso_target_col",    choices = cols, selected = default_col)
      updateSelectInput(session, "subset_source_col",   choices = cols, selected = default_col)
    })

    # ================================================================
    # Subset & Re-cluster — group choices update
    # ================================================================

    # Repopulate group checkboxes when column or seurat object changes
    observe({
      req(myReactives$seurat_object, input$subset_source_col)
      lvls <- get_col_levels(myReactives$seurat_object, input$subset_source_col)
      # Keep current selection if still valid, otherwise select all
      current <- isolate(input$subset_selected_groups)
      selected <- if (length(current) > 0 && all(current %in% lvls)) current else lvls
      updateSelectizeInput(session, "subset_selected_groups",
        choices = lvls, selected = selected, server = TRUE
      )
    })

    observeEvent(input$subset_select_all, {
      req(myReactives$seurat_object, input$subset_source_col)
      lvls <- get_col_levels(myReactives$seurat_object, input$subset_source_col)
      updateSelectizeInput(session, "subset_selected_groups",
        choices = lvls, selected = lvls, server = TRUE
      )
    })

    observeEvent(input$subset_deselect_all, {
      req(myReactives$seurat_object, input$subset_source_col)
      lvls <- get_col_levels(myReactives$seurat_object, input$subset_source_col)
      updateSelectizeInput(session, "subset_selected_groups",
        choices = lvls, selected = character(0), server = TRUE
      )
    })

    # ================================================================
    # Subset & Re-cluster — main logic
    # ================================================================

    observeEvent(input$do_subset_btn, {
      req(myReactives$seurat_object, input$subset_source_col)
      so   <- myReactives$seurat_object
      grps <- input$subset_selected_groups

      if (is.null(grps) || length(grps) == 0) {
        showNotification("Please select at least one group to keep.", type = "error")
        return()
      }

      col_data          <- as.character(so@meta.data[[input$subset_source_col]])
      selected_barcodes <- rownames(so@meta.data)[col_data %in% grps]

      if (length(selected_barcodes) == 0) {
        showNotification("No cells matched the selection.", type = "error")
        return()
      }

      subset_barcodes_pending(selected_barcodes)

      recluster_msg <- if (isTRUE(input$recluster_check)) {
        "<br><b>Re-clustering</b> will be performed after subsetting."
      } else {
        "<br>Existing UMAP / cluster assignments will be kept."
      }

      showModal(modalDialog(
        title = tags$span(icon("cut"), " Confirm Subset & Re-cluster"),
        HTML(paste0(
          "Subset to <b>", length(selected_barcodes), "</b> cells from group(s): <b>",
          paste(grps, collapse = ", "), "</b>.",
          recluster_msg,
          "<br><br>All analysis tabs will update with the subsetted data."
        )),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_new_subset_btn"), "Subset",
            icon = icon("cut"), class = "btn-danger"
          )
        )
      ))
    })

    observeEvent(input$confirm_new_subset_btn, {
      removeModal()
      selected_barcodes <- subset_barcodes_pending()
      if (is.null(selected_barcodes) || length(selected_barcodes) == 0) return()

      so <- myReactives$seurat_object

      # Save full object for reset
      if (is.null(myReactives$full_seurat_object)) {
        myReactives$full_seurat_object <- so
      }

      tryCatch({
        so_sub   <- subset(so, cells = selected_barcodes)
        dims_use <- min(10L, ncol(so_sub) - 1L)
        dims_use <- max(dims_use, 2L)

        if (isTRUE(input$recluster_check)) {
          withProgress(message = "Subsetting & Re-clustering...", value = 0, {
            incProgress(0.15, detail = "Normalizing data")
            so_sub <- NormalizeData(so_sub, verbose = FALSE)
            incProgress(0.30, detail = "Finding variable features")
            so_sub <- FindVariableFeatures(so_sub, verbose = FALSE)
            incProgress(0.45, detail = "Scaling data")
            so_sub <- ScaleData(so_sub, verbose = FALSE)
            incProgress(0.58, detail = "PCA")
            so_sub <- RunPCA(so_sub, verbose = FALSE)
            incProgress(0.70, detail = "Finding neighbors")
            so_sub <- FindNeighbors(so_sub, dims = 1:dims_use, verbose = FALSE)
            incProgress(0.80, detail = "Clustering")
            so_sub <- FindClusters(so_sub, resolution = 0.5, verbose = FALSE)
            incProgress(0.85, detail = "UMAP")
            so_sub <- RunUMAP(so_sub, dims = 1:dims_use, verbose = FALSE)
            incProgress(0.95, detail = "t-SNE")
            so_sub <- RunTSNE(so_sub, dims = 1:dims_use, verbose = FALSE)
          })
          update_group_by_select_input(session, myReactives, selected = "seurat_clusters")
        }

        myReactives$seurat_object <- so_sub
        subset_barcodes_pending(NULL)

        showNotification(
          paste0("\u2705 Subsetted to ", ncol(so_sub), " cells.",
                 if (isTRUE(input$recluster_check)) " Re-clustering complete." else ""),
          type = "message", duration = 5
        )
      }, error = function(e) {
        showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 8)
      })
    }, ignoreInit = TRUE)

    observeEvent(input$reset_subset_btn2, {
      if (!is.null(myReactives$full_seurat_object)) {
        myReactives$seurat_object      <- myReactives$full_seurat_object
        myReactives$full_seurat_object <- NULL
        showNotification("\u2705 Dataset restored to full state.", type = "message", duration = 4)
      } else {
        showNotification("No subset has been performed yet.", type = "info", duration = 4)
      }
    })

    # ================================================================
    # Search & Highlight
    # ================================================================
    observeEvent(input$clonotype_search_toggle, {
      if (!isTRUE(input$clonotype_search_toggle)) highlighted_cells(NULL)
    })

    # Re-run whenever search_source changes OR seurat_object is updated
    # (e.g. new columns added via Lasso Select)
    observe({
      req(input$search_source)
      choices <- if (input$search_source == "metadata") {
        req(myReactives$seurat_object)
        meta    <- myReactives$seurat_object@meta.data
        is_cat  <- sapply(names(meta), function(col) {
          is.character(meta[[col]]) || is.factor(meta[[col]]) || is.logical(meta[[col]])
        })
        names(meta)[is_cat]
      } else if (!is.null(myReactives[[input$search_source]])) {
        names(myReactives[[input$search_source]])
      } else NULL
      updateSelectInput(session, "search_column", choices = choices)
    })

    observeEvent(input$search_column, {
      req(input$search_source, input$search_column, nzchar(input$search_column))
      unique_vals <- if (input$search_source == "metadata") {
        req(myReactives$seurat_object)
        col_data <- myReactives$seurat_object@meta.data[[input$search_column]]
        if (is.logical(col_data)) c("TRUE", "FALSE") else unique(col_data)
      } else if (!is.null(myReactives[[input$search_source]])) {
        unique(myReactives[[input$search_source]][[input$search_column]])
      } else NULL
      updateSelectizeInput(session, "search_values", choices = sort(na.omit(unique_vals)), server = TRUE)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    observeEvent(input$search_highlight_btn, {
      req(input$search_source, input$search_column, input$search_values)
      barcodes_to_highlight <- if (input$search_source == "metadata") {
        meta_df <- myReactives$seurat_object@meta.data
        vals    <- if (is.logical(meta_df[[input$search_column]])) {
          as.logical(input$search_values)
        } else {
          input$search_values
        }
        rownames(meta_df)[meta_df[[input$search_column]] %in% vals]
      } else {
        myReactives[[input$search_source]] %>%
          filter(.data[[input$search_column]] %in% input$search_values) %>%
          pull(barcode)
      }
      highlighted_cells(intersect(barcodes_to_highlight, colnames(myReactives$seurat_object)))
    })

    # ================================================================
    # Main DimPlot
    # ================================================================
    plot <- reactive({
      req(myReactives$seurat_object, input$reduction)
      so <- myReactives$seurat_object

      if (isTRUE(input$clonotype_search_toggle) &&
          !is.null(highlighted_cells()) && length(highlighted_cells()) > 0) {
        p <- DimPlot(so,
          reduction       = input$reduction,
          cells.highlight = highlighted_cells(),
          cols.highlight  = "darkblue",
          pt.size         = input$point_size,
          label           = FALSE
        ) + guides(color = "none") + theme(legend.position = input$legend)
      } else {
        req(input$group_by)
        p <- DimPlot(so,
          reduction  = input$reduction,
          group.by   = input$group_by,
          label      = TRUE,
          pt.size    = input$point_size,
          label.size = input$label_size
        ) + theme(legend.position = input$legend)
      }
      p
    })

    output$dynamic_plot_ui <- renderUI({
      req(input$plot_width, input$plot_height)
      # Switch main plot to interactive plotly when Lasso Select tab is active
      is_lasso <- isTRUE(input$annot_tab == "Manual") &&
                  isTRUE(input$manual_annot_tab == "\u2462 Lasso Select")

      if (is_lasso) {
        tagList(
          div(
            class = "alert alert-warning",
            style = "padding: 8px 14px; font-size: 0.85em; margin-bottom: 8px; display: flex; align-items: center; gap: 8px;",
            icon("hand-pointer"),
            HTML("<b>Lasso Mode:</b>&nbsp; Click the lasso icon in the plot toolbar, draw a selection, then click <b>Assign Selected Cells</b> in the sidebar.")
          ),
          plotlyOutput(session$ns("main_lasso_plot"),
            width  = paste0(input$plot_width,  "px"),
            height = paste0(input$plot_height, "px")
          )
        )
      } else {
        plotOutput(session$ns("plot"),
          width  = paste0(input$plot_width,  "px"),
          height = paste0(input$plot_height, "px"),
          brush  = brushOpts(id = session$ns("plot_brush"), resetOnNew = TRUE)
        )
      }
    })

    output$plot <- renderPlot({ plot() })

    # ================================================================
    # Manual Annotation — Sub-tab ①: Rename Labels (In-place)
    # ================================================================

    # Dynamic text inputs (one per level)
    output$rename_labels_ui <- renderUI({
      req(myReactives$seurat_object, input$rename_col)
      so  <- myReactives$seurat_object
      lvls <- get_col_levels(so, input$rename_col)
      if (length(lvls) == 0) return(p("No labels found.", style = "color:gray;"))

      safe_ids <- make.names(lvls, unique = TRUE)
      tagList(
        p(em(paste(length(lvls), "labels:")),
          style = "font-size: 0.82em; margin: 6px 0 4px;"),
        lapply(seq_along(lvls), function(i) {
          div(style = "display:flex; align-items:center; gap:6px; margin-bottom:3px;",
            div(style = paste0(
              "width:40%; font-size:0.80em; color:#555;",
              "overflow:hidden; text-overflow:ellipsis; white-space:nowrap;"
            ), icon("arrow-right"), " ", lvls[i]),
            div(style = "width:60%;",
              textInput(ns(paste0("inplace_rename_", safe_ids[i])),
                label = NULL, value = lvls[i]
              )
            )
          )
        })
      )
    })

    observeEvent(input$rename_inplace_btn, {
      req(myReactives$seurat_object, input$rename_col)
      so   <- myReactives$seurat_object
      col  <- input$rename_col
      lvls <- get_col_levels(so, col)
      if (length(lvls) == 0) return()

      safe_ids  <- make.names(lvls, unique = TRUE)
      name_map  <- setNames(
        sapply(seq_along(lvls), function(i) {
          val <- trimws(input[[paste0("inplace_rename_", safe_ids[i])]])
          if (is.null(val) || !nzchar(val)) lvls[i] else val
        }),
        lvls
      )

      char_vec          <- as.character(so@meta.data[[col]])
      so@meta.data[[col]] <- factor(unname(name_map[char_vec]))
      myReactives$seurat_object <- so

      update_group_by_select_input(session, myReactives, selected = col)
      showNotification(paste("✅ Labels renamed in-place in '", col, "'"),
        type = "message", duration = 4)
    })

    # ================================================================
    # Manual Annotation — Sub-tab ②: Create New Group
    # ================================================================

    output$newgroup_remap_ui <- renderUI({
      req(myReactives$seurat_object, input$newgroup_source_col)
      so   <- myReactives$seurat_object
      lvls <- get_col_levels(so, input$newgroup_source_col)
      if (length(lvls) == 0) return(p("No labels found.", style = "color:gray;"))

      safe_ids <- make.names(lvls, unique = TRUE)
      tagList(
        p(em(paste(length(lvls), "labels to map:")),
          style = "font-size: 0.82em; margin: 6px 0 4px;"),
        lapply(seq_along(lvls), function(i) {
          div(style = "display:flex; align-items:center; gap:6px; margin-bottom:3px;",
            div(style = paste0(
              "width:40%; font-size:0.80em; color:#555;",
              "overflow:hidden; text-overflow:ellipsis; white-space:nowrap;"
            ), icon("arrow-right"), " ", lvls[i]),
            div(style = "width:60%;",
              textInput(ns(paste0("newgroup_remap_", safe_ids[i])),
                label = NULL, value = lvls[i]
              )
            )
          )
        })
      )
    })

    observeEvent(input$newgroup_save_btn, {
      req(myReactives$seurat_object, input$newgroup_source_col)
      new_col <- gsub("[^A-Za-z0-9_.]", "_", trimws(input$newgroup_col_name))
      if (!nzchar(new_col)) {
        showNotification("Please enter a valid new column name.", type = "error")
        return()
      }

      so   <- myReactives$seurat_object
      col  <- input$newgroup_source_col
      lvls <- get_col_levels(so, col)
      if (length(lvls) == 0) return()

      safe_ids <- make.names(lvls, unique = TRUE)
      name_map <- setNames(
        sapply(seq_along(lvls), function(i) {
          val <- trimws(input[[paste0("newgroup_remap_", safe_ids[i])]])
          if (is.null(val) || !nzchar(val)) lvls[i] else val
        }),
        lvls
      )

      char_vec              <- as.character(so@meta.data[[col]])
      so@meta.data[[new_col]] <- factor(unname(name_map[char_vec]))
      myReactives$seurat_object <- so

      update_group_by_select_input(session, myReactives, selected = new_col)
      showNotification(paste("✅ New annotation column '", new_col, "' created."),
        type = "message", duration = 4)
    })

    # ================================================================
    # Manual Annotation — Sub-tab ③: Lasso Select
    # Rendered as the main plot when lasso tab is active (see dynamic_plot_ui)
    # ================================================================

    output$main_lasso_plot <- renderPlotly({
      req(myReactives$seurat_object, input$reduction)
      so     <- myReactives$seurat_object
      embeds <- Embeddings(so, reduction = input$reduction)

      plot_df <- as.data.frame(embeds[, 1:2])
      colnames(plot_df) <- c("dim1", "dim2")
      plot_df$barcode   <- rownames(plot_df)

      is_highlighted <- isTRUE(input$clonotype_search_toggle) && !is.null(highlighted_cells()) && length(highlighted_cells()) > 0

      if (is_highlighted) {
        plot_df$color_var <- ifelse(plot_df$barcode %in% highlighted_cells(), "Highlighted", "Unselected")
        plot_df$color_var <- factor(plot_df$color_var, levels = c("Unselected", "Highlighted"))
        plot_df <- plot_df[order(plot_df$color_var), ]
        
        colors <- c("Unselected" = "#e0e0e0", "Highlighted" = "darkblue")
        
        plot_ly(
          data      = plot_df,
          x         = ~dim1,
          y         = ~dim2,
          color     = ~color_var,
          colors    = colors,
          type      = "scattergl",
          mode      = "markers",
          marker    = list(size = 4, opacity = 0.8),
          key       = ~barcode,
          source    = "main_lasso_plot",
          text      = ~paste0("<b>Status:</b> ", color_var, "<br><b>Barcode:</b> ", barcode),
          hoverinfo = "text"
        ) %>% layout(dragmode = "lasso", xaxis = list(title = colnames(embeds)[1], zeroline = FALSE), yaxis = list(title = colnames(embeds)[2], zeroline = FALSE), legend = list(itemsizing = "constant"), margin = list(l = 50, r = 20, t = 30, b = 50)) %>% event_register("plotly_selected")
      } else {
        color_col <- if (!is.null(input$group_by) && input$group_by %in% colnames(so@meta.data)) {
          as.character(so@meta.data[plot_df$barcode, input$group_by])
        } else {
          rep("Cells", nrow(plot_df))
        }
        plot_df$color_var <- color_col

        plot_ly(
          data      = plot_df,
          x         = ~dim1,
          y         = ~dim2,
          color     = ~color_var,
          type      = "scattergl",
          mode      = "markers",
          marker    = list(size = 3, opacity = 0.75),
          key       = ~barcode,
          source    = "main_lasso_plot",
          text      = ~paste0("<b>", if (!is.null(input$group_by)) input$group_by else "Group", ":</b> ", color_var, "<br><b>Barcode:</b> ", barcode),
          hoverinfo = "text"
        ) %>% layout(dragmode = "lasso", xaxis = list(title = colnames(embeds)[1], zeroline = FALSE), yaxis = list(title = colnames(embeds)[2], zeroline = FALSE), legend = list(itemsizing = "constant"), margin = list(l = 50, r = 20, t = 30, b = 50)) %>% event_register("plotly_selected")
      }
    })

    observeEvent(input$lasso_assign_btn, {
      req(myReactives$seurat_object, input$lasso_label)

      sel <- plotly::event_data("plotly_selected", source = "main_lasso_plot")
      if (is.null(sel) || length(sel$key) == 0) {
        showNotification(
          "No cells selected. Use the lasso tool in the main UMAP plot to draw a selection first.",
          type = "warning", duration = 6
        )
        return()
      }

      selected_barcodes <- sel$key
      label_name        <- trimws(input$lasso_label)
      if (!nzchar(label_name)) label_name <- "Group_1"

      so <- myReactives$seurat_object

      # Determine target column — new or existing
      if (isTRUE(input$lasso_col_mode == "new")) {
        target_col <- gsub("[^A-Za-z0-9_.]", "_", trimws(input$lasso_new_col_name))
        if (!nzchar(target_col)) target_col <- "LassoGroup"
      } else {
        req(input$lasso_target_col)
        target_col <- input$lasso_target_col
      }

      # Initialize column with "Unassigned" if new
      if (!target_col %in% colnames(so@meta.data)) {
        so@meta.data[[target_col]] <- "Unassigned"
      }

      current_vec <- as.character(so@meta.data[[target_col]])
      current_vec[rownames(so@meta.data) %in% selected_barcodes] <- label_name
      so@meta.data[[target_col]] <- factor(current_vec)

      myReactives$seurat_object <- so
      update_group_by_select_input(session, myReactives, selected = target_col)

      showNotification(
        paste("✅", length(selected_barcodes), "cells labeled as '",
              label_name, "' in '", target_col, "'"),
        type = "message", duration = 5
      )
    })

    # ================================================================
    # Auto Annotation
    # ================================================================
    observeEvent(input$run_annotation_btn, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object

      if (!exists("run_sctype_and_update_seurat")) {
        if (require(shinyalert)) shinyalert("Error", "Annotation functions not found.", type = "error")
        return()
      }

      tryCatch({
        withProgress(message = "Running Annotation...", value = 0, {
          if (input$annotation_method == "sctype") {
            incProgress(0.2, detail = "Running scType cell typing...")
            so <- run_sctype_and_update_seurat(so, tissue = input$sctype_tissue)
          } else if (input$annotation_method == "singler") {
            incProgress(0.1, detail = "Preparing SingleR (memory-efficient mode)...")
            # Downsample for large objects to avoid memory crash
            n_cells <- ncol(so)
            if (n_cells > 10000) {
              showNotification(paste0("Large dataset (", n_cells, " cells). Downsampling to 10,000 for annotation."), type = "warning", duration = 8)
              set.seed(42)
              sampled_cells <- sample(colnames(so), 10000)
              so_small <- subset(so, cells = sampled_cells)
            } else {
              so_small <- so
            }
            incProgress(0.3, detail = "Running SingleR cell typing...")
            so_small <- run_singler_and_update_seurat(so_small, ref_name = input$singler_ref)
            # Transfer labels back to full object
            annot_col <- "singler_celltype"
            if (annot_col %in% names(so_small@meta.data)) {
              labels_full <- rep(NA_character_, ncol(so))
              names(labels_full) <- colnames(so)
              labels_full[colnames(so_small)] <- as.character(so_small@meta.data[[annot_col]])
              so <- AddMetaData(so, labels_full, col.name = annot_col)
            }
          }
          myReactives$seurat_object <- so
          incProgress(1, detail = "Done.")
        })

        selected_col <- if (input$annotation_method == "sctype") "sctype_celltype" else "singler_celltype"
        update_group_by_select_input(session, myReactives, selected = selected_col)

        if (require(shinyalert)) shinyalert("Success", "Annotation completed.", type = "success")
      }, error = function(e) {
        showNotification(paste("Annotation error:", conditionMessage(e)), type = "error", duration = 15)
      })
    })

    # ================================================================
    # Subsetting (Brush on main DimPlot)
    # ================================================================
    coordinates_data <- reactive({
      req(myReactives$seurat_object, input$reduction, input$group_by)
      so             <- myReactives$seurat_object
      reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
      table_data     <- cbind(so@meta.data, reduction_data)
      table_data$Barcode <- rownames(table_data)

      reduction_prefix <- if (grepl("_1$", colnames(reduction_data)[1])) {
        sub("_1$", "", colnames(reduction_data)[1])
      } else {
        "UMAP"
      }
      table_data %>%
        select(Barcode, all_of(input$group_by),
               starts_with(paste0(reduction_prefix, "_")))
    })

    observeEvent(input$subset_btn, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      
      # Determine if we are using Lasso (plotly) or Brush (ggplot)
      is_lasso <- isTRUE(input$annot_tab == "Manual") && isTRUE(input$manual_annot_tab == "\u2462 Lasso Select")

      if (is_lasso) {
        sel <- plotly::event_data("plotly_selected", source = "main_lasso_plot")
        if (is.null(sel) || length(sel$key) == 0) {
          if (require(shinyalert)) shinyalert("Warning", "No cells selected with Lasso. Draw a selection first.", type = "warning")
          return()
        }
        selected_barcodes <- sel$key
      } else {
        req(input$plot_brush)
        coords <- coordinates_data()
        coord_cols <- grep("_[12]$", colnames(coords), value = TRUE)
        if (length(coord_cols) < 2) return()
        selected_rows <- brushedPoints(coords, input$plot_brush, xvar = coord_cols[1], yvar = coord_cols[2])
        if (nrow(selected_rows) == 0) {
          if (require(shinyalert)) shinyalert("Warning", "No cells selected. Please draw a box on the UMAP plot.", type = "warning")
          return()
        }
        selected_barcodes <- selected_rows$Barcode
      }

      showModal(modalDialog(
        title = "Confirm Subset",
        HTML(paste0(
          "You selected <b>", length(selected_barcodes), "</b> cells.<br>",
          "Keep only these cells? <b>All tabs will update immediately.</b>"
        )),
        br(), br(),
        checkboxInput(ns("rerun_umap_checkbox"),
          "Re-run Clustering (Normalization, PCA, UMAP, t-SNE) after subsetting", value = FALSE),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_subset_btn"), "Yes, Subset", class = "btn-danger")
        )
      ))
      
      # Store pending barcodes to subset if clicked
      subset_barcodes_pending(selected_barcodes)
    })
    
    observeEvent(input$subset_highlight_btn, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      selected_barcodes <- highlighted_cells()
      if (is.null(selected_barcodes) || length(selected_barcodes) == 0) {
        showNotification("No highlighted cells found. Please apply highlight first.", type = "error")
        return()
      }
      
      showModal(modalDialog(
        title = "Confirm Highlight Subset",
        HTML(paste0(
          "You have <b>", length(selected_barcodes), "</b> highlighted cells.<br>",
          "Keep only these cells? <b>All tabs will update immediately.</b>"
        )),
        br(), br(),
        checkboxInput(ns("rerun_umap_checkbox"),
          "Re-run Clustering (Normalization, PCA, UMAP, t-SNE) after subsetting", value = FALSE),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_subset_btn"), "Yes, Subset", class = "btn-danger")
        )
      ))
      subset_barcodes_pending(selected_barcodes)
    })

    observeEvent(input$confirm_subset_btn, {
      removeModal()
      selected_barcodes <- subset_barcodes_pending()
      if (is.null(selected_barcodes) || length(selected_barcodes) == 0) return()
      so <- myReactives$seurat_object
      
      if (is.null(myReactives$full_seurat_object)) {
        myReactives$full_seurat_object <- so
      }
      so_sub <- subset(so, cells = selected_barcodes)


        if (isTRUE(input$rerun_umap_checkbox)) {
          withProgress(message = "Re-running pipeline...", value = 0, {
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
            incProgress(0.85, detail = "UMAP")
            so_sub <- RunUMAP(so_sub, dims = 1:10)
            incProgress(0.95, detail = "t-SNE")
            so_sub <- RunTSNE(so_sub, dims = 1:10)
          })
        }

        myReactives$seurat_object <- so_sub
        if (require(shinyalert)) {
          shinyalert("Success",
            paste("Dataset subsetted to", ncol(so_sub), "cells."), type = "success")
        }
      }, ignoreInit = TRUE, once = TRUE)

    observeEvent(input$reset_subset_btn, {
      if (!is.null(myReactives$full_seurat_object)) {
        myReactives$seurat_object      <- myReactives$full_seurat_object
        myReactives$full_seurat_object <- NULL
        if (require(shinyalert)) {
          shinyalert("Reset", "Dataset returned to its full state.", type = "success")
        }
      } else {
        if (require(shinyalert)) {
          shinyalert("Info", "No subset has been performed yet.", type = "info")
        }
      }
    })

    # ================================================================
    # Downloads
    # ================================================================
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("DimPlot_", input$reduction, "_", input$group_by, ".pptx")
      },
      content = function(file) {
        req(plot())
        save_plot_as_pptx(file, plot(), input$plot_width, input$plot_height)
      }
    )

  })
}


# Helper functions moved to utils.R
