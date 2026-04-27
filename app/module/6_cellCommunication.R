cellCommunicationUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      useShinyjs(),
      h4("Cell-Cell Communication (CellChat)"),
      p("Infer signaling networks between cell clusters and compare networks across conditions."),
      if (!requireNamespace("CellChat", quietly = TRUE)) {
        div(class = "alert alert-warning",
          icon("exclamation-triangle"),
          strong(" CellChat not detected."),
          " Note: CellChat is included in the official Docker image (v1.0.0+).",
          tags$code("remotes::install_github('jinworks/CellChat')")
        )
      },

      # --- Analysis Mode ---
      radioButtons(ns("analysis_mode"), "Analysis Mode",
        choices = c("Single Condition" = "single", "Compare Two Conditions" = "compare"),
        selected = "single", inline = TRUE
      ),

      # --- Cluster annotation ---
      selectInput(ns("group_by"), "Cell Annotation (Clusters)", choices = NULL),

      # --- Compare mode: condition selectors ---
      conditionalPanel(
        condition = sprintf("input['%s'] == 'compare'", ns("analysis_mode")),
        selectInput(ns("condition_col"), "Condition Column", choices = NULL),
        selectInput(ns("condition_a"),   "Condition A (e.g. Baseline)",      choices = NULL),
        selectInput(ns("condition_b"),   "Condition B (e.g. Convalescent)",   choices = NULL)
      ),

      selectInput(ns("cellchat_db"), "Signaling Database",
        choices = c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "All"),
        selected = "Secreted Signaling"
      ),
      numericInput(ns("max_cells_per_cond"),
                   "Max cells per condition (memory cap)",
                   value = 3000, min = 500, max = 15000, step = 500),
      tags$p(tags$small(
        "Lower = less RAM & faster. Compare mode processes two conditions, so total RAM ≈ 2× this value. Try 2000 if you hit disconnects."
      ), style = "color: #666;"),
      actionButton(ns("run_cellchat"), "Run CellChat Inference",
                   class = "btn-primary", width = "100%", icon = icon("network-wired")),
      hr(),

      # --- Visualization ---
      h4("Visualization"),

      # Single mode plot type
      conditionalPanel(
        condition = sprintf("input['%s'] == 'single'", ns("analysis_mode")),
        selectInput(ns("plot_type"), "Plot Type",
          choices = c("Network Circle", "Bubble Plot"))
      ),

      # Compare mode plot type
      conditionalPanel(
        condition = sprintf("input['%s'] == 'compare'", ns("analysis_mode")),
        selectInput(ns("compare_plot_type"), "Plot Type",
          choices = c(
            "Pathway Strength Comparison" = "pathway_compare",
            "Differential Interaction Network" = "diff_network",
            "Circle Plot (per condition)"  = "circle_per_cond",
            "Bubble Plot (per condition)"  = "bubble_per_cond"
          )
        )
      ),

      # Pathway selector (circle plots)
      conditionalPanel(
        condition = sprintf(
          "(input['%s'] == 'single' && input['%s'] == 'Network Circle') || (input['%s'] == 'compare' && input['%s'] == 'circle_per_cond')",
          ns("analysis_mode"), ns("plot_type"), ns("analysis_mode"), ns("compare_plot_type")
        ),
        selectInput(ns("pathway_show"), "Signaling Pathway", choices = NULL)
      ),

      # Which condition to show (circle/bubble per condition)
      conditionalPanel(
        condition = sprintf(
          "input['%s'] == 'compare' && (input['%s'] == 'circle_per_cond' || input['%s'] == 'bubble_per_cond')",
          ns("analysis_mode"), ns("compare_plot_type"), ns("compare_plot_type")
        ),
        selectInput(ns("which_cond"), "Show Condition",
          choices = c("Condition A" = "A", "Condition B" = "B"))
      ),

      # Bubble plot: sources / targets
      conditionalPanel(
        condition = sprintf(
          "(input['%s'] == 'single' && input['%s'] == 'Bubble Plot') || (input['%s'] == 'compare' && input['%s'] == 'bubble_per_cond')",
          ns("analysis_mode"), ns("plot_type"), ns("analysis_mode"), ns("compare_plot_type")
        ),
        selectInput(ns("bubble_sources"), "Source Clusters (blank = all)",
          choices = NULL, multiple = TRUE),
        selectInput(ns("bubble_targets"), "Target Clusters (blank = all)",
          choices = NULL, multiple = TRUE)
      ),

      # Differential network measure
      conditionalPanel(
        condition = sprintf("input['%s'] == 'compare' && input['%s'] == 'diff_network'",
          ns("analysis_mode"), ns("compare_plot_type")),
        selectInput(ns("diff_measure"), "Measure",
          choices = c("Number of interactions" = "count", "Interaction strength" = "weight"),
          selected = "count")
      ),

      hr(),
      h5("Plot Options"),
      div(style = "display: flex; gap: 10px;",
        numericInput(ns("plot_width"), "Plot Width (px)", min = 100, max = 2000, value = 700, step = 100),
        numericInput(ns("plot_height"), "Plot Height (px)", min = 100, max = 2000, value = 700, step = 100)
      ),
      numericInput(ns("vertex_size_max"), "Max Vertex Size (Network)", value = 10, min = 2, max = 30, step = 1),
      numericInput(ns("edge_width_max"),  "Max Edge Width (Network)",  value = 8,  min = 1, max = 20, step = 1),
      hr()
    ),

    mainPanel(
      h3("Inferred Communication Networks"),
      p("Note: Running CellChat may take several minutes. Requires the official Docker container (v1.0.0+)."),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns("networkPlot"))
    )
  )
}

# ============================================================
# Helper: run CellChat pipeline on a single Seurat subset
# ============================================================
.log_mem <- function(tag) {
  mem_mb <- tryCatch(round(sum(gc(verbose = FALSE)[, 2]), 1), error = function(e) NA)
  message(sprintf("[CellChat][mem] %s: %s MB used", tag, mem_mb))
}

.slim_seurat_for_cellchat <- function(so_sub) {
  # Drop everything CellChat doesn't need: scale.data, SCT/integrated assays,
  # graphs, commands, reductions. Keep only RNA counts + metadata.
  tryCatch({
    if (length(so_sub@reductions)    > 0) so_sub@reductions    <- list()
    if (length(so_sub@graphs)        > 0) so_sub@graphs        <- list()
    if (length(so_sub@neighbors)     > 0) so_sub@neighbors     <- list()
    if (length(so_sub@commands)      > 0) so_sub@commands      <- list()
    # Drop non-RNA assays (SCT, integrated, etc.)
    other_assays <- setdiff(names(so_sub@assays), "RNA")
    for (a in other_assays) so_sub[[a]] <- NULL
    # Drop scale.data layers (keep counts, drop data if present and recoverable)
    if (inherits(tryCatch(so_sub[["RNA"]], error = function(e) NULL), "Assay5")) {
      existing_layers <- tryCatch(Seurat::Layers(so_sub[["RNA"]]), error = function(e) character(0))
      for (lyr in grep("^scale", existing_layers, value = TRUE)) {
        Seurat::LayerData(so_sub, assay = "RNA", layer = lyr) <- NULL
      }
    }
  }, error = function(e) message("[CellChat] slim warning: ", e$message))
  so_sub
}

.run_single_cellchat <- function(so_sub, group_by_col, db_name, max_cells = 3000) {
  # Label — prepend "Cluster " to numeric labels (CellChat rejects bare "0")
  so_sub@meta.data$cellchat_labels <- as.character(so_sub@meta.data[[group_by_col]])
  so_sub@meta.data$cellchat_labels[is.na(so_sub@meta.data$cellchat_labels) |
                                    !nzchar(so_sub@meta.data$cellchat_labels)] <- "Unknown"
  so_sub@meta.data$cellchat_labels <- ifelse(
    grepl("^[0-9]+$", so_sub@meta.data$cellchat_labels),
    paste0("Cluster ", so_sub@meta.data$cellchat_labels),
    so_sub@meta.data$cellchat_labels
  )

  # Downsample
  if (ncol(so_sub) > max_cells) {
    set.seed(42)
    so_sub <- so_sub[, sample(colnames(so_sub), max_cells)]
    message("[CellChat] Downsampled to ", max_cells, " cells")
  }

  # Remove small groups
  grp_counts <- table(so_sub@meta.data$cellchat_labels)
  small_grps <- names(grp_counts[grp_counts < 3])
  if (length(small_grps) > 0) {
    keep <- colnames(so_sub)[!so_sub@meta.data$cellchat_labels %in% small_grps]
    so_sub <- so_sub[, keep]
  }
  if (ncol(so_sub) == 0) stop("No cells remaining after filtering small groups.")

  # Assay5 → JoinLayers
  if (inherits(tryCatch(so_sub[["RNA"]], error = function(e) NULL), "Assay5")) {
    so_sub <- tryCatch(Seurat::JoinLayers(so_sub), error = function(e) so_sub)
  }

  # Slim down Seurat object to free memory before CellChat's heavy work
  so_sub <- .slim_seurat_for_cellchat(so_sub)
  gc(verbose = FALSE)
  .log_mem("after slim")

  message("[CellChat] Creating object | cells=", ncol(so_sub),
          " groups=", length(unique(so_sub@meta.data$cellchat_labels)))

  cc <- CellChat::createCellChat(object = so_sub, group.by = "cellchat_labels", assay = "RNA")

  # Free the intermediate Seurat subset as CellChat now has its own copy
  rm(so_sub); gc(verbose = FALSE)
  .log_mem("after createCellChat")

  # DB
  CellChatDB <- CellChat::CellChatDB.human
  cc@DB <- if (db_name == "All") CellChatDB else CellChat::subsetDB(CellChatDB, search = db_name)

  # Pipeline
  cc <- CellChat::subsetData(cc)

  # Use presto (fast Wilcoxon) if installed; otherwise fall back to standard Wilcoxon.
  has_presto <- requireNamespace("presto", quietly = TRUE)
  message("[CellChat] presto available: ", has_presto, " → do.fast=", has_presto)
  cc <- CellChat::identifyOverExpressedGenes(cc, do.fast = has_presto)
  .log_mem("after identifyOverExpressedGenes")

  cc <- CellChat::identifyOverExpressedInteractions(cc)
  message("[CellChat] LR pairs: ", tryCatch(nrow(cc@LR$LRsig), error = function(e) "N/A"))

  cc <- CellChat::computeCommunProb(cc, type = "triMean", raw.use = TRUE,
                                    nboot = 5, population.size = FALSE)
  .log_mem("after computeCommunProb")

  n_cells <- nrow(cc@meta)
  min_cells_val <- max(1, floor(n_cells / 1000))
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells_val)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  message("[CellChat] Pathways detected: ", length(cc@netP$pathways))
  gc(verbose = FALSE)
  cc
}

# ============================================================
# Server
# ============================================================
cellCommunicationServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cellchat_results <- reactiveVal(NULL)   # single mode
    cellchat_list    <- reactiveVal(NULL)   # compare mode: list(A=cc_a, B=cc_b)
    cellchat_merged  <- reactiveVal(NULL)   # compare mode: mergeCellChat result

    # --- Update selectors when Seurat object loads ---
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      so <- myReactives$seurat_object
      # Include: character, factor, OR numeric/integer with ≤100 unique values (cluster IDs)
      is_groupable <- sapply(so@meta.data, function(x) {
        is.character(x) || is.factor(x) ||
        ((is.numeric(x) || is.integer(x)) && length(unique(x[!is.na(x)])) <= 100)
      })
      choices <- names(so@meta.data)[is_groupable]
      choices <- choices[!grepl("^(RNA_snn|UMAP|PC_|barcode|orig.ident|cellchat_labels|nCount|nFeature|percent|log10)", choices)]
      # Ensure seurat_clusters is always first if present
      if ("seurat_clusters" %in% choices) {
        choices <- c("seurat_clusters", setdiff(choices, "seurat_clusters"))
      }
      message("[CellChat] group_by choices: ", paste(head(choices, 10), collapse=", "))
      updateSelectInput(session, "group_by",     choices = choices, selected = "seurat_clusters")
      updateSelectInput(session, "condition_col", choices = choices,
                        selected = if ("sample" %in% choices) "sample" else choices[1])
    })

    observeEvent(input$condition_col, {
      req(myReactives$seurat_object, input$condition_col)
      vals <- sort(na.omit(unique(as.character(
        myReactives$seurat_object@meta.data[[input$condition_col]]
      ))))
      updateSelectInput(session, "condition_a", choices = vals,
                        selected = vals[1])
      updateSelectInput(session, "condition_b", choices = vals,
                        selected = if (length(vals) >= 2) vals[2] else vals[1])
    })

    # --- Run button ---
    observeEvent(input$run_cellchat, {
      req(myReactives$seurat_object)
      if (!requireNamespace("CellChat", quietly = TRUE)) {
        shinyalert::shinyalert(title = "CellChat Not Available",
          text = "CellChat is not installed. Rebuild the Docker image (v2.0+).", type = "warning")
        return()
      }
      if (!require("CellChat", quietly = TRUE, character.only = TRUE)) {
        showNotification("CellChat failed to load.", type = "error"); return()
      }

      so <- myReactives$seurat_object

      # Resolve group_by
      group_by_col <- input$group_by
      if (is.null(group_by_col) || !nzchar(group_by_col) || !group_by_col %in% colnames(so@meta.data)) {
        group_by_col <- if ("seurat_clusters" %in% colnames(so@meta.data)) "seurat_clusters" else {
          all_cols <- colnames(so@meta.data)
          all_cols[sapply(so@meta.data[all_cols], function(x) is.character(x) || is.factor(x))][1]
        }
        message("[CellChat] group_by not found, falling back to: ", group_by_col)
      }
      message("[CellChat] group_by=", group_by_col,
              " | type=", class(so@meta.data[[group_by_col]]),
              " | n_unique=", length(unique(so@meta.data[[group_by_col]])),
              " | mode=", input$analysis_mode)

      max_cells_val <- max(500, as.integer(input$max_cells_per_cond %||% 3000))

      if (input$analysis_mode == "single") {
        # ---- Single condition mode ----
        withProgress(message = "Running CellChat (single)...", value = 0, {
          tryCatch({
            incProgress(0.1, detail = "Preparing...")
            cc <- .run_single_cellchat(so, group_by_col, input$cellchat_db, max_cells = max_cells_val)
            incProgress(0.9, detail = "Done")

            cellchat_results(cc)
            cellchat_list(NULL)
            cellchat_merged(NULL)

            pathways <- cc@netP$pathways
            updateSelectInput(session, "pathway_show",
              choices = if (length(pathways) > 0) pathways else "(none detected)")
            all_idents <- levels(cc@idents)
            updateSelectInput(session, "bubble_sources", choices = all_idents, selected = character(0))
            updateSelectInput(session, "bubble_targets", choices = all_idents, selected = character(0))

            showNotification(paste0("CellChat done! ", length(pathways), " pathways."), type = "message")
            incProgress(1)
          }, error = function(e) {
            message("[CellChat] ERROR: ", e$message)
            showNotification(paste("CellChat Error:", e$message), type = "error", duration = 20)
          })
        })

      } else {
        # ---- Compare two conditions mode ----
        cond_col <- input$condition_col
        label_a  <- input$condition_a
        label_b  <- input$condition_b

        if (is.null(cond_col) || !nzchar(cond_col))
          { showNotification("Please select a Condition Column.", type = "error"); return() }
        if (label_a == label_b)
          { showNotification("Condition A and B must be different.", type = "error"); return() }

        withProgress(message = "Running CellChat (compare)...", value = 0, {
          tryCatch({
            .log_mem("compare mode start")
            incProgress(0.05, detail = paste("Preparing Condition A:", label_a))
            so_a <- so[, so@meta.data[[cond_col]] == label_a]
            if (ncol(so_a) == 0) stop(paste("No cells found for Condition A:", label_a))
            message("[CellChat] Condition A cells: ", ncol(so_a))

            cc_a <- .run_single_cellchat(so_a, group_by_col, input$cellchat_db, max_cells = max_cells_val)
            # Explicitly free Condition A's Seurat subset before processing B
            rm(so_a); gc(verbose = FALSE)
            .log_mem("after Condition A complete (so_a freed)")
            incProgress(0.45, detail = paste("Condition A done. Running Condition B:", label_b))

            so_b <- so[, so@meta.data[[cond_col]] == label_b]
            if (ncol(so_b) == 0) stop(paste("No cells found for Condition B:", label_b))
            message("[CellChat] Condition B cells: ", ncol(so_b))

            cc_b <- .run_single_cellchat(so_b, group_by_col, input$cellchat_db, max_cells = max_cells_val)
            rm(so_b); gc(verbose = FALSE)
            .log_mem("after Condition B complete (so_b freed)")
            incProgress(0.85, detail = "Merging conditions...")

            # liftCellChat to unify cluster space
            all_labels <- union(levels(cc_a@idents), levels(cc_b@idents))
            message("[CellChat] Lifting to common label space: ", paste(all_labels, collapse = ", "))
            cc_a <- CellChat::liftCellChat(cc_a, group.new = all_labels)
            cc_b <- CellChat::liftCellChat(cc_b, group.new = all_labels)

            merged <- CellChat::mergeCellChat(
              list(cc_a, cc_b),
              add.names = c(label_a, label_b)
            )
            message("[CellChat] mergeCellChat OK")

            cellchat_list(list(A = cc_a, B = cc_b))
            cellchat_merged(merged)
            cellchat_results(NULL)

            # Update pathway selector from merged
            pathways_a <- cc_a@netP$pathways
            pathways_b <- cc_b@netP$pathways
            shared_pathways <- union(pathways_a, pathways_b)
            updateSelectInput(session, "pathway_show",
              choices = if (length(shared_pathways) > 0) shared_pathways else "(none detected)")

            # Update bubble cluster selectors
            all_idents <- all_labels
            updateSelectInput(session, "bubble_sources", choices = all_idents, selected = character(0))
            updateSelectInput(session, "bubble_targets", choices = all_idents, selected = character(0))

            n_a <- length(pathways_a); n_b <- length(pathways_b)
            showNotification(
              paste0("CellChat comparison done! ", label_a, ": ", n_a,
                     " pathways, ", label_b, ": ", n_b, " pathways."),
              type = "message", duration = 10
            )
            incProgress(1)
          }, error = function(e) {
            message("[CellChat] ERROR: ", e$message)
            showNotification(paste("CellChat Error:", e$message), type = "error", duration = 20)
          })
        })
      }
    })

    # --- Render plot ---
    render_cellchat_plot <- function() {
      v_max    <- input$vertex_size_max %||% 10
      e_max    <- input$edge_width_max  %||% 8
      lbl_sz   <- input$label_font_size %||% 10
      base_fs  <- input$base_font_size  %||% 14
      facet_fs <- input$facet_font_size %||% 12

      # ggplotにfont sizeを適用するヘルパー
      apply_theme <- function(p) {
        p + ggplot2::theme(
          text         = ggplot2::element_text(size = base_fs),
          axis.text    = ggplot2::element_text(size = base_fs),
          legend.text  = ggplot2::element_text(size = base_fs),
          strip.text   = ggplot2::element_text(size = facet_fs)
        )
      }

      if (input$analysis_mode == "single") {
        req(cellchat_results())
        cc <- cellchat_results()

        if (input$plot_type == "Network Circle") {
          req(input$pathway_show)
          if (!input$pathway_show %in% cc@netP$pathways) {
            plot.new(); text(0.5, 0.5, "Pathway not detected.\nTry 'All' database.", cex = 1.2); return()
          }
          CellChat::netVisual_aggregate(cc, signaling = input$pathway_show, layout = "circle",
            vertex.size.max = v_max, edge.width.max = e_max, vertex.label.cex = lbl_sz / 10)

        } else {
          sources <- if (length(input$bubble_sources) > 0) input$bubble_sources else unique(levels(cc@idents))
          targets <- if (length(input$bubble_targets) > 0) input$bubble_targets else unique(levels(cc@idents))
          p <- CellChat::netVisual_bubble(cc, sources.use = sources, targets.use = targets,
                                          remove.isolate = FALSE)
          print(apply_theme(p))
        }

      } else {
        req(cellchat_merged(), cellchat_list())
        merged  <- cellchat_merged()
        cc_list <- cellchat_list()
        ptype   <- input$compare_plot_type

        # condition A/B の表示名を取得
        label_a <- input$condition_a %||% "Condition A"
        label_b <- input$condition_b %||% "Condition B"

        if (ptype == "pathway_compare") {
          p <- CellChat::rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = FALSE)
          # legendラベルをsample番号→条件名に置換（既存のfill colorsを保持）
          current_colors <- ggplot2::ggplot_build(p)$data[[1]]$fill
          unique_colors  <- unique(current_colors)
          p <- p + ggplot2::scale_fill_manual(
            values = setNames(unique_colors, c(label_a, label_b)),
            labels = c(label_a, label_b),
            name   = "Condition"
          )
          print(apply_theme(p))

        } else if (ptype == "diff_network") {
          measure <- input$diff_measure %||% "count"
          CellChat::netVisual_diffInteraction(merged, comparison = c(1, 2), measure = measure,
            vertex.size.max = v_max, edge.width.max = e_max, vertex.label.cex = lbl_sz / 10)

        } else if (ptype == "circle_per_cond") {
          req(input$pathway_show)
          cc_sel <- cc_list[[input$which_cond %||% "A"]]
          if (!input$pathway_show %in% cc_sel@netP$pathways) {
            plot.new(); text(0.5, 0.5, "Pathway not detected in this condition.", cex = 1.2); return()
          }
          CellChat::netVisual_aggregate(cc_sel, signaling = input$pathway_show, layout = "circle",
            vertex.size.max = v_max, edge.width.max = e_max, vertex.label.cex = lbl_sz / 10)

        } else if (ptype == "bubble_per_cond") {
          cc_sel  <- cc_list[[input$which_cond %||% "A"]]
          all_ids <- unique(levels(cc_sel@idents))
          sources <- if (length(input$bubble_sources) > 0) input$bubble_sources else all_ids
          targets <- if (length(input$bubble_targets) > 0) input$bubble_targets else all_ids
          p <- CellChat::netVisual_bubble(cc_sel, sources.use = sources, targets.use = targets,
                                          remove.isolate = FALSE)
          print(apply_theme(p))
        }
      }
    }

    output$networkPlot <- renderPlot({
      render_cellchat_plot()
    }, width = reactive(input$plot_width %||% 700), height = reactive(input$plot_height %||% 700))

    output$download_plot <- downloadHandler(
      filename = function() {
        mode  <- if (input$analysis_mode == "single") input$plot_type else input$compare_plot_type
        paste0("CellChat_", gsub(" ", "_", mode), "_", Sys.Date(), ".pptx")
      },
      content = function(file) {
        is_base_plot <- (input$analysis_mode == "single" && input$plot_type == "Network Circle") ||
                        (input$analysis_mode == "compare" && input$compare_plot_type %in% c("diff_network", "circle_per_cond"))
        if (is_base_plot) {
          save_baseplot_as_pptx(file, render_cellchat_plot, input$plot_width, input$plot_height)
        } else {
          p <- render_cellchat_plot()
          save_plot_as_pptx(file, p, input$plot_width, input$plot_height)
        }
      }
    )
  })
}
