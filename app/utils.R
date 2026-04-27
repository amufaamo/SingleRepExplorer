# SingleRepExplorer - Utility Functions
# Contains helper functions for UI updates and data processing

# --- Helper Functions for UI Updates ---

# Updates a SelectInput with categorical metadata columns from Seurat object
update_group_by_select_input <- function(session, myReactives, selected = NULL) {
  req(myReactives$seurat_object)
  
  meta_data <- myReactives$seurat_object@meta.data

  # 1. Filter for categorical types (character, factor, logical)
  is_categorical <- sapply(names(meta_data), function(col) {
    is.character(meta_data[[col]]) || is.factor(meta_data[[col]]) || is.logical(meta_data[[col]])
  })
  potential_choices <- names(meta_data)[is_categorical]

  # 2. Exclude specific columns that are not suitable for grouping
  minus_column <- c("orig.ident", "barcode") 
  choices <- setdiff(potential_choices, minus_column)
  
  # Exclude technical or sequence-heavy columns
  choices <- choices[!grepl(
    "^(RNA_snn_res\\.|TCR_CT|BCR_CT|TCR_clonotype_id|BCR_clonotype_id|TCR_TRA_|TCR_TRB_|BCR_IGH_|BCR_IGL_|BCR_IGK_)",
    choices
  )]

  shiny::validate(shiny::need(length(choices) > 0, "No suitable metadata columns found for grouping."))

  # 3. Determine default selection
  selected_value <- if (!is.null(selected) && selected %in% choices) {
    selected
  } else if ("sample" %in% choices) {
    "sample"
  } else {
    choices[1]
  }
  
  # 4. Update the UI
  updateSelectInput(session, "group_by", choices = choices, selected = selected_value)
}

# Updates a SelectInput with available dimensionality reductions
update_reduction_choices <- function(session, myReactives) {
  req(myReactives$seurat_object)
  reduction_names   <- names(myReactives$seurat_object@reductions)
  choices           <- stats::setNames(reduction_names, toupper(reduction_names))
  default_selection <- if ("umap" %in% reduction_names) "umap" else reduction_names[1]
  updateSelectInput(session, "reduction", choices = choices, selected = default_selection)
}

# Updates CheckboxGroupInput with unique groups for sub-selection
update_unique_group_choices <- function(session, myReactives, group_by_col) {
  req(myReactives$seurat_object, group_by_col, nzchar(group_by_col))
  shiny::validate(shiny::need(group_by_col %in% names(myReactives$seurat_object@meta.data),
                paste("Column '", group_by_col %||% "", "' not found in Seurat metadata.")))

  unique_groups <- tryCatch({
      unique(myReactives$seurat_object@meta.data[[group_by_col]])
  }, error = function(e) {
      warning("Error getting unique groups for '", group_by_col, "': ", e$message); NULL
  })
  req(unique_groups)

  # Remove NAs
  unique_groups <- na.omit(unique_groups)
  shiny::validate(shiny::need(length(unique_groups) > 0, paste("No valid unique groups found in column '", group_by_col %||% "", "'.")))

  # Sorting logic (numeric vs character)
  if (is.factor(unique_groups)) { unique_groups <- levels(unique_groups) }
  if (all(!is.na(suppressWarnings(as.numeric(as.character(unique_groups)))))) {
     numeric_representation <- suppressWarnings(as.numeric(as.character(unique_groups)))
     if(all(!is.na(numeric_representation))){
         unique_groups <- as.character(sort(numeric_representation))
     } else {
         unique_groups <- sort(as.character(unique_groups))
     }
  } else { unique_groups <- sort(as.character(unique_groups)) }

  # Select first by default
  selected_values <- if (length(unique_groups) > 0) unique_groups[1] else NULL

  updateCheckboxGroupInput(session, "unique_group", choices = unique_groups, selected = selected_values, inline = TRUE)
}

# Updates Group By and Cluster choices specifically for Differential Expression/Marker identification
update_group_by_for_marker <- function(session, input, myReactives) {
  req(myReactives$seurat_object)
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt")

  metadatas <- tryCatch({
      myReactives$seurat_object@meta.data %>%
        dplyr::select(-any_of(minus_column), -starts_with("RNA_snn_res.")) %>%
        dplyr::select(-starts_with("TCR_CT"), -starts_with("BCR_CT")) 
    }, error = function(e){ warning("Error selecting metadata in update_group_by_for_marker: ", e$message); NULL })
  req(metadatas)

  metadata_cols <- names(metadatas)
  shiny::validate(shiny::need(length(metadata_cols) > 0, "No suitable metadata columns found."))

  group_cols <- setNames(as.list(metadata_cols), metadata_cols)
  selected_group <- if ("seurat_clusters" %in% metadata_cols) "seurat_clusters" else metadata_cols[1]

  updateSelectInput(session, "group_by", choices = group_cols, selected = selected_group)

  # Update dependent cluster selections
  current_group_by <- input$group_by %||% selected_group 

  if (!is.null(current_group_by) && current_group_by %in% colnames(myReactives$seurat_object@meta.data)) {
      cluster_choices <- tryCatch(unique(myReactives$seurat_object@meta.data[[current_group_by]]), error = function(e) NULL)
      if (!is.null(cluster_choices)) {
          cluster_choices <- sort(na.omit(cluster_choices))
           updateSelectInput(session, "target_cluster", choices = cluster_choices, selected = cluster_choices[1])
           updateSelectInput(session, "reference_cluster", choices = cluster_choices, selected = NULL)
      }
  }
}

# --- UI Component Helper Functions ---

# Legend position selector
legendPositionInput <- function(ns, selected = "right") {
  selectInput(ns("legend"), "Legend Position", choices = c("right", "left", "bottom", "top", "none"), selected = selected)
}

# Plot width input
plotWidthInput <- function(ns, value = 500, min = 100, max = 2000, step = 100) {
  numericInput(ns("plot_width"), "Plot Width (px)", min = min, max = max, value = value, step = step)
}

# Plot height input
plotHeightInput <- function(ns, value = 500, min = 100, max = 2000, step = 100) {
  numericInput(ns("plot_height"), "Plot Height (px)", min = min, max = max, value = value, step = step)
}

# Grouped UI for common plot options (Legend, Width, Height)
commonPlotOptions <- function(ns, legend_selected = "right", width_value = 500, height_value = 500) {
  tagList(
    legendPositionInput(ns, selected = legend_selected),
    div(style = "display: flex; gap: 10px;",
      plotWidthInput(ns, value = width_value),
      plotHeightInput(ns, value = height_value)
    )
  )
}

# Generic reduction input placeholder
reductionInput <- function(ns) {
  selectInput(ns("reduction"), "Reduction Method", choices = list("Loading..." = ""), selected = "")
}

# Generic grouping input
groupByInput <- function(ns, choices = c("sample", "seurat_clusters"), selected = "sample") {
  selectInput(ns("group_by"), "Group by", choices = choices, selected = selected)
}

# Point size selection
pointSizeInput <- function(ns, value = 0.1, min = 0.01, max = 10, step = 0.01) {
  numericInput(ns("point_size"), "Point Size", min = min, max = max, value = value, step = step)
}

# Label size selection
labelSizeInput <- function(ns, value = 10, min = 0, max = 20, step = 1) {
  numericInput(ns("label_size"), "Label Size", min = min, max = max, value = value, step = step)
}

# Base font size input
baseFontSizeInput <- function(ns, value = 12, min = 6, max = 24, step = 1) {
  numericInput(ns("base_font_size"), "Base Font Size", min = min, max = max, value = value, step = step)
}

# X-axis label angle selector
xAxisAngleInput <- function(ns, selected = "45") {
  selectInput(ns("x_axis_angle"), "X-axis Label Angle",
              choices = c("0°" = "0", "45°" = "45", "90°" = "90"), selected = selected)
}

# Jitter points toggle + size
jitterPointsInput <- function(ns) {
  tagList(
    checkboxInput(ns("show_jitter"), "Show Individual Points (Jitter)", value = FALSE),
    conditionalPanel(
      condition = sprintf("input['%s'] == true", ns("show_jitter")),
      numericInput(ns("jitter_size"), "Jitter Point Size", min = 0.1, max = 5, value = 0.5, step = 0.1)
    )
  )
}

# VDJ type selection (TCR/BCR)
vdjType <- function(ns){
  selectInput(ns('vdj_type'), label = 'VDJ Type', choices = c("TCR" = "tcr", "BCR" = "bcr"), selected = 'tcr')
}

# --- Data Processing and Integration ---

# Loads H5 file and converts to Seurat object
h5_to_seurat_object <- function(myReactives) {
  req(myReactives$h5_path)
  
  # Use pre-loaded raw Seurat if available (from QC step)
  if (!is.null(myReactives$raw_qc_seurat)) {
    message("Using pre-loaded raw Seurat object from QC preview.")
    myReactives$seurat_object <- myReactives$raw_qc_seurat
    
    # Keeping raw_qc_seurat in memory allows the user to adjust filtering and re-run!
    # myReactives$raw_qc_seurat <- NULL
    # gc()
    
    return(myReactives)
  }
  
  h5 <- Seurat::Read10X_h5(myReactives$h5_path)
  if (is.list(h5) && !is.data.frame(h5) && ("Gene Expression" %in% names(h5))) {
    h5 <- h5[["Gene Expression"]]
  }
  raw_so <- Seurat::CreateSeuratObject(h5)
  raw_so[["percent.mt"]] <- Seurat::PercentageFeatureSet(raw_so, pattern = "^MT-")
  myReactives$seurat_object <- raw_so
  return(myReactives)
}

# Processes TCR CSV into scRepertoire list (scRepertoire required)
dataframe_tcr <- function(myReactives){
  req(myReactives$tcr_path)
  df <- tryCatch(read.csv(myReactives$tcr_path), error = function(e) {warning("Error reading TCR CSV: ", e$message); NULL})
  req(df)
  
  if (!requireNamespace("scRepertoire", quietly = TRUE)) {
     warning("scRepertoire package is needed for dataframe_tcr function.")
     myReactives$tcr_df <- df
     return(myReactives)
  }
  
  df_list <- df %>%
      mutate(sample_id = str_remove(barcode, "^.+-")) %>%
      group_split(sample_id) %>%
      setNames(unique(.$sample_id))

  combined_list <- tryCatch(scRepertoire::combineTCR(df_list, removeNA = TRUE, filterMulti = TRUE),
                           error = function(e) { warning("combineTCR failed: ", e$message); NULL })
  req(combined_list)
  myReactives$tcr_list <- combined_list
  return(myReactives)
}

# Processes BCR CSV into scRepertoire list
dataframe_bcr <- function(myReactives){
  req(myReactives$bcr_path)
  df <- tryCatch(read.csv(myReactives$bcr_path), error = function(e) {warning("Error reading BCR CSV: ", e$message); NULL})
  req(df)
  if (!requireNamespace("scRepertoire", quietly = TRUE)) {
     warning("scRepertoire package is needed for dataframe_bcr function.")
     myReactives$bcr_df <- df
     return(myReactives)
  }
  
  df_list <- df %>%
      mutate(sample_id = str_remove(barcode, "^.+-")) %>%
      group_split(sample_id) %>%
      setNames(unique(.$sample_id))

  combined_list <- tryCatch(scRepertoire::combineBCR(df_list, removeNA = TRUE, filterMulti = TRUE),
                           error = function(e) { warning("combineBCR failed: ", e$message); NULL })
  req(combined_list)
  myReactives$bcr_list <- combined_list
  return(myReactives)
}

# Join TCR clonotype IDs to Seurat metadata
addTCRClonotypeIdToSeuratObject <- function(myReactives){
  req(myReactives$tcr_path, myReactives$seurat_object)
  tcr <- tryCatch(read.csv(myReactives$tcr_path), error=function(e) {warning("Failed to read TCR file in addTCRClonotypeId: ", e$message); NULL})
  req(tcr)
  
  tcr_select <- tcr %>%
    dplyr::select(any_of(c("barcode", "raw_clonotype_id", "exact_subclonotype_id"))) %>%
    distinct()
  
  req("barcode" %in% names(tcr_select), "raw_clonotype_id" %in% names(tcr_select))

  tcr_renamed <- tcr_select %>%
     rename(TCR_raw_clonotype_id = raw_clonotype_id)
  if ("exact_subclonotype_id" %in% names(tcr_renamed)) {
      tcr_renamed <- tcr_renamed %>% rename(TCR_exact_subclonotype_id = exact_subclonotype_id)
  }

  seurat_meta <- myReactives$seurat_object@meta.data %>%
    tibble::rownames_to_column("barcode_rownames")

  join_col <- if ("barcode" %in% names(seurat_meta)) "barcode" else "barcode_rownames"
  if (!"barcode" %in% names(tcr_renamed)) {
      warning("TCR data does not have 'barcode' column for joining.")
      return(myReactives)
  }

  tcr_renamed$barcode <- as.character(tcr_renamed$barcode)
  seurat_meta[[join_col]] <- as.character(seurat_meta[[join_col]])

  joined_meta <- dplyr::left_join(seurat_meta, tcr_renamed, by = setNames("barcode", join_col))

  rownames(joined_meta) <- joined_meta$barcode_rownames
  joined_meta$barcode_rownames <- NULL

  myReactives$seurat_object@meta.data <- joined_meta
  return(myReactives)
}

# Join BCR clonotype IDs to Seurat metadata
addBCRClonotypeIdToSeuratObject <- function(myReactives){
  req(myReactives$bcr_path, myReactives$seurat_object)
  bcr <- tryCatch(read.csv(myReactives$bcr_path), error=function(e) {warning("Failed to read BCR file in addBCRClonotypeId: ", e$message); NULL})
  req(bcr)
  
  bcr_select <- bcr %>%
    dplyr::select(any_of(c("barcode", "raw_clonotype_id", "exact_subclonotype_id"))) %>%
    distinct()
  
  req("barcode" %in% names(bcr_select), "raw_clonotype_id" %in% names(bcr_select))

  bcr_renamed <- bcr_select %>%
     rename(BCR_raw_clonotype_id = raw_clonotype_id)
  if ("exact_subclonotype_id" %in% names(bcr_renamed)) {
      bcr_renamed <- bcr_renamed %>% rename(BCR_exact_subclonotype_id = exact_subclonotype_id)
  }

  seurat_meta <- myReactives$seurat_object@meta.data %>%
    tibble::rownames_to_column("barcode_rownames")

  join_col <- if ("barcode" %in% names(seurat_meta)) "barcode" else "barcode_rownames"
  if (!"barcode" %in% names(bcr_renamed)) {
      warning("BCR data does not have 'barcode' column for joining.")
      return(myReactives)
  }

  bcr_renamed$barcode <- as.character(bcr_renamed$barcode)
  seurat_meta[[join_col]] <- as.character(seurat_meta[[join_col]])

  joined_meta <- dplyr::left_join(seurat_meta, bcr_renamed, by = setNames("barcode", join_col))

  rownames(joined_meta) <- joined_meta$barcode_rownames
  joined_meta$barcode_rownames <- NULL

  myReactives$seurat_object@meta.data <- joined_meta
  return(myReactives)
}

# Extracts Seurat default colors for a given metadata column
get_seurat_colors <- function(seurat_obj, group_by_col) {
  if (is.null(seurat_obj) || is.null(group_by_col) || !group_by_col %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  meta_col <- seurat_obj@meta.data[[group_by_col]]
  
  if (is.factor(meta_col)) {
    all_groups <- levels(meta_col)
  } else {
    all_groups <- sort(unique(meta_col))
    if (all(!is.na(suppressWarnings(as.numeric(all_groups))))) {
      all_groups <- as.character(sort(as.numeric(all_groups)))
    }
  }
  
  colors <- scales::hue_pal()(length(all_groups))
  names(colors) <- as.character(all_groups)
  return(colors)
}

# Enrich TCR data with Seurat metadata and external attributes
addMetadataToTCR <- function(myReactives){
  req(myReactives$seurat_object, myReactives$tcr_df, myReactives$metadata)
  
  metadata_seurat <- myReactives$seurat_object@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::select(any_of(c("sample", "seurat_clusters", "barcode")))

  validate(need(is.data.frame(myReactives$tcr_df), "myReactives$tcr_df is not a data frame."))
  validate(need("barcode" %in% names(myReactives$tcr_df), "tcr_df requires 'barcode' column."))

  S1 <- myReactives$tcr_df %>%
    dplyr::filter(barcode %in% metadata_seurat$barcode) %>%
    dplyr::left_join(metadata_seurat, by = 'barcode') %>%
    dplyr::left_join(myReactives$metadata, by = 'sample')

  myReactives$tcr_df <- S1
  return(myReactives)
}

# Enrich BCR data with Seurat metadata and external attributes
addMetadataToBCR <- function(myReactives){
  req(myReactives$seurat_object, myReactives$bcr_df, myReactives$metadata)
  
  metadata_seurat <- myReactives$seurat_object@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::select(any_of(c("sample", "seurat_clusters", "barcode")))

  validate(need(is.data.frame(myReactives$bcr_df), "myReactives$bcr_df is not a data frame."))
  validate(need("barcode" %in% names(myReactives$bcr_df), "bcr_df requires 'barcode' column."))

  S1 <- myReactives$bcr_df %>%
    dplyr::filter(barcode %in% metadata_seurat$barcode) %>%
    dplyr::left_join(metadata_seurat, by = 'barcode') %>%
    dplyr::left_join(myReactives$metadata, by = 'sample')

  myReactives$bcr_df <- S1
  return(myReactives)
}

# Extracts sample metadata from VDJ barcode format
addMetadata <- function(myReactives, vdj = 'tcr'){
  req(myReactives, vdj)
  file_path <- if(vdj == 'tcr') myReactives$tcr_path else myReactives$bcr_path
  req(file_path)

  df <- tryCatch(read.csv(file_path), error = function(e){ warning("Error reading VDJ file in addMetadata: ", e$message); NULL })
  req(df)
  validate(need("barcode" %in% names(df), "VDJ file requires 'barcode' column."))

  df_meta_extracted <- df %>%
    mutate(sample = str_remove(barcode, "^.+-")) %>%
    dplyr::select(any_of(c("barcode", "sample"))) %>%
    distinct()

  validate(need("sample" %in% names(df_meta_extracted), "Could not extract 'sample' column."))

  myReactives$metadata <- df_meta_extracted
  return(myReactives)
}

# Value type selection (Count vs Percentage)
valueType <- function(ns){
  radioButtons(ns("value_type"),
    label = "Value Type:",
    choices = c("Count" = "count", "Percentage (%)" = "percentage"),
    selected = "count")
}

# --- PowerPoint (Editable) Download Helpers ---
# Uses officer + rvg to create editable/ungroupable PowerPoint files.

#' Save a ggplot object to an editable PowerPoint file
#' @param file Output file path (.pptx)
#' @param plot_obj A ggplot object
#' @param width_px Plot width in pixels (will be converted to inches)
#' @param height_px Plot height in pixels (will be converted to inches)
save_plot_as_pptx <- function(file, plot_obj, width_px = 700, height_px = 500) {
  width_in  <- width_px / 72
  height_in <- height_px / 72
  pptx <- officer::read_pptx()
  pptx <- officer::add_slide(pptx, layout = "Blank")
  
  # Use code = print(plot_obj) instead of ggobj = plot_obj to safely handle patchwork/complex objects
  pptx <- officer::ph_with(
    pptx,
    rvg::dml(code = print(plot_obj)),
    location = officer::ph_location(
      left = 0.5, top = 0.5,
      width = width_in, height = height_in
    )
  )
  print(pptx, target = file)
}

save_baseplot_as_pptx <- function(file, plot_fn, width_px = 700, height_px = 700) {
  width_in  <- width_px / 72
  height_in <- height_px / 72
  pptx <- officer::read_pptx()
  pptx <- officer::add_slide(pptx, layout = "Blank")
  pptx <- officer::ph_with(
    pptx,
    rvg::dml(code = plot_fn()),
    location = officer::ph_location(
      left = 0.5, top = 0.5,
      width = width_in, height = height_in
    )
  )
  print(pptx, target = file)
}

#' Save a ggplot object as a high-res image embedded in an editable PowerPoint file
#' Used as a fallback for massive plots (FeaturePlot, DimPlot) where rvg would OOM
#' @param file Output file path (.pptx)
#' @param plot_obj A ggplot object
#' @param width_px Plot width in pixels
#' @param height_px Plot height in pixels
save_plot_as_pptx_image <- function(file, plot_obj, width_px = 700, height_px = 500) {
  width_in  <- width_px / 72
  height_in <- height_px / 72
  
  # Create temporary PNG file
  tmp_img <- tempfile(fileext = ".png")
  suppressWarnings(
    ggplot2::ggsave(tmp_img, plot = plot_obj, width = width_in, height = height_in, dpi = 300, bg = "white")
  )
  
  pptx <- officer::read_pptx()
  pptx <- officer::add_slide(pptx, layout = "Blank")
  pptx <- officer::ph_with(
    pptx,
    value = officer::external_img(src = tmp_img, width = width_in, height = height_in),
    location = officer::ph_location(
      left = 0.5, top = 0.5,
      width = width_in, height = height_in
    )
  )
  print(pptx, target = file)
}