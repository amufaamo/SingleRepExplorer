#  download .rds section was commented out.

uploadCellrangerUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Upload cellranger output files"),
      p('Upload these files and press "Run".'),
      h4("1. count file"),
      p("This is the file for clustering."),
      p(".../outs/count/filtered_feature_bc_matrix.h5"),
      h4("2. TCR file (if exists)"),
      p("This is the file for TCR. If you analyzed TCR, upload the csv file"),
      p(".../outs/vdj_t/filtered_contig_annotations.csv"),
      h4("3. BCR file (if exists)"),
      p("This is the file for BCR. If you analyzed BCR, upload the csv file"),
      p(".../outs/vdj_b/filtered_contig_annotations.csv")
    ),
    mainPanel(
      fluidRow(
        column(6,
          fileInput(ns("h5"),  "1. Choose .h5  file"),
          fileInput(ns("tcr"), "2. Choose .tcr file (optional)"),
          fileInput(ns("bcr"), "3. Choose .bcr file (optional)"),
          actionButton(ns("run"), "Run")
        ),
        )
      )
    )
}


uploadCellrangerServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$run, {
      myReactives <- fileUpload(input, myReactives)
    })
    
    observeEvent(myReactives$h5_path, {
      myReactives <- h5_to_seurat_object(myReactives)
      myReactives <- run_seurat_object(myReactives)
#      myReactives <- add_celltype(myReactives)
    })
    
    observeEvent(myReactives$tcr_path, {
      req(myReactives$seurat_object)
      myReactives <- dataframe_tcr(myReactives)
      myReactives <- addMetadata(myReactives, 'tcr')
      myReactives <- addClonotypeIdToTCR(myReactives)
      myReactives <- addMetadataToTCR(myReactives)
      myReactives$seurat_object@meta.data <- dplyr::left_join(myReactives$seurat_object@meta.data, myReactives$metadata, by='sample')
      addTCRClonotypeIdToSeuratObject(myReactives)
      rownames(myReactives$seurat_object@meta.data) <- myReactives$seurat_object@meta.data$barcode
      myReactives$seurat_object <- combineExpression(myReactives$tcr_df, myReactives$seurat_object)
      current_colnames <- names(myReactives$seurat_object@meta.data)
      current_colnames <- ifelse(current_colnames %in% c("CTgene", "CTnt", "CTaa", "CTstrict", "clonalProportion", "clonalFrequency","cloneSize"), paste0("TCR_", current_colnames),current_colnames)
      names(myReactives$seurat_object@meta.data) <- current_colnames
      myReactives$seurat_object@meta.data <- myReactives$seurat_object@meta.data %>% mutate(TCR = case_when(
        !is.na(TCR_clonalFrequency) ~ "1.TRUE",
        TRUE ~ "2. FALSE"
      ))
      saveRDS(myReactives$tcr_df,'tcr_df.rds')
      
    })
    
    
    observeEvent(myReactives$bcr_path, {
      req(myReactives$seurat_object)
      myReactives <- dataframe_bcr(myReactives)
      myReactives$bcr_df$S1$sample <- NULL
      if (is.null(myReactives$tcr_path)) {
        myReactives <- addMetadata(myReactives, 'bcr')
        myReactives$seurat_object@meta.data <- dplyr::left_join(myReactives$seurat_object@meta.data, myReactives$metadata, by='sample')
        rownames(myReactives$seurat_object@meta.data) <- myReactives$seurat_object@meta.data$barcode
      }
      myReactives <- addClonotypeIdToBCR(myReactives)
      myReactives <- addMetadataToBCR(myReactives)
      myReactives$seurat_object <- combineExpression(myReactives$bcr_df, myReactives$seurat_object)
      current_colnames <- names(myReactives$seurat_object@meta.data)
      current_colnames <- ifelse(current_colnames %in% c("CTgene", "CTnt", "CTaa", "CTstrict", "clonalProportion", "clonalFrequency","cloneSize"), paste0("BCR_", current_colnames),current_colnames)
      names(myReactives$seurat_object@meta.data) <- current_colnames
      myReactives$seurat_object@meta.data <- myReactives$seurat_object@meta.data %>% mutate(BCR = case_when(
        !is.na(BCR_clonalFrequency) ~ "1.TRUE",
        TRUE ~ "2. FALSE"
      ))
      saveRDS(myReactives$bcr_df,'bcr_df.rds')
      saveRDS(myReactives$seurat_object, 'seurat_object.rds')
    })
    
  })
}

fileUpload <- function(input, myReactives) {
  myReactives$h5_path <- input$h5$datapath
  myReactives$tcr_path <- input$tcr$datapath
  myReactives$bcr_path <- input$bcr$datapath
  return(myReactives)
}

h5_to_seurat_object <- function(myReactives) {
  h5 <- Read10X_h5(myReactives$h5_path)
  if (!is.null(names(h5)) & ("Gene Expression" %in% names(h5))) { # names(h5) are NULL OR c('Gene Expression','Antibody Capture')
    h5 <- h5[["Gene Expression"]]
  }
  
  myReactives$seurat_object <- CreateSeuratObject(h5)
  return(myReactives)
}

run_seurat_object <- function(myReactives) {
  seurat_object <- myReactives$seurat_object
  seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    mutate(sample = str_remove(barcode, "^.+-"))
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), npcs = 50)
  seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  seurat_object <- RunUMAP(seurat_object, dims = 1:10)
  seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  myReactives$seurat_object <- seurat_object
  myReactives$meta.data <- seurat_object@meta.data
  return(myReactives)
}

celltyping <- function(myReactives){
  
}

dataframe_tcr <- function(myReactives){
  df <- read.csv(myReactives$tcr_path)
  myReactives$tcr_df <- combineTCR(df, removeNA = FALSE, filterMulti = TRUE)
  return(myReactives)
}

dataframe_bcr <- function(myReactives){
  df <- read.csv(myReactives$bcr_path)
  myReactives$bcr_df <- combineBCR(df, samples = 'S1', removeNA = FALSE, filterMulti = TRUE)
  return(myReactives)
}

addClonotypeIdToTCR <- function(myReactives){
  S1 <- myReactives$tcr_df$S1
  tcr <- read.csv(myReactives$tcr_path)
  tcr <- tcr %>% select(barcode, raw_clonotype_id, exact_subclonotype_id) %>% distinct()
  S1 <- dplyr::left_join(S1, tcr, by = 'barcode')
  myReactives$tcr_df$S1 <- S1
  return(myReactives)
}

addTCRClonotypeIdToSeuratObject <- function(myReactives){
  tcr <- read.csv(myReactives$tcr_path)
  tcr <- tcr %>% select(barcode, raw_clonotype_id, exact_subclonotype_id) %>% distinct()
  names(tcr) <- c("barcode", "TCR_raw_clonotype_id", "TCR_exact_subclonotype_id")
  myReactives$seurat_object@meta.data <- dplyr::left_join(myReactives$seurat_object@meta.data, tcr, by = 'barcode')
  rownames(myReactives$seurat_object@meta.data) <- myReactives$seurat_object@meta.data$barcode
  return(myReactives)
}

addClonotypeIdToBCR <- function(myReactives){
  S1 <- myReactives$bcr_df$S1
  bcr <- read.csv(myReactives$bcr_path)
  bcr <- bcr %>% select(barcode, raw_clonotype_id, exact_subclonotype_id) %>% distinct()
  S1 <- dplyr::left_join(S1, bcr, by = 'barcode')
  myReactives$bcr_df$S1 <- S1
  return(myReactives)
}

addMetadataToTCR <- function(myReactives){
  metadata <- myReactives$seurat_object@meta.data %>% select(sample, seurat_clusters, barcode)
  S1 <- myReactives$tcr_df$S1
  S1 <- dplyr::filter(S1, barcode %in% myReactives$seurat_object@meta.data$barcode) 
  S1 <- dplyr::left_join(S1, metadata, by = 'barcode')
  S1 <- dplyr::left_join(S1, myReactives$metadata, by = 'sample')
  myReactives$tcr_df$S1 <- S1
  return(myReactives)
}

addMetadataToBCR <- function(myReactives){
  metadata <- myReactives$seurat_object@meta.data %>% select(sample, seurat_clusters, barcode)
  S1 <- myReactives$bcr_df$S1
  S1$barcode <- gsub("^S1_", "", S1$barcode)
  print(S1)
  S1 <- dplyr::filter(S1, barcode %in% myReactives$seurat_object@meta.data$barcode) 
  S1 <- dplyr::left_join(S1, metadata, by = 'barcode')
  S1 <- dplyr::left_join(S1, myReactives$metadata, by = 'sample')
  myReactives$bcr_df$S1 <- S1
  return(myReactives)
}

addMetadata <- function(myReactives, vdj = 'tcr'){
  delete_column <- c('barcode', 
                     'is_cell', 
                     'contig_id', 
                     'high_confidence', 
                     'length', 
                     'chain',
                     'v_gene', 
                     'd_gene', 
                     'j_gene', 
                     'c_gene', 
                     'full_length', 
                     'productive', 
                     'fwr1', 'fwr1_nt', 
                     'cdr1', 'cdr1_nt', 
                     'fwr2', 'fwr2_nt', 
                     'cdr2', 'cdr2_nt', 
                     'fwr3', 'fwr3_nt', 
                     'cdr3', 'cdr3_nt', 
                     'fwr4', 'fwr4_nt', 
                     'reads', 
                     'umis', 
                     'raw_clonotype_id', 
                     'raw_consensus_id', 
                     'exact_subclonotype_id')
  
  if(vdj == 'tcr'){
    df <- read.csv(myReactives$tcr_path)
  } else if (vdj == 'bcr'){
    df <- read.csv(myReactives$bcr_path)
  }
  df <- df %>% mutate(sample = str_remove(barcode, "^.+-"))
  df <- df %>% select(-delete_column, -starts_with("TCR_"), -starts_with("BCR_"),) %>% distinct()
  print(df)
  myReactives$metadata <- df
  return(myReactives)
}

