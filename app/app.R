# SingleRepExplorer - Main Application
# Web Application for Single Cell Repertoire Analysis

# Load required libraries
library(ggsci)
library(RColorBrewer)
library(tidyverse)
library(Seurat)
library(shiny)
library(scRepertoire)
library(DT)
library(immunarch)
library(ggpointdensity)
library(viridis)
library(ggExtra)
library(ggpubr)
library(shinyWidgets)
library(patchwork)
library(plotly)
library(shinyjs)
library(bslib)
library(scater)
library(scran)

library(vegan)
library(ggalluvial)
library(tools)
library(shinyalert)
library(openxlsx)    # Required for reading scType marker files
library(HGNChelper)  # Required for scType gene name conversion
library(enrichR)     # Required for Pathway Enrichment Analysis

# Set maximum upload size
options(shiny.maxRequestSize = 50 * 1024^2 * 1000)

# Load function files
function_file_lst <- list.files("functions", pattern = ".R$", full.names = TRUE)
for (i in function_file_lst) {
  source(i, local = TRUE)
}

# Load module files
source('module/1_uploadCellranger.R', local = TRUE)
source('module/2_reductionPlot.R', local = TRUE)
source('module/3_1_geneExpression.R', local = TRUE)
source('module/3_2_differentialGeneExpression.R', local = TRUE)
source('module/4_1_clonotypeInformation.R', local = TRUE)
source('module/4_3_uniqueClones.R', local = TRUE)
source('module/4_4_frequencyPlot.R', local = TRUE)
source('module/4_5_clonalLength.R', local = TRUE)
source('module/4_6_diversity.R', local = TRUE)
source('module/4_7_trackClonotype.R', local = TRUE)
source('module/4_8_clonalOverlap.R', local = TRUE)
source('module/4_11_phylogeneticTree.R', local = TRUE)
source('module/4_14_shmIsotype.R', local = TRUE)
source('module/4_12_clonalAbundance.R', local = TRUE)

source('module/4_13_publicClonotype.R', local = TRUE)
source('module/5_trajectoryInference.R', local = TRUE)
source('module/3_3_phenotypeRepertoire.R', local = TRUE)
source('module/6_cellCommunication.R', local = TRUE)
source('module/subset.R', local = TRUE)
source('module/4_15_sequenceSimilarity.R', local = TRUE)
source('module/7_llmAnnotation.R', local = TRUE)
source('module/3_4_congaIntegration.R', local = TRUE)
source('module/4_16_cdr3Properties.R', local = TRUE)






source('utils.R', local = TRUE)

# Custom theme definition
app_theme <- bs_theme(
  version = 5,
  bootswatch = "litera",
  primary = "#2B5B84",
  secondary = "#8F9CA4",
  success = "#208b5e",
  base_font = font_google("Inter"),
  heading_font = font_google("Outfit")
) %>%
  bs_add_rules(
    ".navbar { box-shadow: 0 4px 12px rgba(0,0,0,0.05); }
    .nav-link { font-weight: 500 !important; }
    .card { border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.03); border: none; margin-bottom: 20px; }
    .btn { border-radius: 8px; font-weight: 500; }
    .sidebar { background-color: #f8fbff; border-right: 1px solid #e1e8ed; padding: 20px; }"
  )

# UI Definition
ui <- page_navbar(
  theme = app_theme,
  title = "SingleRepExplorer v2.1",
  window_title = "SingleRepExplorer v2.1",
  id = "main_navbar",
  header = tags$head(
    tags$meta(name = "description", content = "SingleRepExplorer v2.1: A premium web application for integrated single-cell RNA-seq and TCR/BCR repertoire analysis.")
  ),
  
  navbarMenu(
    "Data Input",
    icon = icon("upload"),
    tabPanel(
      "Upload and Run",
      uploadCellrangerUI("upload")
    )
  ),
  tabPanel(
    "Dimensional Plot",
    icon = icon("project-diagram"),
    reductionPlotUI("reduction_plot")
  ),
  navbarMenu(
    "Transcriptome",
    icon = icon("dna"),
    tabPanel(

      "Expression Visualization",
      geneExpressionUI("gene_expression")
    ),
    tabPanel(
      "Differential Expression",
      differentialGeneExpressionUI("differential_gene_expression")
    ),
    tabPanel(
      "Trajectory & Pseudotime",
      trajectoryInferenceUI("trajectory_inference")
    ),
    tabPanel(
      "Phenotype-Repertoire Correlation",
      phenotypeRepertoireUI("phenotype_repertoire")
    ),
    tabPanel(
      "Cell-Cell Communication",
      cellCommunicationUI("cell_communication")
    ),
    tabPanel(
      "Interactive Subsetting",
      subsetUI("subset")
    ),
    tabPanel(
      "LLM Auto-Annotation",
      llmAnnotationUI("llm_annotation")
    ),
    tabPanel(
      "Multimodal Integration (CoNGA)",
      congaIntegrationUI("conga_integration")
    )
  ),

  navbarMenu(
    "Repertoire",
    icon = icon("viruses"),
    # Basic Properties group
    "Basic Properties",

    tabPanel(
      "Clonotype Information",
      clonotypeInformationUI("clonotype_information")
    ),
    tabPanel(
      "Frequency Analysis",
      frequencyPlotUI("frequency_plot")
    ),
    tabPanel(
      "Clonal Abundance",
      clonalAbundancePlotUI("clonal_abundance_plot")
    ),
    tabPanel(
      "Diversity",
      diversityAnalysisUI("diversity_analysis")
    ),
    tabPanel(
      "Clonal Length",
      cdrLengthDistUI("cdr_length_dist")
    ),
    tabPanel(
      "CDR3 Physicochemical Profiling",
      cdr3PropertiesUI("cdr3_properties")
    ),
    # Separator line between groups

    "----",
    # Clonotype comparison group
    "Clonotype comparison",
    tabPanel(
      "Clonal Overlap",
      clonalOverlapUI("clonal_overlap")
    ),
    tabPanel(
      "Public Clonotype",
      publicClonotypeUI("public_clonotype")
    ),
    tabPanel(
      "Sequence Similarity",
      sequenceSimilarityUI("sequence_similarity")
    ),
    tabPanel(
      "Track Clonotype",
      trackClonotypeUI("track_clonotype")
    ),
    # Evolutionary analysis group
    "----",
    "Evolutionary Analysis",
    tabPanel(
      "BCR Lineage Tree",
      phylogeneticTreeUI("phylogenetic_tree")
    ),
    tabPanel(
      "BCR SHM & Isotype",
      shmIsotypeUI("shm_isotype")
    )
  ),
  footer = shiny::HTML("<div style='text-align: center; padding: 20px; background: #f8fbff; border-top: 1px solid #e1e8ed; color: #666; font-size: 0.9em;'>© 2026 SingleRepExplorer v2.1 - Developed for Premium scRNA-seq & VDJ Analysis.</div>")
)

# Server logic definition

server <- function(input, output, session) {
  # Initialize reactive values
  myReactives <- reactiveValues(grouping_updated = 0)

  # Enable data upload server module
  uploadCellrangerServer("upload", myReactives)

  # Load filename.RData to set initial data if exists (Disabled: causes cache/clutter issues)
  # if (file.exists('filename.RData')) {
  #   load('filename.RData', envir = .GlobalEnv)
  # }

  # Call server modules for each analysis tab
  reductionPlotServer("reduction_plot", myReactives)
  geneExpressionServer("gene_expression", myReactives)
  differentialGeneExpressionServer('differential_gene_expression', myReactives)
  trajectoryInferenceServer("trajectory_inference", myReactives)
  phenotypeRepertoireServer("phenotype_repertoire", myReactives)
  cellCommunicationServer("cell_communication", myReactives)
  clonotypeInformationServer("clonotype_information", myReactives)
  uniqueClonesServer('unique_clones', myReactives)
  clonalAbundancePlotServer('clonal_abundance_plot', myReactives)
  cdrLengthDistServer('cdr_length_dist', myReactives)
  diversityAnalysisServer('diversity_analysis', myReactives)
  trackClonotypeServer('track_clonotype', myReactives)
  clonalOverlapServer("clonal_overlap", myReactives)
  publicClonotypeServer("public_clonotype", myReactives)
  sequenceSimilarityServer("sequence_similarity", myReactives)
  frequencyPlotServer('frequency_plot', myReactives)
  phylogeneticTreeServer("phylogenetic_tree", myReactives)

  shmIsotypeServer("shm_isotype", myReactives)
  subsetServer("subset", myReactives)
  llmAnnotationServer("llm_annotation", myReactives)
  congaIntegrationServer("conga_integration", myReactives)
  cdr3PropertiesServer("cdr3_properties", myReactives)
}






# Run the application
shinyApp(ui, server)