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
library(scater)
library(scran)
library(vegan)
library(ggalluvial)
library(tools)
library(shinyalert)
library(openxlsx)    # Required for reading scType marker files
library(HGNChelper)  # Required for scType gene name conversion

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
source('module/4_12_clonalAbundance.R', local = TRUE)
source('module/4_13_publicClonotype.R', local = TRUE)

source('utils.R', local = TRUE)

# UI Definition
ui <- navbarPage(
  includeCSS("style.css"),
  title = "SingleRepExplorer: Web Application for Single Cell Repertoire Analysis",
  footer = shiny::HTML("<div style='text-align: center; padding: 10px;'>© 2025 SingleRepExplorer. All rights reserved.</div>"),
  
  tabPanel(
    "Upload and Run",
    uploadCellrangerUI("upload")
  ),
  tabPanel(
    "Dimensional Plot",
    reductionPlotUI("reduction_plot")
  ),
  navbarMenu(
    "Transcriptome",
    tabPanel(
      "Expression Visualization",
      geneExpressionUI("gene_expression")
    ),
    tabPanel(
      "Differential Expression",
      differentialGeneExpressionUI("differential_gene_expression")
    )
  ),
  navbarMenu(
    "Repertoire",
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
      "Track Clonotype",
      trackClonotypeUI("track_clonotype")
    )
  )
)

# Server logic definition
server <- function(input, output, session) {
  # Initialize reactive values
  myReactives <- reactiveValues(grouping_updated = 0)

  # Enable data upload server module
  uploadCellrangerServer("upload", myReactives)

  # Load filename.RData to set initial data if exists
  if (file.exists('filename.RData')) {
    load('filename.RData', envir = .GlobalEnv)
  }

  # Call server modules for each analysis tab
  reductionPlotServer("reduction_plot", myReactives)
  geneExpressionServer("gene_expression", myReactives)
  differentialGeneExpressionServer('differential_gene_expression', myReactives)
  clonotypeInformationServer("clonotype_information", myReactives)
  uniqueClonesServer('unique_clones', myReactives)
  clonalAbundancePlotServer('clonal_abundance_plot', myReactives)
  cdrLengthDistServer('cdr_length_dist', myReactives)
  diversityAnalysisServer('diversity_analysis', myReactives)
  trackClonotypeServer('track_clonotype', myReactives)
  clonalOverlapServer("clonal_overlap", myReactives)
  publicClonotypeServer("public_clonotype", myReactives)
  frequencyPlotServer('frequency_plot', myReactives)
}

# Run the application
shinyApp(ui, server)