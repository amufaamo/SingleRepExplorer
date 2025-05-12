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
library(shinyjs) # ステータス表示更新のため
library(scater) # combineExpressionなどで必要になる可能性
library(scran)  # combineExpressionなどで必要になる可能性
library(vegan)    # 多様性指標の計算
library(ggalluvial)
library(tools)

options(shiny.maxRequestSize = 50 * 1024^2 * 1000)
options(shiny.port = 8103)

function_file_lst <- list.files("functions", pattern = ".R$", full.names = TRUE)
for (i in function_file_lst) {
  source(i)
}

source('module/1_uploadCellranger.R')
source('module/2_reductionPlot.R')
source('module/3_1_geneExpression.R')
source('module/3_2_differentialGeneExpression.R')
source('module/4_1_clonotypeInformation.R')
source('module/4_2_summaryPlot.R')
source('module/4_3_uniqueClones.R')
source('module/4_4_geneUsage.R')
source('module/4_5_clonalLength.R')
source('module/4_6_diversity.R')
source('module/4_7_trackClonotype.R')
source('module/4_8_clonalOVerlap.R')
source('module/4_9_clonalProportion.R')
source('module/4_10_antigenPrediction.R')
source('module/4_11_phylogeneticTree.R')


# module_file_lst <- list.files("module", pattern = ".R$", full.names = TRUE)
# for (i in module_file_lst) {
#   source(i)
# }

source('utils.R')


ui <- navbarPage(

  includeCSS("style.css"),
  title = "SiCR: Web Application for Single Cell Repertoire Analysis",
  tags$footer("© 2023 Your Company. All rights reserved."),
  tabPanel(
    "Upload and Run",
   uploadCellrangerUI("upload"),
  ),
  tabPanel(
    "Dimentional Plot",
    reductionPlotUI("reduction_plot"),
  ),
  tabPanel(
    "Transcriptome",
    tabsetPanel(
      tabPanel(
        "Expression Visualization",
        geneExpressionUI("gene_expression"),
      ),
      tabPanel(
        "Differential Expression",
        differentialGeneExpressionUI("differential_gene_expression"),
      ),
    ),
  ),
  tabPanel(
    "Repertoire",
    tabsetPanel(
      tabPanel(
        "Clonotype Information",
        clonotypeInformationUI("clonotype_information"),
      ),
      tabPanel(
        "summaryPlot",
        summaryPlotUI("summary_plot"),
      ),
      tabPanel(
        "Unique Clones",
        uniqueClonesUI("unique_clones"),
      ),
      tabPanel(
        "Gene Usage",
        geneUsageUI("gene_usage"),
      ),
      tabPanel(
        "Clonal Length",
        cdrLengthDistUI("cdr_length_dist"),
      ),
      tabPanel(
        "Diversity",
        diversityAnalysisUI("diversity_analysis"),
      ),
      tabPanel(
        "Track Clonotype",
        trackClonotypeUI("track_clonotype"),
      ),
      tabPanel(
        "Clonal Overlap",
        clonalOverlapUI("clonal_overlap"),
      ),
      tabPanel(
        "Clonal Proportion",
        clonalProportionUI("clonal_proportion"),
      ),
      tabPanel(
        "Antigen Prediction",
        antigenPredictionUI("antigen_prediction"),
      ),
      tabPanel(
        "Phylogenetic Tree",
        phylogeneticTreeUI("phylogenetic_tree"),
      ),
    ),
),
)

server <- function(input, output) {
#  myReactives <- reactiveValues()
  load('filename.RData')
#   uploadCellrangerServer("upload", myReactives)
   reductionPlotServer("reduction_plot", myReactives)
   geneExpressionServer("gene_expression", myReactives)
   differentialGeneExpressionServer('differential_gene_expression', myReactives)
   clonotypeInformationServer("clonotype_information", myReactives)
   summaryPlotServer('summary_plot', myReactives)
   uniqueClonesServer('unique_clones', myReactives)
   geneUsageServer('gene_usage', myReactives)
   cdrLengthDistServer('cdr_length_dist', myReactives)
   diversityAnalysisServer('diversity_analysis', myReactives)
   trackClonotypeServer('track_clonotype', myReactives)
   clonalOverlapServer("clonal_overlap", myReactives)
   clonalProportionServer('clonal_proportion', myReactives)
   diversityAnalysisServer("diversity", myReactives)
   antigenPredictionServer("antigen_prediction", myReactives)
   phylogeneticTreeServer("phylogenetic_tree", myReactives)
#   trackClonotypeServer("track_clonotype_tcr", myReactives, "tcr")
# #  diversityServer("diversity_tcr", myReactives, "tcr")
# #  diversityServer("diversity_bcr", myReactives, "bcr")
#   # clonotypeInformationServer("clonotype_information_bcr", myReactives, "bcr")
#   # clonotypeInformationServer("clonotype_information_tcr", myReactives, "tcr")
#   # clonotypeInformationServer("clonotype_information_bcr", myReactives, "bcr")
#   clonotypeExpandServer("clonotype_expand_tcr", myReactives, "tcr")
#   clonotypeExpandServer("clonotype_expand_bcr", myReactives, "bcr")
#   uniqueClonesServer("unique_clones_tcr", myReactives, "tcr")
#   uniqueClonesServer("unique_clones_bcr", myReactives, "bcr")
#   clonalAbundanceServer("clonal_abundance_tcr", myReactives, "tcr")
#   clonalAbundanceServer("clonal_abundance_bcr", myReactives, "bcr")
#   clonalLengthServer("clonal_length_tcr", myReactives, "tcr")
#   clonalLengthServer("clonal_length_bcr", myReactives, "bcr")
#   clonalHomeostasisServer("clonal_homeostasis_tcr", myReactives, "tcr")
#   clonalHomeostasisServer("clonal_homeostasis_bcr", myReactives, "bcr")
#   clonalProportionServer("clonal_proportion_tcr", myReactives, "tcr")
#   clonalProportionServer("clonal_proportion_bcr", myReactives, "bcr")
#   vizGenesServer("vizgenes_tcr", myReactives, "tcr")
#   vizGenesServer("vizgenes_bcr", myReactives, "bcr")
#   clonalOverlapServer("clonal_overlap_tcr", myReactives, "tcr")
#   clonalOverlapServer("clonal_overlap_bcr", myReactives, "bcr")
#   integratedQCSubsetServer('reanalyze', myReactives)


  #  dataList <- uploadCellrangerServer("upload_cellranger")
  #  dataList <- uploadCellranger_mainServer("upload")

  #  seuratObject <- reactive({ dataList()[["seurat_object"]] })
  #  tcrList      <- reactive({ dataList()[["tcr_list"]] })
  #  bcrList      <- reactive({ dataList()[["bcr_list"]] })
  #  groupCols    <- reactive({ dataList()[["group_cols"]] })

  #  seuratObject()@misc$BCR <- bcrList()
  #  seuratObject()@misc$TCR <- tcrList()
  #  seuratObject()@misc$group <- groupCols()

  # output$downloaddata = downloadHandler(
  #   filename = "data.rds",
  #   content = function(file) {
  #     saveRDS(dataList, file)
  #   }
  # )

  #  geneExpressionUmap_mainServer("gene_expression_umap", seuratObject(), groupCols())
  #  geneExpressionBarplot_mainServer("gene_expression_barplot", seuratObject()@meta.data, groupCols())

  # alphaDiversity_mainServer("tcr_alpha_diversity", tcrList()[["TRB"]], groupCols())
  # alphaDiversity_mainServer("bcr_alpha_diversity", bcrList()[["IGH"]], groupCols())
  #
  # clonalAbundance_mainServer("tcr_clonal_abundance", tcrList()[["TRB"]], groupCols())
  # clonalAbundance_mainServer("bcr_clonal_abundance", bcrList()[["IGH"]], groupCols())
  #
  # clonotypeExpand_mainServer("tcr_clonotype_expand", tcrList()[["TRB"]], groupCols())
  # clonotypeExpand_mainServer("bcr_clonotype_expand", bcrList()[["IGH"]], groupCols())
  #
  # geneUsage_mainServer("tcr_gene_usage", tcrList(), groupCols())
  # geneUsage_mainServer("bcr_gene_usage", bcrList(), groupCols())
  #
  # antigenPrediction_mainServer("tcr_antigen_prediction", tcrList()[["TRB"]], "data/230323_ruft_TCR_antigen_database.tsv")
  # antigenPrediction_mainServer("bcr_antigen_prediction", bcrList()[["IGH"]], "data/230323_ruft_BCR_antigen_database.tsv")
  #
  # phylogeneticTree_mainServer("bcr_phylogenetic_tree", bcrList()[["IGH"]])
}



shinyApp(ui, server)
