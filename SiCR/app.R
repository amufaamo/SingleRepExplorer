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
library(openxlsx) # scTypeのマーカーファイル読み込みに必要
library(HGNChelper) # scTypeの遺伝子名変換に必要

options(shiny.maxRequestSize = 50 * 1024^2 * 1000)

# 関数の読み込み
function_file_lst <- list.files("functions", pattern = ".R$", full.names = TRUE)
for (i in function_file_lst) {
  source(i, local = TRUE)
}

# モジュールの読み込み
source('module/1_uploadCellranger.R', local = TRUE)
source('module/2_reductionPlot.R', local = TRUE)
source('module/3_1_geneExpression.R', local = TRUE)
source('module/3_2_differentialGeneExpression.R', local = TRUE)
source('module/4_1_clonotypeInformation.R', local = TRUE)
#source('module/4_2_summaryPlot.R', local = TRUE)
source('module/4_3_uniqueClones.R', local = TRUE)
# source('module/4_4_geneUsage.R', local = TRUE)
source('module/4_4_frequencyPlot.R', local = TRUE)
source('module/4_5_clonalLength.R', local = TRUE)
source('module/4_6_diversity.R', local = TRUE)
source('module/4_7_trackClonotype.R', local = TRUE)
source('module/4_8_clonalOverlap.R', local = TRUE)
source('module/4_9_clonalProportion.R', local = TRUE)
source('module/4_10_antigenPrediction.R', local = TRUE)
source('module/4_11_phylogeneticTree.R', local = TRUE)
source('module/4_12_clonalAbundance.R', local = TRUE)
source('module/4_13_publicClonotype.R', local = TRUE)

source('utils.R', local = TRUE)


# # UI定義
# ui <- navbarPage(
#   includeCSS("style.css"),
#   title = "SingleRepExplorer: Web Application for Single Cell Repertoire Analysis",
#   footer = shiny::HTML("<div style='text-align: center; padding: 10px;'>© 2025 SiCR. All rights reserved.</div>"),
  
#   tabPanel(
#     "Upload and Run",
#     uploadCellrangerUI("upload"),
#   ),
#   tabPanel(
#     "Dimensional Plot",
#     reductionPlotUI("reduction_plot"),
#   ),
#   tabPanel(
#     "Transcriptome",
#     tabsetPanel(
#       tabPanel(
#         "Expression Visualization",
#         geneExpressionUI("gene_expression"),
#       ),
#       tabPanel(
#         "Differential Expression",
#         differentialGeneExpressionUI("differential_gene_expression"),
#       ),
#     ),
#   ),
#   tabPanel(
#     "Repertoire",
#     tabsetPanel(
#       tabPanel(
#         "Clonotype Information",
#         clonotypeInformationUI("clonotype_information"),
#       ),
#       tabPanel(
#         "Summary Plot",
#         summaryPlotUI("summary_plot"),
#       ),
#       tabPanel(
#         "Clonal Abundance",
#         clonalAbundancePlotUI("clonal_abundance_plot"),
#       ),
#       # tabPanel(
#       #   "Unique Clones",
#       #   uniqueClonesUI("unique_clones"),
#       # ),
#       # tabPanel(
#       #   "Gene Usage",
#       #   geneUsageUI("gene_usage"),
#       # ),
#       tabPanel(
#         "Analysis Plot",
#         analysisPlotUI("analysis_plot"),
#       ),      
#       tabPanel(
#         "Analysis Plot",
#         analysisPlotUI("analysis_plot"),
#       ),      
#       tabPanel(
#         "Clonal Length",
#         cdrLengthDistUI("cdr_length_dist"),
#       ),
#       tabPanel(
#         "Diversity",
#         diversityAnalysisUI("diversity_analysis"),
#       ),
#       tabPanel(
#         "Track Clonotype",
#         trackClonotypeUI("track_clonotype"),
#       ),
#       tabPanel(
#         "Clonal Overlap",
#         clonalOverlapUI("clonal_overlap"),
#       ),
#       tabPanel(
#         "Public Clonotype",
#         publicClonotypeUI("public_clonotype"),
#       ),
#     #   tabPanel(
#     #     "Clonal Proportion",
#     #     clonalProportionUI("clonal_proportion"),
#     #   ),
#     #   tabPanel(
#     #     "Antigen Prediction",
#     #     antigenPredictionUI("antigen_prediction"),
#     #   ),
#     #   tabPanel(
#     #     "Phylogenetic Tree",
#     #     phylogeneticTreeUI("phylogenetic_tree"),
#     #   ),
#     ),
#   )
# )

ui <- navbarPage(
  includeCSS("style.css"),
  title = "SingleRepExplorer: Web Application for Single Cell Repertoire Analysis",
  footer = shiny::HTML("<div style='text-align: center; padding: 10px;'>© 2025 SiCR. All rights reserved.</div>"),
  
  tabPanel(
    "Upload and Run",
    uploadCellrangerUI("upload"),
  ),
  tabPanel(
    "Dimensional Plot",
    reductionPlotUI("reduction_plot"),
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
  
  # --- ここからが変更点です！ ---
  navbarMenu(
    "Repertoire",
    # "Basic Properties" グループ
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
    # tabPanel(
    #   "Gene Usage",
    #   # 元のコードでコメントアウトされていたものを有効にしました
    #   geneUsageUI("gene_usage")
    # ),
    
    # --- グループを分けるための区切り線 ---
    "----", 
    
    # "Clonotype comparison" グループ
    "Clonotype comparison",
    tabPanel(
      "Clonal Overlap", # ご希望の "overlap" に対応
      clonalOverlapUI("clonal_overlap")
    ),
    tabPanel(
      "Public Clonotype", # ご希望の "public clonotypes" に対応
      publicClonotypeUI("public_clonotype")
    ),
    tabPanel(
      "Track Clonotype",
      trackClonotypeUI("track_clonotype")
    )
  )
  # --- ここまでが変更点です！ ---
)

# サーバーロジック定義
server <- function(input, output, session) {
  # myReactivesを初期化
  myReactives <- reactiveValues()

  # データアップロード用のサーバーモジュールを有効化
  uploadCellrangerServer("upload", myReactives)

  # filename.RData をロードして初期データを設定する
  if (file.exists('filename.RData')) {
    load('filename.RData', envir = .GlobalEnv)
  }

  # 各解析タブのサーバーモジュールを呼び出す
  reductionPlotServer("reduction_plot", myReactives)
  geneExpressionServer("gene_expression", myReactives)
  differentialGeneExpressionServer('differential_gene_expression', myReactives)
  clonotypeInformationServer("clonotype_information", myReactives)
#  summaryPlotServer('summary_plot', myReactives)
  uniqueClonesServer('unique_clones', myReactives)
  clonalAbundancePlotServer('clonal_abundance_plot', myReactives)
#  geneUsageServer('gene_usage', myReactives)
#  analysisPlotServer('analysis_plot', myReactives)
  cdrLengthDistServer('cdr_length_dist', myReactives)
  diversityAnalysisServer('diversity_analysis', myReactives)
  trackClonotypeServer('track_clonotype', myReactives)
  clonalOverlapServer("clonal_overlap", myReactives)
  publicClonotypeServer("public_clonotype", myReactives)
  frequencyPlotServer('frequency_plot', myReactives)
  # clonalProportionServer('clonal_proportion', myReactives)
  # antigenPredictionServer("antigen_prediction", myReactives)
  # phylogeneticTreeServer("phylogenetic_tree", myReactives)
}

# アプリケーションの実行
shinyApp(ui, server)