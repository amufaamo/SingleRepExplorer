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
library(stringdist)  # Required for sequence similarity
library(msigdbr)    # Required for GSEA
library(fgsea)      # Required for GSEA
library(ggvenn)     # Required for public clonotypes
library(igraph)     # Required for network analysis
library(ggnetwork)  # Required for network analysis
library(scales)     # Required for plotting
library(officer)    # Required for PowerPoint export (editable/ungroupable)
library(rvg)        # Required for PowerPoint vector graphics (DrawingML)
if (!requireNamespace("qs2", quietly = TRUE)) install.packages("qs2", repos = "https://cloud.r-project.org/")
library(qs2)        # Required for fast session save/load (.qs2)

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
source('module/3_4_congaIntegration.R', local = TRUE)
source('module/4_16_cdr3Properties.R', local = TRUE)
source('module/4_10_antigenPrediction.R', local = TRUE)
source('module/4_9_clonalProportion.R', local = TRUE)

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
    ".app-navbar { box-shadow: 0 4px 12px rgba(0,0,0,0.06); background: white; position: sticky; top: 0; z-index: 1030; padding: 0 28px; height: 56px; display: flex; align-items: center; gap: 16px; border-bottom: 1px solid #e8eef4; }
    .app-navbar-brand { font-size: 1.3rem; font-weight: 800; color: #2B5B84; letter-spacing: -0.3px; text-decoration: none; display: flex; align-items: center; gap: 10px; }
    .app-navbar-subtitle { color: #aab4be; font-size: 0.85rem; margin-left: 4px; }
    .navbar-save-btn { margin-left: auto; background: #2B5B84 !important; color: white !important; border: none !important; font-size: 0.82rem !important; padding: 5px 13px !important; border-radius: 6px !important; white-space: nowrap; }
    .navbar-save-btn:hover { background: #1e4166 !important; }
    .navbar-save-btn:disabled { background: #8F9CA4 !important; cursor: not-allowed; }
    #qs2_save_spinner { font-size: 0.8rem; color: #2B5B84; white-space: nowrap; }
    .section-divider { margin: 36px 0 6px 0; padding: 13px 22px; background: linear-gradient(135deg, #374151 0%, #1f2937 100%); border-radius: 10px; color: white; font-size: 1.15rem; font-weight: 700; display: flex; align-items: center; gap: 12px; box-shadow: 0 4px 14px rgba(31,41,55,0.28); }
    .section-divider i { opacity: 0.9; }
    .module-card { margin-bottom: 16px; border-radius: 12px; box-shadow: 0 8px 24px rgba(0,0,0,0.08); border: none; overflow: hidden; transition: all 0.3s ease; }
    .module-card-header { background-color: #2B5B84; color: white; font-weight: 600; font-size: 1.15rem; padding: 15px 25px; display: flex; align-items: center; gap: 10px; cursor: pointer; position: relative; user-select: none; }
    .module-card-header i:first-child { margin-right: 8px; }
    .module-card-header::after {
      content: '\\f077';
      font-family: 'Font Awesome 5 Free', 'Font Awesome 6 Free', 'FontAwesome';
      font-weight: 900;
      margin-left: auto;
      transition: transform 0.3s ease;
    }
    .module-card-header[aria-expanded='false']::after {
      transform: rotate(180deg);
    }
    .card-body { padding: 30px; }
    .sidebar { background-color: #f8fbff; border-right: 1px solid #e1e8ed; padding: 20px; }
    /* Enhance sidebar panels used by modules */
    .well { background-color: #f8fbff !important; border: 1px solid #e1e8ed !important; border-radius: 8px !important; box-shadow: inset 0 1px 3px rgba(0,0,0,0.02) !important; padding: 20px !important; }

    @keyframes pulse-btn {
      0% { transform: scale(1); box-shadow: 0 4px 15px rgba(32, 139, 94, 0.4); }
      50% { transform: scale(1.02); box-shadow: 0 6px 20px rgba(32, 139, 94, 0.6); }
      100% { transform: scale(1); box-shadow: 0 4px 15px rgba(32, 139, 94, 0.4); }
    }
    .run-btn-pulse {
      animation: pulse-btn 2s infinite;
      background: linear-gradient(135deg, #208b5e 0%, #155d3e 100%) !important;
      border: none !important;
      color: white !important;
    }
    .run-btn-pulse:hover {
      background: linear-gradient(135deg, #1e7d54 0%, #114c32 100%) !important;
      transform: translateY(-2px) !important;
      box-shadow: 0 8px 25px rgba(32, 139, 94, 0.7) !important;
      animation: none !important;
    }
    .run-btn-pulse:active {
      transform: translateY(1px) !important;
      box-shadow: 0 2px 10px rgba(32, 139, 94, 0.4) !important;
    }
    /* Fix: bslib html-fill-container sets display:flex which overrides Bootstrap collapse display:none */
    .module-card .card-body.collapse:not(.show) {
      display: none !important;
      padding: 0 !important;
    }
    /* Fix: bslib sets grid-auto-rows:1fr which fixes each row height even after card collapse */
    .bslib-grid {
      grid-auto-rows: auto !important;
    }"
  )

# Helper: colored section divider
section_divider <- function(label, icon_name) {
  tags$div(
    class = "section-divider",
    icon(icon_name),
    label
  )
}

# UI Definition
ui <- page_fluid(
  theme = app_theme,
  lang = "en",

  tags$head(
    tags$title("SingleRepExplorer v1.1.0 - scRNA-seq & Immune Repertoire"),
    tags$meta(name = "description", content = "SingleRepExplorer v1.1.0: A web application for integrated single-cell RNA-seq and TCR/BCR repertoire analysis."),
    tags$script(HTML("
      // Bootstrap 5 native collapse — wire up each .module-card on page load
      $(function() {
        var idx = 0;
        document.querySelectorAll('.module-card').forEach(function(card) {
          // header: direct child with .card-header class (bslib adds it)
          var header = card.querySelector(':scope > .card-header');
          // body: first direct child that is NOT a .card-header
          var body = null;
          var children = card.children;
          for (var i = 0; i < children.length; i++) {
            if (!children[i].classList.contains('card-header')) {
              body = children[i];
              break;
            }
          }
          if (!header || !body) return;

          var id = 'mc-body-' + idx++;
          body.id = id;

          // Mark as Bootstrap 5 collapse (initially shown)
          body.classList.add('collapse', 'show');

          // Wire toggle attributes on header
          header.setAttribute('data-bs-toggle', 'collapse');
          header.setAttribute('data-bs-target', '#' + id);
          header.setAttribute('aria-expanded', 'true');
          header.setAttribute('aria-controls', id);

          // After collapse: reset card AND bslib-grid-item wrapper height
          body.addEventListener('hidden.bs.collapse', function() {
            card.style.height = 'auto';
            card.style.flex = '0 0 auto';
            card.style.minHeight = '0';
            var gridItem = card.closest('.bslib-grid-item');
            if (gridItem) {
              gridItem.style.height = 'auto';
              gridItem.style.minHeight = '0';
            }
          });

          // After expand: remove overrides and trigger plot redraw
          body.addEventListener('shown.bs.collapse', function() {
            card.style.height = '';
            card.style.flex = '';
            card.style.minHeight = '';
            var gridItem = card.closest('.bslib-grid-item');
            if (gridItem) {
              gridItem.style.height = '';
              gridItem.style.minHeight = '';
            }
            $(window).trigger('resize');
          });
        });
      });
    "))
  ),

  # ── Sticky top navbar ──────────────────────────────────────────────────────
  tags$div(
    class = "app-navbar",
    tags$span(class = "app-navbar-brand",
      icon("flask"), "SingleRepExplorer v1.1.0"
    ),
    tags$span(class = "app-navbar-subtitle",
      "Single-Cell RNA-seq & Immune Repertoire Analysis"
    ),
    tags$div(
      style = "margin-left: auto; display: flex; align-items: center; gap: 10px;",
      uiOutput("qs2_save_spinner"),
      downloadButton("download_session_qs2",
        label = tagList(icon("save"), " Save Session (.qs2)"),
        class = "navbar-save-btn"
      )
    )
  ),

  # ── Main content ───────────────────────────────────────────────────────────
  tags$div(
    style = "padding: 8px 28px 40px 28px;",

    layout_column_wrap(
      width = 1,

      # ── 1. Data Processing & UMAP ─────────────────────────────────────────
      section_divider("1. Data Processing & UMAP", "flask"),

      # uploadCellrangerUI generates Step 1 / Step 2 / Step 3 cards internally
      uploadCellrangerUI("upload"),

      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("project-diagram"), "Dimensional Plot & Subsetting"),
        reductionPlotUI("reduction_plot")
      ),

      # ── 2. Transcriptome Analysis ─────────────────────────────────────────
      section_divider("2. Transcriptome Analysis", "dna"),

      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("chart-bar"), "Expression Visualization"),
        geneExpressionUI("gene_expression")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("balance-scale"), "Differential Expression"),
        differentialGeneExpressionUI("differential_gene_expression")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("route"), "Trajectory & Pseudotime"),
        trajectoryInferenceUI("trajectory_inference")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("link"), "Phenotype-Repertoire Correlation"),
        phenotypeRepertoireUI("phenotype_repertoire")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("network-wired"), "Cell-Cell Communication"),
        cellCommunicationUI("cell_communication")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("layer-group"), "Multimodal Integration (CoNGA)"),
        congaIntegrationUI("conga_integration")
      ),

      # ── 3. Repertoire Analysis ────────────────────────────────────────────
      section_divider("3. Repertoire Analysis", "viruses"),

      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("info-circle"), "Clonotype Information"),
        clonotypeInformationUI("clonotype_information")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("fingerprint"), "Unique Clonotypes"),
        uniqueClonesUI("unique_clones")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("chart-pie"), "Frequency Analysis"),
        frequencyPlotUI("frequency_plot")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("sort-amount-down"), "Clonal Abundance"),
        clonalAbundancePlotUI("clonal_abundance_plot")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("leaf"), "Diversity"),
        diversityAnalysisUI("diversity_analysis")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("ruler-horizontal"), "Clonal Length"),
        cdrLengthDistUI("cdr_length_dist")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("flask"), "CDR3 Physicochemical Profiling"),
        cdr3PropertiesUI("cdr3_properties")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("object-group"), "Clonal Overlap"),
        clonalOverlapUI("clonal_overlap")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("users"), "Public Clonotype"),
        publicClonotypeUI("public_clonotype")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("project-diagram"), "Sequence Similarity"),
        sequenceSimilarityUI("sequence_similarity")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("route"), "Track Clonotype"),
        trackClonotypeUI("track_clonotype")
      ),
      card(
        class = "module-card",
        card_header(class = "module-card-header", icon("tree"), "BCR Lineage Tree"),
        phylogeneticTreeUI("phylogenetic_tree")
      )
      # BCR SHM & Isotype module hidden for v2.2 release (column-structure issues with
      # combineBCR output). Re-enable once bcr_df column naming is fully stabilized.
      # card(
      #   class = "module-card",
      #   card_header(class = "module-card-header", icon("code-branch"), "BCR SHM & Isotype"),
      #   shmIsotypeUI("shm_isotype")
      # )
    )
  ),

  # ── Footer ─────────────────────────────────────────────────────────────────
  tags$div(
    style = "text-align: center; padding: 20px; background: #f8fbff; border-top: 1px solid #e1e8ed; color: #666; font-size: 0.9em;",
    "© 2026 SingleRepExplorer v1.1.0 - Open-source integrated scRNA-seq & immune repertoire analysis for immunologists."
  )
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

  # ── qs2 Session Download ────────────────────────────────────────────────────
  output$qs2_save_spinner <- renderUI({ NULL })

  output$download_session_qs2 <- downloadHandler(
    filename = function() {
      paste0("singlerepexplorer_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".qs2")
    },
    content = function(file) {
      validate(need(!is.null(myReactives$seurat_object),
        "No analysis data found. Please run the analysis pipeline first."))
      notif_id <- showNotification(
        tagList(icon("spinner", class = "fa-spin"), " Saving session to .qs2 ... This may take a moment."),
        duration = NULL, type = "message", id = "qs2_save_notif"
      )
      on.exit(removeNotification("qs2_save_notif"), add = TRUE)
      session_data <- reactiveValuesToList(myReactives)
      session_data$raw_qc_seurat <- NULL  # 未処理QCデータは除外してファイルサイズを削減
      qs2::qs_save(session_data, file)
    }
  )

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
  # shmIsotypeServer("shm_isotype", myReactives)  # hidden for v2.2 release
  # subsetServer is integrated into reductionPlotServer (no standalone UI needed)
  congaIntegrationServer("conga_integration", myReactives)
  cdr3PropertiesServer("cdr3_properties", myReactives)
}

# Run the application
shinyApp(ui, server)
