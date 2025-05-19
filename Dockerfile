FROM rocker/shiny:latest

# システムライブラリのインストール（Rパッケージに必要）
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libglpk-dev \
    libhdf5-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libxt-dev \
    libudunits2-dev \
    libfftw3-dev \
    libtiff-dev \
    libpng-dev \
    libjpeg-dev \
    libx11-dev \
    libglu1-mesa-dev \
    libfreetype6-dev \
    && rm -rf /var/lib/apt/lists/*

# 必要パッケージのインストール（CRANとBioconductor）
RUN R -e "install.packages(c('shiny', 'ggplot2', 'dplyr', 'shinythemes', 'ggsci', 'RColorBrewer', 'tidyverse', 'DT', 'ggpointdensity', 'viridis', 'ggExtra', 'ggpubr', 'shinyWidgets', 'patchwork', 'plotly', 'shinyjs', 'vegan', 'ggalluvial', 'tools'), repos = 'https://cloud.r-project.org')"

# Bioconductorパッケージのインストール
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); \
          BiocManager::install(c('Seurat', 'scRepertoire', 'immunarch', 'scater', 'scran'))"

# アプリケーションファイルをコピー
COPY . /srv/shiny-server/

# 権限の修正
RUN chown -R shiny:shiny /srv/shiny-server

# ポートを公開
EXPOSE 3838

# 起動は rocker/shiny ベースで自動

