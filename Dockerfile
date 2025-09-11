# STEP 1: ベースイメージの選択
FROM --platform=linux/amd64 rocker/shiny-verse:latest

# STEP 2: システムライブラリのインストール
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    bzip2 \
    git \
    libgsl-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# STEP 3: Rパッケージのインストール

# --- 3.1 & 3.2 (変更なし) ---
# ★★★★★ ここに 'openxlsx' を追加 ★★★★★
RUN echo "---> Installing Shiny, Tidyverse, and plotting packages..." && \
    R -e "install.packages(c('shiny', 'shinyjs', 'shinyalert', 'shinyBS', 'shinythemes', 'shinybusy', 'shinyWidgets', 'DT', 'dplyr', 'tidyr', 'ggplot2', 'ggsci', 'ggpubr', 'ggrepel', 'patchwork', 'plotly', 'RColorBrewer', 'scales', 'stringr', 'forcats', 'data.table', 'readr', 'tibble', 'magrittr', 'htmltools', 'openxlsx'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')" 

# --- 3.3: バイオインフォマティクス関連のCRANパッケージ ---
RUN echo "---> Installing Bioinformatics related packages from CRAN..." && \
    R -e "install.packages(c('ape', 'phangorn', 'circlize', 'vegan', 'HGNChelper', 'hdf5r'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# --- 3.4: ggplot拡張とその他のユーティリティ ---
# ggtreeを削除
RUN echo "---> Installing ggplot extensions and other utilities..." && \
    R -e "install.packages(c('ggvenn', 'ggnewscale', 'ggtreeExtra', 'ggpointdensity', 'viridis', 'ggExtra', 'ggalluvial', 'ComplexHeatmap'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# --- 3.5: BioconductorとGitHubパッケージのインストーラー (変更なし) ---
RUN echo "---> Installing installer packages (BiocManager, remotes)..." && \
    R -e "install.packages(c('BiocManager', 'devtools', 'remotes', 'tools'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# --- 3.6: Bioconductorパッケージのインストール ---
RUN echo "---> Installing Bioconductor packages..." && \
    R -e "BiocManager::install(c('Seurat', 'msa', 'Biostrings', 'scater', 'scran', 'ggtree', 'dowser'), update = FALSE, ask = FALSE, dependencies = TRUE)"

# --- 3.7: GitHubパッケージのインストール (presto を追加) ---
RUN echo "---> Installing GitHub packages (scRepertoire, immunarch, presto)..." && \
    R -e "remotes::install_github('ncborcherding/scRepertoire', dependencies = TRUE, upgrade = 'never')" && \
    R -e "remotes::install_github('immunomind/immunarch', dependencies = TRUE, upgrade = 'never')" && \
    R -e "remotes::install_github('immunogenomics/presto', dependencies = TRUE, upgrade = 'never')"

# --- 3.8: PhantomJSのインストール (変更なし) ---
RUN echo "---> Installing phantomjs for webshot..." && \
    R -e "webshot::install_phantomjs()"

# STEP 4: アプリケーションファイルのコピー
# COPY SiCR/ /srv/shiny-server/SiCR/
COPY SiCR/ /srv/shiny-server/

# STEP 5: ポートの開放
EXPOSE 3838

# STEP 6: 起動コマンド
CMD ["/usr/bin/shiny-server"]