# STEP 1: ベースイメージの選択
FROM --platform=linux/amd64 rocker/shiny-verse:latest

# STEP 2: システムライブラリのインストール
USER root
RUN apt-get update -o APT::Sandbox::User=root && \
    apt-get install -y --no-install-recommends \
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
    libglpk-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# STEP 3: Rパッケージのインストール

# --- 3.1 & 3.2 (変更なし) ---
# ★★★★★ ここに 'openxlsx' を追加 ★★★★★
RUN echo "---> Installing Shiny, Tidyverse, and plotting packages..." && \
    R -e "install.packages(c('shiny', 'shinyjs', 'shinyalert', 'shinyBS', 'shinythemes', 'shinybusy', 'shinyWidgets', 'DT', 'dplyr', 'tidyr', 'ggplot2', 'ggsci', 'ggpubr', 'ggrepel', 'patchwork', 'plotly', 'RColorBrewer', 'scales', 'stringr', 'forcats', 'data.table', 'readr', 'tibble', 'magrittr', 'htmltools', 'openxlsx', 'bslib', 'officer', 'rvg'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')" 


# --- 3.3: バイオインフォマティクス関連のCRANパッケージ ---
RUN echo "---> Installing Bioinformatics related packages from CRAN..." && \
    R -e "install.packages(c('ape', 'phangorn', 'circlize', 'vegan', 'HGNChelper', 'hdf5r', 'enrichR'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# --- 3.4: ggplot拡張とその他のユーティリティ ---
# ggtreeを削除
RUN echo "---> Installing ggplot extensions and other utilities..." && \
    R -e "install.packages(c('ggvenn', 'ggnewscale', 'ggtreeExtra', 'ggpointdensity', 'viridis', 'ggExtra', 'ggalluvial', 'stringdist', 'igraph', 'ggnetwork', 'jsonlite', 'httr', 'qs2', 'RANN'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')"


# --- 3.5: BioconductorとGitHubパッケージのインストーラー (変更なし) ---
RUN echo "---> Installing installer packages (BiocManager, remotes)..." && \
    R -e "install.packages(c('BiocManager', 'devtools', 'remotes', 'tools'), dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# --- 3.6: Bioconductor packages ---
RUN echo "---> Installing Bioconductor packages..." && \
    R -e "BiocManager::install(c('Seurat', 'harmony', 'SingleR', 'celldex', 'fgsea', 'slingshot', 'msigdbr', 'msa', 'Biostrings', 'scater', 'scran', 'ggtree', 'dowser', 'ComplexHeatmap'), update = FALSE, ask = FALSE, dependencies = TRUE)"

# --- 3.7: GitHubパッケージのインストール ---
RUN echo "---> Installing GitHub packages (scRepertoire, immunarch, presto)..." && \
    R -e "remotes::install_github('ncborcherding/scRepertoire', dependencies = TRUE, upgrade = 'never')" && \
    R -e "remotes::install_github('immunomind/immunarch', dependencies = TRUE, upgrade = 'never')" && \
    R -e "remotes::install_github('immunogenomics/presto', dependencies = TRUE, upgrade = 'never')"

# --- 3.8: CellChat v2 (jinworks/CellChat) ---
# ComplexHeatmap を Bioconductor から入れた後に CellChat をインストールすること
# sqjin/CellChat (v1) ではなく jinworks/CellChat (v2) を使う
RUN echo "---> Installing CellChat v2..." && \
    R -e "remotes::install_github('jinworks/CellChat', upgrade = 'never')"

# --- 3.9: PhantomJSのインストール ---
RUN echo "---> Installing phantomjs for webshot..." && \
    R -e "webshot::install_phantomjs()"

# STEP 4: アプリケーションファイルのコピー
# COPY app/ /srv/shiny-server/SiCR/
COPY app/ /srv/shiny-server/

# STEP 4b: Shiny Server タイムアウト設定（起動に時間がかかるため延長）
RUN printf 'run_as shiny;\n\napp_init_timeout 300;\napp_idle_timeout 0;\n\nserver {\n  listen 3838;\n  location / {\n    site_dir /srv/shiny-server;\n    log_dir /var/log/shiny-server;\n    directory_index on;\n  }\n}\n' > /etc/shiny-server/shiny-server.conf

# STEP 5: ポートの開放
EXPOSE 3838

# STEP 6: 起動コマンド
CMD ["/usr/bin/shiny-server"]