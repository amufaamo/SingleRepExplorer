# --- 必要なライブラリ (dplyrが読み込まれている前提) ---
# library(dplyr) # スクリプトの先頭やglobal.Rなどで読み込む

# --- ヘルパー関数定義 ---

# Group By を RadioButtons で更新する関数 (もし必要なら)
update_group_by <- function(session, myReactives) {
  print('update_group_by (for RadioButtons)') # 関数名を明確化
  # ★ dplyr::select と all_of() を使用 ★
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")
  metadatas <- tryCatch({ # エラーハンドリング追加
    myReactives$seurat_object@meta.data %>%
      dplyr::select(-all_of(minus_column)) %>%
      dplyr::select(!starts_with("TCR")) %>%
      dplyr::select(!starts_with("BCR"))
  }, error = function(e){
    warning("Error selecting metadata in update_group_by: ", e$message)
    return(NULL)
  })
  req(metadatas) # metadatasがNULLなら停止

  metadata_cols <- names(metadatas)
  validate(need(length(metadata_cols) > 0, "No suitable metadata columns found."))

  # リスト形式で選択肢を作成 (名前と値が同じ)
  group_cols <- setNames(as.list(metadata_cols), metadata_cols)

  # デフォルト選択を決定
  selected_value <- if ("sample" %in% metadata_cols) "sample" else metadata_cols[1]

  updateRadioButtons(session, "group_by", choices = group_cols, selected = selected_value)
}

# Group By を SelectInput で更新する関数
update_group_by_select_input <- function(session, myReactives) {
  print('update_group_by_select_input') # 関数名を明確化
  # ★ dplyr::select と all_of() を使用 ★
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")
  metadatas <- tryCatch({ # エラーハンドリング追加
    myReactives$seurat_object@meta.data %>%
      dplyr::select(-all_of(minus_column)) %>%
      dplyr::select(!starts_with("TCR")) %>%
      dplyr::select(!starts_with("BCR"))
   }, error = function(e){
    warning("Error selecting metadata in update_group_by_select_input: ", e$message)
    return(NULL)
  })
  req(metadatas) # metadatasがNULLなら停止

  metadata_cols <- names(metadatas)
  validate(need(length(metadata_cols) > 0, "No suitable metadata columns found."))

  # リスト形式で選択肢を作成 (名前と値が同じ)
  group_cols <- setNames(as.list(metadata_cols), metadata_cols)

  # デフォルト選択を決定
  selected_value <- if ("sample" %in% metadata_cols) "sample" else metadata_cols[1]

  updateSelectInput(session, "group_by", choices = group_cols, selected = selected_value)
}

# DimPlot 用に Group By を SelectInput で更新する関数 (TCR/BCR列を追加)
# ★ 重複していた定義を削除し、こちらに統一 ★
update_group_by_for_dimplot <- function(session, myReactives) {
  print('update_group_by_for_dimplot (for SelectInput)') # 関数名を明確化
  req(myReactives$seurat_object) # Seuratオブジェクトが必要

  # 除外する列名
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")

  # ★ dplyr::select と all_of() を使用 ★
  metadatas <- tryCatch({ # エラーハンドリング追加
      myReactives$seurat_object@meta.data %>%
        dplyr::select(-all_of(minus_column)) %>%
        dplyr::select(-starts_with("TCR_"), -starts_with("BCR_")) # TCR/BCR列を一旦除外（後で追加するため）
    }, error = function(e){
      warning("Error selecting metadata in update_group_by_for_dimplot: ", e$message)
      return(NULL) # エラー時は NULL を返す
    })
  req(metadatas) # metadatasがNULLなら停止

  metadata_cols <- names(metadatas)

  # TCR/BCRデータが存在する場合、特定の列を選択肢に追加
  # 注意: ここで追加する列名は myReactives$seurat_object@meta.data に実際に存在する必要があります。
  # もし addTCRClonotypeIdToSeuratObject などで追加されるならOK。
  tcr_cols_to_add <- c()
  bcr_cols_to_add <- c()
  if (!is.null(myReactives$tcr_df)) {
    # myReactives$tcr_df 自体からではなく、Seuratメタデータに追加されている想定の列名
    potential_tcr_cols <- c("TCR", "TCR_clonalFrequency", "TCR_cloneSize", "TCR_raw_clonotype_id") # 例
    tcr_cols_to_add <- intersect(potential_tcr_cols, names(myReactives$seurat_object@meta.data))
  }
  if (!is.null(myReactives$bcr_df)) {
    # myReactives$bcr_df 自体からではなく、Seuratメタデータに追加されている想定の列名
    potential_bcr_cols <- c("BCR", "BCR_clonalFrequency", "BCR_cloneSize", "BCR_raw_clonotype_id") # 例
    bcr_cols_to_add <- intersect(potential_bcr_cols, names(myReactives$seurat_object@meta.data))
  }
  metadata_cols <- unique(c(metadata_cols, tcr_cols_to_add, bcr_cols_to_add)) # 重複を除いて結合

  validate(need(length(metadata_cols) > 0, "No suitable metadata columns found after processing."))

  # グループ列をリスト形式で作成 (名前=値)
  group_cols <- setNames(as.list(metadata_cols), metadata_cols)

  # デフォルトの選択肢を設定
  default_selection <- if ("sample" %in% metadata_cols) "sample" else metadata_cols[1]
  if (!default_selection %in% metadata_cols) default_selection <- metadata_cols[1] # さらにフォールバック

  # SelectInput を更新
  updateSelectInput(session, "group_by", choices = group_cols, selected = default_selection)
}


# チェックボックスグループの選択肢を更新する関数
update_unique_group_choices <- function(session, myReactives, group_by_col) {
  req(myReactives$seurat_object, group_by_col, nzchar(group_by_col)) # group_by_col が空でないことも確認
  # group_by_col が実際にメタデータに存在するか確認
  validate(need(group_by_col %in% names(myReactives$seurat_object@meta.data),
                paste("Column '", group_by_col, "' not found in Seurat metadata.")))

  unique_groups <- tryCatch({
      unique(myReactives$seurat_object@meta.data[[group_by_col]])
  }, error = function(e) {
      warning("Error getting unique groups for '", group_by_col, "': ", e$message); NULL
  })
  req(unique_groups) # NULL なら停止

  # NAを除外
  unique_groups <- na.omit(unique_groups)
  validate(need(length(unique_groups) > 0, paste("No valid unique groups found in column '", group_by_col, "'.")))

  # 数値/文字列ソート (変更なし)
  if (is.factor(unique_groups)) { unique_groups <- levels(unique_groups) }
  if (all(!is.na(suppressWarnings(as.numeric(as.character(unique_groups)))))) { # Check if all can be numeric without NA warnings
    # Handle potential non-numeric coercion if original was factor/character mixing numbers and text
     numeric_representation <- suppressWarnings(as.numeric(as.character(unique_groups)))
     if(all(!is.na(numeric_representation))){ # If all successfully converted to numeric
         unique_groups <- as.character(sort(numeric_representation))
     } else { # Mixed types or non-numeric, sort as character
         unique_groups <- sort(as.character(unique_groups))
     }
  } else { unique_groups <- sort(as.character(unique_groups)) }


  # デフォルトで最初の値を選択
  selected_values <- if (length(unique_groups) > 0) unique_groups[1] else NULL

  updateCheckboxGroupInput(session, "unique_group", choices = unique_groups, selected = selected_values, inline = TRUE)
}


# マーカー遺伝子検索用に Group By と クラスター選択肢を更新する関数
update_group_by_for_marker <- function(session, input, myReactives) {
  req(myReactives$seurat_object)
  minus_column <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "barcode", "percent.mt", "RNA_snn_res.0.5")
  # ★ dplyr::select と all_of() を使用 ★
  metadatas <- tryCatch({
      myReactives$seurat_object@meta.data %>%
        dplyr::select(-all_of(minus_column)) %>%
        dplyr::select(!starts_with("TCR")) %>%
        dplyr::select(!starts_with("BCR"))
    }, error = function(e){ warning("Error selecting metadata in update_group_by_for_marker: ", e$message); NULL })
  req(metadatas)

  metadata_cols <- names(metadatas)
  validate(need(length(metadata_cols) > 0, "No suitable metadata columns found."))

  # 通常、クラスタリング結果 (seurat_clustersなど) やサンプルで比較する
  # ここでは全てのカテゴリカル/ファクター列を候補とする
  group_cols <- setNames(as.list(metadata_cols), metadata_cols)

  # デフォルトで 'seurat_clusters' があればそれを選択
  selected_group <- if ("seurat_clusters" %in% metadata_cols) "seurat_clusters" else metadata_cols[1]

  updateSelectInput(session, "group_by", choices = group_cols, selected = selected_group)

  # ★ input$group_by の値を使ってクラスター選択肢を更新 ★
  # observeEvent などで input$group_by の変更をトリガーする必要がある場合があるが、
  # この関数が呼ばれる時点で input$group_by が更新されている前提
  # req(input$group_by) # input$group_by をここで req すると循環参照になる可能性
  # 代わりに、更新された selected_group を直接使う
  current_group_by <- session$input$group_by %||% selected_group # session$input を使うか、引数で渡す

  if (!is.null(current_group_by) && current_group_by %in% colnames(myReactives$seurat_object@meta.data)) {
      cluster_choices <- tryCatch(unique(myReactives$seurat_object@meta.data[[current_group_by]]), error = function(e) NULL)
      if (!is.null(cluster_choices)) {
          cluster_choices <- sort(na.omit(cluster_choices)) # ソートしてNA除去
           # target_cluster と reference_cluster の選択肢を更新
           updateSelectInput(session, "target_cluster", choices = cluster_choices, selected = cluster_choices[1]) # 最初のクラスタを選択
           updateSelectInput(session, "reference_cluster", choices = cluster_choices, selected = NULL) # 最初は選択しない
      } else {
          warning("Could not get unique choices for target/reference cluster from: ", current_group_by)
           updateSelectInput(session, "target_cluster", choices = list("N/A" = ""), selected = "")
           updateSelectInput(session, "reference_cluster", choices = list("N/A" = ""), selected = "")
      }
  } else {
       warning("Selected group_by column '", current_group_by, "' not found or is NULL.")
       updateSelectInput(session, "target_cluster", choices = list("Select Group By First" = ""), selected = "")
       updateSelectInput(session, "reference_cluster", choices = list("Select Group By First" = ""), selected = "")
  }
}


# --- UIコンポーネント定義関数 ---

legendPositionInput <- function(ns, selected = "right") {
  selectInput(ns("legend"), "Legend Position", choices = c("right", "left", "bottom", "top", "none"), selected = selected)
}

plotWidthInput <- function(ns, value = 500, min = 100, max = 2000, step = 100) {
  # numericInputの方が使いやすい場合がある
  numericInput(ns("plot_width"), "Plot Width (px)", min = min, max = max, value = value, step = step)
}

plotHeightInput <- function(ns, value = 500, min = 100, max = 2000, step = 100) {
  numericInput(ns("plot_height"), "Plot Height (px)", min = min, max = max, value = value, step = step)
}

commonPlotOptions <- function(ns, legend_selected = "right", width_value = 500, height_value = 500) {
  tagList(
    legendPositionInput(ns, selected = legend_selected),
    plotWidthInput(ns, value = width_value),
    plotHeightInput(ns, value = height_value)
    # 他の共通オプションがあればここに追加
  )
}

reductionInput <- function(ns) {
  # 初期状態では選択肢を空にする。サーバー側 (update_reduction_choices) で動的に更新されるため。
  selectInput(ns("reduction"), "Reduction Method", choices = list("Loading..." = ""), selected = "")
}

groupByInput <- function(ns, choices = c("sample", "seurat_clusters"), selected = "sample") {
  # 初期選択肢 (サーバー側で更新される)
  selectInput(ns("group_by"), "Group by", choices = choices, selected = selected)
}

pointSizeInput <- function(ns, value = 0.1, min = 0.01, max = 10, step = 0.01) {
  numericInput(ns("point_size"), "Point Size", min = min, max = max, value = value, step = step)
}

labelSizeInput <- function(ns, value = 10, min = 0, max = 20, step = 1) {
  numericInput(ns("label_size"), "Label Size", min = min, max = max, value = value, step = step)
}

vdjType <- function(ns){
  # radioButtons の方が視覚的にわかりやすい場合も
  selectInput(ns('vdj_type'), label = 'VDJ Type', choices = c("TCR" = "tcr", "BCR" = "bcr"), selected = 'tcr')
}

valueType <- function(ns){ # clonalProportion などで使う想定
  radioButtons(ns("value_type"), label = "Value Type:", choices = c("Count" = "count", "Percentage (%)" = "percentage"), selected = "count")
}


# --- データ処理・結合関数 ---



h5_to_seurat_object <- function(myReactives) {
  req(myReactives$h5_path) # パスが必要
  h5 <- Seurat::Read10X_h5(myReactives$h5_path) # Seurat:: を明示
  # Gene Expression matrix がリストに含まれているかチェック
  if (is.list(h5) && !is.data.frame(h5) && ("Gene Expression" %in% names(h5))) {
    h5 <- h5[["Gene Expression"]]
  }
  myReactives$seurat_object <- Seurat::CreateSeuratObject(h5)
  return(myReactives)
}



# ★ dataframe_tcr/bcr は scRepertoire に依存 ★
# scRepertoire なしで実装する場合は、read.csv と dplyr で必要な列を整形する必要がある
# ここでは元の関数のまま残すが、scRepertoire がロードされている必要がある
dataframe_tcr <- function(myReactives){
  req(myReactives$tcr_path)
  df <- tryCatch(read.csv(myReactives$tcr_path), error = function(e) {warning("Error reading TCR CSV: ", e$message); NULL})
  req(df)
  # combineTCR は scRepertoire の関数
  if (!requireNamespace("scRepertoire", quietly = TRUE)) {
     warning("scRepertoire package is needed for dataframe_tcr function.")
     myReactives$tcr_df <- df # または最低限の処理
     return(myReactives)
  }
  # CombineTCR に渡す前に、サンプル名を抽出・追加する必要があるかもしれない
  # （元のコードでは addMetadataToTCR で後から結合していた）
  df_list <- df %>%
      mutate(sample_id = str_remove(barcode, "^.+-")) %>% # サンプル名を推定
      group_split(sample_id) %>% # サンプルごとに分割
      setNames(unique(.$sample_id)) # リストに名前を付ける

  # combineTCR にリストを渡す
  combined_list <- tryCatch(scRepertoire::combineTCR(df_list, removeNA = TRUE, filterMulti = TRUE),
                           error = function(e) { warning("combineTCR failed: ", e$message); NULL })
  req(combined_list)
  # 結果を myReactives に格納（形式は combineTCR の返り値による）
  myReactives$tcr_list <- combined_list # 例: リストとして格納
  # 必要であれば、ここからデータフレーム tcr_df を作成する処理を追加
  # myReactives$tcr_df <- ...
  return(myReactives)
}

dataframe_bcr <- function(myReactives){
  req(myReactives$bcr_path)
  df <- tryCatch(read.csv(myReactives$bcr_path), error = function(e) {warning("Error reading BCR CSV: ", e$message); NULL})
  req(df)
  if (!requireNamespace("scRepertoire", quietly = TRUE)) {
     warning("scRepertoire package is needed for dataframe_bcr function.")
     myReactives$bcr_df <- df
     return(myReactives)
  }
  # CombineBCR に渡す前にサンプル名を抽出・追加
  df_list <- df %>%
      mutate(sample_id = str_remove(barcode, "^.+-")) %>%
      group_split(sample_id) %>%
      setNames(unique(.$sample_id))

  # combineBCR にリストを渡す
  combined_list <- tryCatch(scRepertoire::combineBCR(df_list, removeNA = TRUE, filterMulti = TRUE),
                           error = function(e) { warning("combineBCR failed: ", e$message); NULL })
  req(combined_list)
  myReactives$bcr_list <- combined_list
  # 必要であれば bcr_df を作成
  # myReactives$bcr_df <- ...
  return(myReactives)
}

# --- メタデータ結合関数 (★ dplyr::select を使用 ★) ---

addTCRClonotypeIdToSeuratObject <- function(myReactives){
  req(myReactives$tcr_path, myReactives$seurat_object)
  tcr <- tryCatch(read.csv(myReactives$tcr_path), error=function(e) {warning("Failed to read TCR file in addTCRClonotypeId: ", e$message); NULL})
  req(tcr)
  # ★ dplyr::select を使用 ★
  tcr_select <- tcr %>%
    dplyr::select(any_of(c("barcode", "raw_clonotype_id", "exact_subclonotype_id"))) %>% # any_of で存在しない列でもエラーにしない
    distinct()
  # 必要な列が存在するか確認
  req("barcode" %in% names(tcr_select), "raw_clonotype_id" %in% names(tcr_select))

  # 列名を変更（存在しない列は変更されない）
  tcr_renamed <- tcr_select %>%
     rename(TCR_raw_clonotype_id = raw_clonotype_id)
  if ("exact_subclonotype_id" %in% names(tcr_renamed)) {
      tcr_renamed <- tcr_renamed %>% rename(TCR_exact_subclonotype_id = exact_subclonotype_id)
  }

  # 元のメタデータを行名基準で結合できるように準備
  seurat_meta <- myReactives$seurat_object@meta.data %>%
    tibble::rownames_to_column("barcode_rownames") # 元のrownamesを保持

  # barcode 列で結合 (barcode列がない場合は rownames で代用)
  join_col <- if ("barcode" %in% names(seurat_meta)) "barcode" else "barcode_rownames"
  if (!"barcode" %in% names(tcr_renamed)) {
      warning("TCR data does not have 'barcode' column for joining.")
      return(myReactives) # 結合できないのでそのまま返す
  }

  # 結合実行前に型を合わせる（推奨）
  tcr_renamed$barcode <- as.character(tcr_renamed$barcode)
  seurat_meta[[join_col]] <- as.character(seurat_meta[[join_col]])

  # left_join
  joined_meta <- dplyr::left_join(seurat_meta, tcr_renamed, by = setNames("barcode", join_col))

  # rownamesを戻す
  rownames(joined_meta) <- joined_meta$barcode_rownames
  joined_meta$barcode_rownames <- NULL # 不要になった列を削除

  myReactives$seurat_object@meta.data <- joined_meta
  return(myReactives)
}

# addClonotypeIdToBCR は bcr_df が combineBCR のリスト形式である前提のようだが、
# 他のコードでは tcr_df/bcr_df がデータフレームであることを期待しているため、
# この関数の役割と bcr_df の期待される形式を見直す必要がある。
# 一旦コメントアウトするか、bcr_df がデータフレームである前提で書き直す。
# addClonotypeIdToBCR <- function(myReactives){ ... }

addMetadataToTCR <- function(myReactives){
  req(myReactives$seurat_object, myReactives$tcr_df, myReactives$metadata) # 依存関係を明確に
  # ★ dplyr::select を使用 ★
  metadata_seurat <- myReactives$seurat_object@meta.data %>%
    tibble::rownames_to_column("barcode") %>% # rownamesをbarcode列に
    dplyr::select(any_of(c("sample", "seurat_clusters", "barcode"))) # 存在しない列は無視

  # tcr_df がデータフレームであることを確認
  validate(need(is.data.frame(myReactives$tcr_df), "myReactives$tcr_df is not a data frame."))
  validate(need("barcode" %in% names(myReactives$tcr_df), "tcr_df requires 'barcode' column."))

  S1 <- myReactives$tcr_df %>%
    # barcode が Seurat オブジェクトに存在するもののみにフィルタ
    dplyr::filter(barcode %in% metadata_seurat$barcode) %>%
    # Seurat メタデータ (sample, seurat_clusters) を結合
    dplyr::left_join(metadata_seurat, by = 'barcode') %>%
    # さらに外部メタデータ (myReactives$metadata) を sample で結合
    dplyr::left_join(myReactives$metadata, by = 'sample')

  myReactives$tcr_df <- S1
  return(myReactives)
}

addMetadataToBCR <- function(myReactives){
  req(myReactives$seurat_object, myReactives$bcr_df, myReactives$metadata)
  # ★ dplyr::select を使用 ★
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

addMetadata <- function(myReactives, vdj = 'tcr'){
  req(myReactives, vdj)
  file_path <- if(vdj == 'tcr') myReactives$tcr_path else myReactives$bcr_path
  req(file_path)

  df <- tryCatch(read.csv(file_path), error = function(e){ warning("Error reading VDJ file in addMetadata: ", e$message); NULL })
  req(df)
  validate(need("barcode" %in% names(df), "VDJ file requires 'barcode' column."))

  # ★ 元の delete_column は多くの情報を含むため、ここではサンプル抽出と重複削除に留める方が安全かも ★
  # metadata として、例えばサンプルごとの追加情報（治療前/後など）を想定
  # ここでは、元のコードの意図を尊重しつつ select を修正
  # ただし、この関数は 'sample' 列を生成し、他の情報を削除して myReactives$metadata に上書きするため、
  # 使い方によっては他の関数 (addMetadataToTCR/BCR) と競合する可能性あり。
  # 関数の目的（サンプルごとの外部メタデータを読み込む or VDJファイルからサンプル情報を抽出する）を明確にする必要がある。

  # VDJファイルからサンプル情報だけ抜き出す場合：
  df_meta_extracted <- df %>%
    mutate(sample = str_remove(barcode, "^.+-")) %>%
    # ★ dplyr::select と all_of(), starts_with() を使用 ★
    # select で指定した列 *以外* を削除する意図と解釈
    # dplyr::select(-all_of(delete_column), -starts_with("TCR_"), -starts_with("BCR_"))
    # ここではサンプル情報だけを抽出する方が安全か？
    dplyr::select(any_of(c("barcode", "sample"))) %>% # barcode と sample 列のみを選択 (存在すれば)
    distinct() # 重複削除

  validate(need("sample" %in% names(df_meta_extracted), "Could not extract 'sample' column."))
  print("Extracted metadata (head):")
  print(head(df_meta_extracted))

  # myReactives$metadata に格納（上書き）
  # もし既存のメタデータと結合したい場合は、このロジックを変更する必要がある
  myReactives$metadata <- df_meta_extracted
  return(myReactives)
}

# --- valueType UI 関数 ---
valueType <- function(ns){
  radioButtons(ns("value_type"),
    label = "Value Type:",
    choices = c("Count" = "count", "Percentage (%)" = "percentage"),
    selected = "count")
}