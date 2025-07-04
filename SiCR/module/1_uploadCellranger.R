uploadCellrangerUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(6,
      offset = 1,
      # shinyjs を有効にします
      shinyjs::useShinyjs(),
      
      tags$div(
        style = "border: 0.5px solid #000; padding: 10px;",
        tags$p("Please upload the filtered_feature_bc_matrix.h5 file generated by cellranger. You can optionally upload the T cell receptor file (filtered_contig_annotations.csv) and/or the B cell receptor file (filtered_contig_annotations.csv)."),
        tags$p("After uploading, click the Run button.")
      ),
      tags$div(style = "margin-top: 10px;"),
      
      # --- ファイルアップロードUIを改善 ---
      # 横並びにして、ステータス表示用のuiOutputを追加
      fluidRow(
        column(8, fileInput(ns("h5"), "1. Choose .h5 file")),
        column(4, style = "margin-top: 25px;", uiOutput(ns("h5_status_ui"))) # ステータス表示
      ),
      fluidRow(
        column(8, fileInput(ns("tcr"), "2. Choose .tcr file (optional)")),
        column(4, style = "margin-top: 25px;", uiOutput(ns("tcr_status_ui"))) # ステータス表示
      ),
      fluidRow(
        column(8, fileInput(ns("bcr"), "3. Choose .bcr file (optional)")),
        column(4, style = "margin-top: 25px;", uiOutput(ns("bcr_status_ui"))) # ステータス表示
      ),
      
      actionButton(ns("run"), "Run"),
      
      # --- 実行ステータスを表示するエリア ---
      # 最初は隠しておきます
      shinyjs::hidden(
        tags$div(id = ns("status_message"),
                 style = "margin-top: 15px; padding: 10px; border-radius: 5px; background-color: #f0f0f0;")
      )
    )
  )
}


uploadCellrangerServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- 1. ファイルアップロード完了を通知 ---
    output$h5_status_ui <- renderUI({
      req(input$h5) # ファイルがアップロードされたら...
      tags$p("✅ Uploaded!", style = "color: green;")
    })
    
    output$tcr_status_ui <- renderUI({
      req(input$tcr) # ファイルがアップロードされたら...
      tags$p("✅ Uploaded!", style = "color: green;")
    })
    
    output$bcr_status_ui <- renderUI({
      req(input$bcr) # ファイルがアップロードされたら...
      tags$p("✅ Uploaded!", style = "color: green;")
    })


    # --- 2. Runボタンが押された後の処理をここに集約 ---
    observeEvent(input$run, {
      
      # h5ファイルは必須なので、アップロードされているかチェック
      if (is.null(input$h5$datapath)) {
        shinyalert::shinyalert("Oops!", "Please upload the .h5 file first.", type = "error")
        return() # 処理を中断
      }
      
      # --- 3. 解析開始をユーザーに通知 ---
      shinyjs::disable("run") # Runボタンを無効化
      shinyjs::show("status_message") # ステータス表示エリアを表示
      # HTMLを使ってメッセージを動的に変更！
      shinyjs::html(id = "status_message", 
                    html = "<p style='color: blue;'>🏃 Run button clicked. Starting analysis...</p>")
      
      # --- 時間のかかる解析処理 ---
      # withProgressで進行状況を表示します
      withProgress(message = 'Analysis in progress', value = 0, {
        
        # 処理を進めるごとにプログレスバーを進めます
        incProgress(0.1, detail = "Assigning file paths...")
        myReactives <- fileUpload(input, myReactives)

        incProgress(0.2, detail = "Loading .h5 file...")
        myReactives <- h5_to_seurat_object(myReactives)
        
        incProgress(0.2, detail = "Running Seurat object...")
        myReactives <- run_seurat_object(myReactives)
        
        incProgress(0.2, detail = "Adding cell types...")
        myReactives <- addCelltype(myReactives)
        
        # TCRファイルの処理 (ファイルがアップロードされていれば)
        if (!is.null(myReactives$tcr_path)) {
          incProgress(0.1, detail = "Processing TCR data...")
          req(myReactives$seurat_object)
          myReactives$tcr_df <- tcr_csv_to_dataframe(myReactives$tcr_path)
          # --- DEBUG START (1_uploadCellranger.R - TCR processing) ---
          print(paste0("DEBUG: TCR: tcr_df after tcr_csv_to_dataframe - nrow: ", nrow(myReactives$tcr_df)))
          print(paste0("DEBUG: TCR: tcr_df after tcr_csv_to_dataframe - names: ", paste(names(myReactives$tcr_df), collapse = ", ")))
          
          # Get Seurat object barcodes (cell names)
          seurat_barcodes <- rownames(myReactives$seurat_object@meta.data)

          if (nrow(myReactives$tcr_df) > 0) {
            print(paste0("DEBUG: TCR: First 5 barcodes from tcr_df: ", paste(head(myReactives$tcr_df$barcode, 5), collapse = ", ")))
          }
          print(paste0("DEBUG: TCR: First 5 barcodes from seurat_object (rownames): ", paste(head(seurat_barcodes, 5), collapse = ", ")))
          
          # Standardize barcodes by removing '-1' suffix for matching
          # This is a common point of mismatch between Cell Ranger outputs and Seurat objects
          tcr_barcodes_clean <- stringr::str_replace(myReactives$tcr_df$barcode, "-1$", "")
          seurat_barcodes_clean <- stringr::str_replace(seurat_barcodes, "-1$", "")
          
          initial_overlap_count <- length(intersect(tcr_barcodes_clean, seurat_barcodes_clean))
          print(paste0("DEBUG: TCR: Number of overlapping barcodes (cleaned) before filter: ", initial_overlap_count))
          
          # Filter tcr_df using cleaned barcodes
          # Temporarily add cleaned barcode to tcr_df for filtering
          myReactives$tcr_df <- myReactives$tcr_df %>%
            dplyr::mutate(barcode_clean = stringr::str_replace(barcode, "-1$", "")) %>%
            dplyr::filter(barcode_clean %in% seurat_barcodes_clean) %>%
            dplyr::select(-barcode_clean) # Remove temporary column

          print(paste0("DEBUG: TCR: tcr_df after filter (cleaned barcodes) - nrow: ", nrow(myReactives$tcr_df)))

          metadata <- myReactives$seurat_object@meta.data %>%
            tibble::rownames_to_column(var = "barcode") %>% 
            dplyr::select(sample, seurat_clusters, barcode)
          myReactives$tcr_df <- dplyr::left_join(myReactives$tcr_df, metadata, by = "barcode")
          # --- DEBUG END (1_uploadCellranger.R - TCR processing) ---
          print(paste0("DEBUG: TCR: tcr_df after filter and left_join - nrow: ", nrow(myReactives$tcr_df))) # This line was already there, keeping it.
          print(paste0("DEBUG: TCR: tcr_df after filter and left_join - names: ", paste(names(myReactives$tcr_df), collapse = ", "))) # This line was already there, keeping it.
          saveRDS(myReactives$tcr_df, 'tcr_df.rds')
          write.csv(myReactives$tcr_df, 'tcr_df.csv', row.names = FALSE, quote = FALSE)
        }
        
        # BCRファイルの処理 (ファイルがアップロードされていれば)
        if (!is.null(myReactives$bcr_path)) {
          incProgress(0.1, detail = "Processing BCR data...")
          req(myReactives$seurat_object)
          myReactives$bcr_df <- bcr_csv_to_dataframe(myReactives$bcr_path)
          # --- DEBUG START (1_uploadCellranger.R - BCR processing) ---
          print(paste0("DEBUG: BCR: bcr_df after bcr_csv_to_dataframe - nrow: ", nrow(myReactives$bcr_df)))
          print(paste0("DEBUG: BCR: bcr_df after bcr_csv_to_dataframe - names: ", paste(names(myReactives$bcr_df), collapse = ", ")))

          # Get Seurat object barcodes (cell names)
          seurat_barcodes <- rownames(myReactives$seurat_object@meta.data)

          if (nrow(myReactives$bcr_df) > 0) {
            print(paste0("DEBUG: BCR: First 5 barcodes from bcr_df: ", paste(head(myReactives$bcr_df$barcode, 5), collapse = ", ")))
          }
          print(paste0("DEBUG: BCR: First 5 barcodes from seurat_object (rownames): ", paste(head(seurat_barcodes, 5), collapse = ", ")))

          # Standardize barcodes by removing '-1' suffix for matching
          # This is a common point of mismatch between Cell Ranger outputs and Seurat objects
          bcr_barcodes_clean <- stringr::str_replace(myReactives$bcr_df$barcode, "-1$", "")
          seurat_barcodes_clean <- stringr::str_replace(seurat_barcodes, "-1$", "")

          initial_overlap_count_bcr <- length(intersect(bcr_barcodes_clean, seurat_barcodes_clean))
          print(paste0("DEBUG: BCR: Number of overlapping barcodes (cleaned) before filter: ", initial_overlap_count_bcr))

          # Filter bcr_df using cleaned barcodes
          # Temporarily add cleaned barcode to bcr_df for filtering
          myReactives$bcr_df <- myReactives$bcr_df %>%
            dplyr::mutate(barcode_clean = stringr::str_replace(barcode, "-1$", "")) %>%
            dplyr::filter(barcode_clean %in% seurat_barcodes_clean) %>%
            dplyr::select(-barcode_clean) # Remove temporary column

          print(paste0("DEBUG: BCR: bcr_df after filter (cleaned barcodes) - nrow: ", nrow(myReactives$bcr_df)))
          # --- DEBUG END (1_uploadCellranger.R - BCR processing) ---
          metadata <- myReactives$seurat_object@meta.data %>%
            tibble::rownames_to_column(var = "barcode") %>%
            dplyr::select(sample, seurat_clusters, barcode)
          myReactives$bcr_df <- dplyr::left_join(myReactives$bcr_df, metadata, by = "barcode")
          saveRDS(myReactives$bcr_df, 'bcr_df.rds')
          write.csv(myReactives$bcr_df, 'bcr_df.csv',  row.names = FALSE, quote = FALSE)
          save(myReactives, file = "filename.RData")
        }
        
        incProgress(0.1, detail = "Finalizing...")
        Sys.sleep(1) # ちょっと待って完了メッセージを見せる演出
      }) # withProgress 終了

      # --- 4. 解析完了をユーザーに通知 ---
      shinyjs::html(id = "status_message", 
                    html = "<p style='color: green;'>🎉 Analysis finished successfully!</p>")
      shinyjs::enable("run") # Runボタンを再度有効化
      
      # 5秒後に完了メッセージを消す（お好みで）
      shinyjs::delay(5000, shinyjs::hide("status_message"))
    })
    
    # 元のリアクティブな値の変更を監視するobserveEventは不要になります
    # なぜなら、全ての処理が "Run" ボタン起点で順番に実行されるようになったからです！
    # observeEvent(myReactives$h5_path, { ... })
    # observeEvent(myReactives$tcr_path, { ... })
    # observeEvent(myReactives$bcr_path, { ... })
    
  })
}

# この関数は変更なしでOKです！
fileUpload <- function(input, myReactives) {
  myReactives$h5_path <- input$h5$datapath
  myReactives$tcr_path <- input$tcr$datapath
  myReactives$bcr_path <- input$bcr$datapath
  return(myReactives)
}


# (パッケージとしてビルドする場合の推奨記述)
#' @importFrom dplyr %>% mutate case_when relocate starts_with coalesce across ends_with na_if select all_of any_of filter distinct full_join rename_with rename everything
#' @importFrom stringr str_c str_replace_all
#' @importFrom tidyselect everything
#' @importFrom readr read_csv
#' @importFrom scRepertoire combineTCR
#' @importFrom stats setNames
#' @importFrom rlang sym .data

# --- 定数定義 ---
PREFIX_TCR_PAIR  <- "TCR_pair_"
PREFIX_TCR_ALPHA <- "TCR_TRA_" # Alpha鎖 (TRA)
PREFIX_TCR_BETA  <- "TCR_TRB_" # Beta鎖 (TRB)

PREFIX_BCR_PAIR  <- "BCR_pair_" # Assuming this constant is needed for BCR processing


# --- ヘルパー関数定義 ---
# (csv_to_tcr_pair_dataframe, csv_to_tcr_tra_dataframe, csv_to_tcr_trb_dataframe, merge_tcr, paste_existing_cols は変更なし)
# (省略のため、ここでは再掲しません)
# --- (省略されたヘルパー関数) ---
#' @title CSVからTCRペア鎖データを抽出・整形
#' @description 指定されたCSVファイルからscRepertoire::combineTCRを使用して
#'              ペア鎖データ (主にTRA/TRB) を抽出し、整形します。
#'              同一barcodeに複数α/β鎖がある場合はフィルタリングします (filterMulti=TRUE)。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @param sample_name `character(1)`. combineTCRで使用する一時的なサンプル名。デフォルトは "temp_sample"。
#' @param id_name `character(1)`. combineTCRで使用する一時的なID名。デフォルトは "temp_id"。
#' @return `data.frame`. 整形されたTCRペア鎖情報。列名は `PREFIX_TCR_PAIR` 付与。
csv_to_tcr_pair_dataframe <- function(csv_path,
                                      sample_name = "temp_sample",
                                      id_name = "temp_id") {

  barcode_prefix <- paste0(sample_name, "_", id_name, "_")
  tcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  if (is.null(tcr_raw)) return(NULL)

  pair_list <- tryCatch({
    scRepertoire::combineTCR(tcr_raw, samples = sample_name, ID = id_name, filterMulti = TRUE, removeNA = TRUE)
  }, error = function(e) {
    stop("Error in combineTCR: ", e$message)
    return(NULL)
  })

  if (is.null(pair_list) || length(pair_list) == 0 || is.null(pair_list[[1]])) {
    warning("combineTCR did not return valid pair data.")
    return(data.frame()) # 空のデータフレームを返す
  }
  pair <- pair_list[[1]]

  pair <- pair %>%
    dplyr::mutate(barcode = stringr::str_replace_all(barcode, pattern = barcode_prefix, replacement = "")) %>%
    dplyr::select(-dplyr::any_of(c("sample", "ID"))) %>%
    dplyr::rename_with(~ stringr::str_c(PREFIX_TCR_PAIR, .), dplyr::everything())

  return(pair)
}

# Assuming csv_to_bcr_pair_dataframe exists in this file,
# analogous to csv_to_tcr_pair_dataframe, and contains the problematic call.
# The following is a hypothetical reconstruction to show the fix.
#' @title CSVからBCRペア鎖データを抽出・整形
#' @description 指定されたCSVファイルからscRepertoire::combineBCRを使用して
#'              ペア鎖データ (主にIGH/IGL/IGK) を抽出し、整形します。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @param sample_name `character(1)`. combineBCRで使用する一時的なサンプル名。デフォルトは "temp_sample"。
#' @param id_name `character(1)`. combineBCRで使用する一時的なID名。デフォルトは "temp_id"。
#' @return `data.frame`. 整形されたBCRペア鎖情報。列名は `PREFIX_BCR_PAIR` 付与。
csv_to_bcr_pair_dataframe <- function(csv_path,
                                      sample_name = "temp_sample",
                                      id_name = "temp_id") {

  barcode_prefix <- paste0(sample_name, "_", id_name, "_")
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  if (is.null(bcr_raw)) return(NULL)

  pair_list <- tryCatch({
    # FIX: combineBCR does NOT have filterMulti = TRUE. Removed it.
    scRepertoire::combineBCR(bcr_raw, samples = sample_name, ID = id_name, removeNA = TRUE)
  }, error = function(e) {
    stop("Error in combineBCR: ", e$message)
    return(NULL)
  })

  if (is.null(pair_list) || length(pair_list) == 0 || is.null(pair_list[[1]])) {
    warning("combineBCR did not return valid pair data.")
    return(data.frame()) # 空のデータフレームを返す
  }
  pair <- pair_list[[1]]

  pair <- pair %>%
    dplyr::mutate(barcode = stringr::str_replace_all(barcode, pattern = barcode_prefix, replacement = "")) %>%
    dplyr::select(-dplyr::any_of(c("sample", "ID"))) %>%
    dplyr::rename_with(~ stringr::str_c(PREFIX_BCR_PAIR, .), dplyr::everything())

  return(pair)
}

#' @title CSVからTCR Alpha鎖(TRA)データを抽出・整形
#' @description 指定されたCSVファイルからTRA鎖のデータを抽出し、
#'              barcodeごとにユニークにして整形します。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 整形されたTCR TRA情報。列名は `PREFIX_TCR_ALPHA` 付与。
csv_to_tcr_tra_dataframe <- function(csv_path){

  tcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  if (is.null(tcr_raw) || !"chain" %in% names(tcr_raw) || !"barcode" %in% names(tcr_raw)) {
    warning("CSV file must contain 'chain' and 'barcode' columns.")
    return(data.frame()) # エラーではなく空を返す
  }

  TRA <- tcr_raw %>%
    dplyr::filter(chain == 'TRA') %>%
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    dplyr::rename_with(~ stringr::str_c(PREFIX_TCR_ALPHA, .), dplyr::everything())

  return(TRA)
}

#' @title CSVからTCR Beta鎖(TRB)データを抽出・整形
#' @description 指定されたCSVファイルからTRB鎖のデータを抽出し、
#'              barcodeごとにユニークにして整形します。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 整形されたTCR TRB情報。列名は `PREFIX_TCR_BETA` 付与。
csv_to_tcr_trb_dataframe <- function(csv_path){

  tcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  if (is.null(tcr_raw) || !"chain" %in% names(tcr_raw) || !"barcode" %in% names(tcr_raw)) {
    warning("CSV file must contain 'chain' and 'barcode' columns.")
    return(data.frame()) # エラーではなく空を返す
  }

  TRB <- tcr_raw %>%
    dplyr::filter(chain == 'TRB') %>%
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    dplyr::rename_with(~ stringr::str_c(PREFIX_TCR_BETA, .), dplyr::everything())

  return(TRB)
}

#' @title TCR各鎖データを結合
#' @description ペア鎖、Alpha鎖(TRA)、Beta鎖(TRB)のデータフレームをbarcodeをキーにして結合します。
#' @param pair `data.frame`. `csv_to_tcr_pair_dataframe` の出力。
#' @param TRA `data.frame`. `csv_to_tcr_tra_dataframe` の出力。
#' @param TRB `data.frame`. `csv_to_tcr_trb_dataframe` の出力。
#' @return `data.frame`. 結合されたTCRデータフレーム。
merge_tcr <- function(pair, TRA, TRB){

  pair_barcode_col <- paste0(PREFIX_TCR_PAIR, "barcode")
  tra_barcode_col  <- paste0(PREFIX_TCR_ALPHA, "barcode")
  trb_barcode_col  <- paste0(PREFIX_TCR_BETA, "barcode")

  # pair が空なら空を返す
  if (nrow(pair) == 0) {
      warning("Input 'pair' data frame is empty in merge_tcr. Returning empty data frame.")
      return(data.frame())
  }
  # TRA または TRB が空の場合でも結合は可能 (full_joinのため)
  if (nrow(TRA) == 0) warning("Input 'TRA' data frame is empty in merge_tcr.")
  if (nrow(TRB) == 0) warning("Input 'TRB' data frame is empty in merge_tcr.")


  tcr_merged <- dplyr::full_join(pair, TRA, by = stats::setNames(tra_barcode_col, pair_barcode_col))
  tcr_merged <- dplyr::full_join(tcr_merged, TRB, by = stats::setNames(trb_barcode_col, pair_barcode_col))

  return(tcr_merged)
}

#' @title 存在する列の値を連結する内部ヘルパー関数
#' @description データフレームと列名のベクトルを受け取り、存在する列の値を行ごとに連結。
#' @param df `data.frame`. 対象のデータフレーム。
#' @param cols `character`. 連結したい列名のベクトル。
#' @return `character`. 各行について連結された文字列ベクトル。
paste_existing_cols <- function(df, cols) {
  existing_cols <- intersect(cols, names(df))
  if (length(existing_cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  # apply を使うより rowwise + paste0/unite の方が dplyr 的かもしれないが、
  # NA処理を含めると apply の方が簡潔な場合もある
  apply(df[, existing_cols, drop = FALSE], 1, function(row_values) {
      # NA を空文字に変換してから結合
      row_values_no_na <- ifelse(is.na(row_values), "", row_values)
      paste0(row_values_no_na, collapse = "")
  })
}
# --- メイン関数定義 (修正版) ---

#' @title CSVからTCRデータを処理するメイン関数 (代表ID修正・ペア必須版)
#' @description 指定されたCSVファイルからTCRデータ (主にTRA/TRB) を読み込み、
#'              ペア鎖、Alpha鎖、Beta鎖に分割・整形した後、結合し、TRB鎖のIDを
#'              元に代表ID (`raw_clonotype_id` と `exact_subclonotype_id`) を格納し、
#'              全長配列を生成し、列順序を整え、barcode列名を変更します。
#'              `raw_consensus_id` は削除されます。
#'              `exact_subclonotype_id` は TRB の `raw_clonotype_id` と
#'              `exact_subclonotype_id` を `_` で連結した値になります。
#'              **最終的にTRA鎖とTRB鎖の両方が存在する行のみを保持します。**
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 全ての処理が完了し、TRA/TRBペアが存在する行のみを含むTCRデータフレーム。
#'         処理中にエラーが発生した場合や、有効なペアデータが見つからない場合は、
#'         空のデータフレームまたは警告付きの部分的な結果を返すことがあります。
tcr_csv_to_dataframe <- function(csv_path){

  # --- 1. 入力チェック ---
  if (!file.exists(csv_path)) {
    stop("Input CSV file not found: ", csv_path)
    return(NULL)
  }

  # --- 2. データ抽出・整形 ---
  pair_data <- csv_to_tcr_pair_dataframe(csv_path)
  tra_data  <- csv_to_tcr_tra_dataframe(csv_path)
  trb_data  <- csv_to_tcr_trb_dataframe(csv_path)

  if (is.null(pair_data) || nrow(pair_data) == 0 || !paste0(PREFIX_TCR_PAIR, "barcode") %in% names(pair_data)) {
      warning("No valid paired TCR data found or extracted. Returning an empty data frame.")
      return(data.frame())
  }
  if (is.null(tra_data) || nrow(tra_data) == 0) {
      warning("No TRA data extracted. Final result might be empty after filtering.")
  }
  if (is.null(trb_data) || nrow(trb_data) == 0) {
      warning("No TRB data extracted. Final result might be empty after filtering.")
  }

  # --- 3. データ結合 ---
  merged_df <- merge_tcr(pair_data, tra_data, trb_data)

  pair_barcode_col_check <- paste0(PREFIX_TCR_PAIR, "barcode")
  if (nrow(merged_df) == 0 || !pair_barcode_col_check %in% names(merged_df)) {
    warning("Merging TCR data resulted in an empty data frame or missing key barcode column. Returning the empty frame.")
    return(merged_df)
  }

  # --- 4. 代表ID生成と不要列削除 (★修正箇所★) ---
  # TRB鎖のIDを元に代表IDを生成
  trb_clonotype_col <- paste0(PREFIX_TCR_BETA, "raw_clonotype_id")
  trb_subclonotype_col <- paste0(PREFIX_TCR_BETA, "exact_subclonotype_id")
  # raw_consensus_id はここで生成しない

  # 削除対象の元のID列名 (TRA, TRB, Pair) - consensusも含む
  cols_to_remove <- c(
    paste0(PREFIX_TCR_ALPHA, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")),
    paste0(PREFIX_TCR_BETA, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")),
    paste0(PREFIX_TCR_PAIR, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id"))
  )
  # データフレームに実際に存在する列のみを削除対象とする
  cols_to_remove_existing <- intersect(cols_to_remove, names(merged_df))

  processed_df <- merged_df %>%
    dplyr::mutate(
      # raw_clonotype_id は TRB由来のものをそのまま使用
      raw_clonotype_id = if (trb_clonotype_col %in% names(.)) .data[[trb_clonotype_col]] else NA_character_,
      # exact_subclonotype_id は TRB由来の raw_clonotype_id と exact_subclonotype_id を連結
      exact_subclonotype_id = dplyr::case_when(
        # 両方の列が存在し、かつ両方の値が NA でない場合のみ連結
        trb_clonotype_col %in% names(.) & trb_subclonotype_col %in% names(.) &
          !is.na(.data[[trb_clonotype_col]]) & !is.na(.data[[trb_subclonotype_col]])
          ~ stringr::str_c(.data[[trb_clonotype_col]], .data[[trb_subclonotype_col]], sep = "_"),
        # それ以外の場合は NA
        TRUE ~ NA_character_
      )
      # raw_consensus_id はここで生成しない
    ) %>%
    # 不要になった元のTRA/TRB/PairのID列(+ consensus 列)を削除
    dplyr::select(-dplyr::any_of(cols_to_remove_existing)) %>%
    # もしマージや他の処理で raw_consensus_id 列が残ってしまった場合に備えて明示的に削除
    dplyr::select(-dplyr::any_of("raw_consensus_id"))


  # --- 5. 列順序の整理 (★修正箇所★) ---
  # 代表ID列、ペア列、Alpha列、Beta列の順に並べ替え (consensus を除外)
  processed_df <- processed_df %>%
    dplyr::relocate(
      # 代表ID列を先頭に (raw_consensus_id を除外)
      dplyr::any_of(c("raw_clonotype_id", "exact_subclonotype_id")),
      # 次に各プレフィックスを持つ列
      dplyr::starts_with(PREFIX_TCR_PAIR),
      dplyr::starts_with(PREFIX_TCR_ALPHA),
      dplyr::starts_with(PREFIX_TCR_BETA),
      .before = tidyselect::everything() # その他の列の前に配置
    )

  # --- 6. 全長配列の生成 ---
  regions_nt <- c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt")
  regions_aa <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")

  tra_nt_cols <- paste0(PREFIX_TCR_ALPHA, regions_nt)
  trb_nt_cols <- paste0(PREFIX_TCR_BETA, regions_nt)
  tra_aa_cols <- paste0(PREFIX_TCR_ALPHA, regions_aa)
  trb_aa_cols <- paste0(PREFIX_TCR_BETA, regions_aa)

  processed_df <- processed_df %>%
    dplyr::mutate(
      !!paste0(PREFIX_TCR_ALPHA, "full_length_nt") := paste_existing_cols(., tra_nt_cols),
      !!paste0(PREFIX_TCR_BETA, "full_length_nt") := paste_existing_cols(., trb_nt_cols),
      !!paste0(PREFIX_TCR_ALPHA, "full_length_aa") := paste_existing_cols(., tra_aa_cols),
      !!paste0(PREFIX_TCR_BETA, "full_length_aa") := paste_existing_cols(., trb_aa_cols)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::ends_with("full_length_nt") | dplyr::ends_with("full_length_aa"),
        ~ dplyr::na_if(., "")
      )
    )

  # --- 7. 最終的な列整理 (barcode) ---
  pair_barcode_col_rename <- paste0(PREFIX_TCR_PAIR, "barcode")
  if (!"barcode" %in% names(processed_df) && pair_barcode_col_rename %in% names(processed_df)) {
    processed_df <- processed_df %>%
      dplyr::rename(barcode = !!rlang::sym(pair_barcode_col_rename))
  } else if (!"barcode" %in% names(processed_df)) {
      warning(paste("Column", pair_barcode_col_rename, "not found for renaming and 'barcode' column also does not exist."))
  }

  if ("barcode" %in% names(processed_df)) {
    processed_df <- processed_df %>% dplyr::relocate(barcode)
  }

  # --- 8. TRA/TRBペアが存在する行のみフィルタリング ---
  tra_key_col <- paste0(PREFIX_TCR_ALPHA, "chain")
  trb_key_col <- paste0(PREFIX_TCR_BETA, "chain")

  if (tra_key_col %in% names(processed_df) && trb_key_col %in% names(processed_df)) {
      final_df <- processed_df %>%
          dplyr::filter(!is.na(.data[[tra_key_col]]) & !is.na(.data[[trb_key_col]]))

      if(nrow(final_df) == 0 && nrow(processed_df) > 0) {
          warning("Filtering removed all rows. No barcodes with both TRA and TRB chains found after merging.")
      }

  } else {
      missing_cols <- c()
      if (!tra_key_col %in% names(processed_df)) missing_cols <- c(missing_cols, tra_key_col)
      if (!trb_key_col %in% names(processed_df)) missing_cols <- c(missing_cols, trb_key_col)
      warning(paste("Key columns for filtering (", paste(missing_cols, collapse=", "), ") not found in the merged data frame. Skipping the final filtering step."))
      final_df <- processed_df
  }

  # --- 9. 結果を返す ---
  return(final_df)
}
