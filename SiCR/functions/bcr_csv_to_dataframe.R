# --- ライブラリ読み込み ---
# (必要に応じてコメント解除)
# library(dplyr)
# library(stringr)
# library(tidyselect)
# library(readr)
# library(scRepertoire)
# library(rlang)

# (パッケージとしてビルドする場合の推奨記述)
#' @importFrom dplyr %>% mutate case_when relocate starts_with coalesce across ends_with na_if select all_of any_of filter distinct full_join rename_with rename everything
#' @importFrom stringr str_c str_replace_all
#' @importFrom tidyselect everything
#' @importFrom readr read_csv
#' @importFrom scRepertoire combineBCR
#' @importFrom stats setNames
#' @importFrom rlang sym .data

# --- 定数定義 ---
PREFIX_PAIR <- "BCR_pair_"
PREFIX_IGH <- "BCR_IGH_"
PREFIX_IGL <- "BCR_IGL_" # IGK/IGL共通のプレフィックス

# --- ヘルパー関数定義 ---
# (csv_to_bcr_pair_dataframe, csv_to_bcr_igh_dataframe, csv_to_bcr_igl_dataframe, merge_bcr, paste_existing_cols は変更なし)
# (省略のため、ここでは再掲しません)
# --- (省略されたヘルパー関数) ---
#' @title CSVからBCRペア鎖データを抽出・整形
#' @description (内容は変更なし)
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
    scRepertoire::combineBCR(bcr_raw, samples = sample_name, ID = id_name, removeNA = TRUE, filterMulti = TRUE)
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
    dplyr::rename_with(~ stringr::str_c(PREFIX_PAIR, .), dplyr::everything())
  
  return(pair)
}

#' @title CSVからBCR重鎖(IGH)データを抽出・整形
#' @description (内容は変更なし)
csv_to_bcr_igh_dataframe <- function(csv_path){
  
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  if (is.null(bcr_raw) || !"chain" %in% names(bcr_raw) || !"barcode" %in% names(bcr_raw)) {
    warning("CSV file must contain 'chain' and 'barcode' columns.")
    return(data.frame()) # エラーではなく空を返す
  }
  
  IGH <- bcr_raw %>%
    dplyr::filter(chain == 'IGH') %>%
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    dplyr::rename_with(~ stringr::str_c(PREFIX_IGH, .), dplyr::everything())
  
  return(IGH)
}

#' @title CSVからBCR軽鎖(IGK/IGL)データを抽出・整形
#' @description (内容は変更なし)
csv_to_bcr_igl_dataframe <- function(csv_path){
  
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  if (is.null(bcr_raw) || !"chain" %in% names(bcr_raw) || !"barcode" %in% names(bcr_raw)) {
    warning("CSV file must contain 'chain' and 'barcode' columns.")
    return(data.frame()) # エラーではなく空を返す
  }
  
  IGL <- bcr_raw %>%
    dplyr::filter(chain %in% c("IGL", "IGK")) %>%
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    dplyr::rename_with(~ stringr::str_c(PREFIX_IGL, .), dplyr::everything())
  
  return(IGL)
}

#' @title BCR各鎖データを結合
#' @description (内容は変更なし)
merge_bcr <- function(pair, IGH, IGL){
  
  pair_barcode_col <- paste0(PREFIX_PAIR, "barcode")
  igh_barcode_col <- paste0(PREFIX_IGH, "barcode")
  igl_barcode_col <- paste0(PREFIX_IGL, "barcode")
  
  if (nrow(pair) == 0) {
    warning("Input 'pair' data frame is empty in merge_bcr. Returning empty data frame.")
    return(data.frame())
  }
  if (nrow(IGH) == 0) warning("Input 'IGH' data frame is empty in merge_bcr.")
  if (nrow(IGL) == 0) warning("Input 'IGL' data frame is empty in merge_bcr.")
  
  bcr_merged <- dplyr::full_join(pair, IGH, by = stats::setNames(igh_barcode_col, pair_barcode_col))
  bcr_merged <- dplyr::full_join(bcr_merged, IGL, by = stats::setNames(igl_barcode_col, pair_barcode_col))
  
  return(bcr_merged)
}

#' @title 存在する列の値を連結する内部ヘルパー関数
#' @description (内容は変更なし)
paste_existing_cols <- function(df, cols) {
  existing_cols <- intersect(cols, names(df))
  if (length(existing_cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  apply(df[, existing_cols, drop = FALSE], 1, function(row_values) {
    row_values_no_na <- ifelse(is.na(row_values), "", row_values)
    paste0(row_values_no_na, collapse = "")
  })
}
# --- メイン関数定義 (修正版) ---

#' @title CSVからBCRデータを処理するメイン関数 (ペア必須・ID修正版)
#' @description 指定されたCSVファイルからBCRデータを読み込み、ペア鎖、重鎖、軽鎖に
#'              分割・整形した後、結合し、IGH鎖のIDを元に代表ID
#'              (`raw_clonotype_id`, `exact_subclonotype_id`) を格納し、
#'              全長配列を生成し、列順序を整え、barcode列名を変更します。
#'              `raw_consensus_id` は削除されます。
#'              `exact_subclonotype_id` は IGH の `raw_clonotype_id` と
#'              `exact_subclonotype_id` を `_` で連結した値になります。
#'              **最終的に重鎖 (IGH) と軽鎖 (IGK/IGL) の両方が存在する行のみを保持します。**
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 全ての処理が完了し、IGH/IGK+IGLペアが存在する行のみを含むBCRデータフレーム。
#'         処理中にエラーが発生した場合や、有効なペアデータが見つからない場合は、
#'         空のデータフレームまたは警告付きの部分的な結果を返すことがあります。
bcr_csv_to_dataframe <- function(csv_path){

  # --- 1. 入力チェック ---
  if (!file.exists(csv_path)) {
    stop("Input CSV file not found: ", csv_path)
    return(NULL)
  }

  # --- 2. データ抽出・整形 ---
  pair_data <- csv_to_bcr_pair_dataframe(csv_path)
  igh_data <- csv_to_bcr_igh_dataframe(csv_path)
  igl_data <- csv_to_bcr_igl_dataframe(csv_path)

  if (is.null(pair_data) || nrow(pair_data) == 0 || !paste0(PREFIX_PAIR, "barcode") %in% names(pair_data)) {
    warning("No valid paired BCR data found or extracted by combineBCR. Returning an empty data frame.")
    return(data.frame())
  }
  if (is.null(igh_data) || nrow(igh_data) == 0) {
    warning("No IGH data extracted. Final result might be empty after filtering.")
    if (is.null(igh_data)) igh_data <- data.frame()
  }
  if (is.null(igl_data) || nrow(igl_data) == 0) {
    warning("No IGK/IGL data extracted. Final result might be empty after filtering.")
    if (is.null(igl_data)) igl_data <- data.frame()
  }

  # --- 3. データ結合 ---
  merged_df <- merge_bcr(pair_data, igh_data, igl_data)

  pair_barcode_col_check <- paste0(PREFIX_PAIR, "barcode")
  if (nrow(merged_df) == 0 || !pair_barcode_col_check %in% names(merged_df)) {
    warning("Merging BCR data resulted in an empty data frame or missing key barcode column. Returning the empty frame.")
    return(merged_df)
  }

  # --- 4. 代表ID生成と不要列削除 (★修正箇所★) ---
  # 代表IDとしてIGH鎖のIDを使用
  igh_clonotype_col <- paste0(PREFIX_IGH, "raw_clonotype_id")
  # igh_consensus_col <- paste0(PREFIX_IGH, "raw_consensus_id") # consensus は使用しない
  igh_subclonotype_col <- paste0(PREFIX_IGH, "exact_subclonotype_id")

  # 削除対象の元のID列名 (IGH, IGL, Pair) - consensusも含む
  cols_to_remove <- c(
    paste0(PREFIX_IGH, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")),
    paste0(PREFIX_IGL, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")),
    paste0(PREFIX_PAIR, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id"))
  )
  cols_to_remove_existing <- intersect(cols_to_remove, names(merged_df))

  processed_df <- merged_df %>%
    dplyr::mutate(
      # raw_clonotype_id は IGH 由来のものをそのまま使用
      raw_clonotype_id = if (igh_clonotype_col %in% names(.)) .data[[igh_clonotype_col]] else NA_character_,
      # exact_subclonotype_id は IGH 由来の raw_clonotype_id と exact_subclonotype_id を連結
      exact_subclonotype_id = dplyr::case_when(
        # 両方の列が存在し、かつ両方の値が NA でない場合のみ連結
        igh_clonotype_col %in% names(.) & igh_subclonotype_col %in% names(.) &
          !is.na(.data[[igh_clonotype_col]]) & !is.na(.data[[igh_subclonotype_col]])
          ~ stringr::str_c(.data[[igh_clonotype_col]], .data[[igh_subclonotype_col]], sep = "_"),
        # それ以外の場合は NA
        TRUE ~ NA_character_
      )
      # raw_consensus_id は生成しない
    ) %>%
    # 不要になった元のID列(+ consensus 列)を削除
    dplyr::select(-dplyr::any_of(cols_to_remove_existing)) %>%
    # もしマージや他の処理で raw_consensus_id 列が残ってしまった場合に備えて明示的に削除
    dplyr::select(-dplyr::any_of("raw_consensus_id"))

  # --- 5. 列順序の整理 (★修正箇所★) ---
  processed_df <- processed_df %>%
    dplyr::relocate(
      # 代表ID列 (consensus を除外)
      dplyr::any_of(c("raw_clonotype_id", "exact_subclonotype_id")),
      # ペア鎖、重鎖、軽鎖の順
      dplyr::starts_with(PREFIX_PAIR),
      dplyr::starts_with(PREFIX_IGH),
      dplyr::starts_with(PREFIX_IGL),
      .before = tidyselect::everything()
    )

  # --- 6. 全長配列の生成 ---
  regions_nt <- c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt")
  regions_aa <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")

  igh_nt_cols <- paste0(PREFIX_IGH, regions_nt)
  igl_nt_cols <- paste0(PREFIX_IGL, regions_nt)
  igh_aa_cols <- paste0(PREFIX_IGH, regions_aa)
  igl_aa_cols <- paste0(PREFIX_IGL, regions_aa)

  processed_df <- processed_df %>%
    dplyr::mutate(
      !!paste0(PREFIX_IGH, "full_length_nt") := paste_existing_cols(., igh_nt_cols),
      !!paste0(PREFIX_IGL, "full_length_nt") := paste_existing_cols(., igl_nt_cols),
      !!paste0(PREFIX_IGH, "full_length_aa") := paste_existing_cols(., igh_aa_cols),
      !!paste0(PREFIX_IGL, "full_length_aa") := paste_existing_cols(., igl_aa_cols)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::ends_with("full_length_nt") | dplyr::ends_with("full_length_aa"),
        ~ dplyr::na_if(., "")
      )
    )

  # --- 7. 最終的な列整理 (barcode) ---
  pair_barcode_col <- paste0(PREFIX_PAIR, "barcode")
  if (!"barcode" %in% names(processed_df) && pair_barcode_col %in% names(processed_df)) {
    processed_df <- processed_df %>%
      dplyr::rename(barcode = !!rlang::sym(pair_barcode_col))
  } else if (!"barcode" %in% names(processed_df)) {
    warning(paste("Column", pair_barcode_col, "not found and 'barcode' column also does not exist."))
  }

  if ("barcode" %in% names(processed_df)) {
    processed_df <- processed_df %>% dplyr::relocate(barcode)
  }

  # --- 8. IGH/IGLペアが存在する行のみフィルタリング ---
  igh_key_col <- paste0(PREFIX_IGH, "chain")
  igl_key_col <- paste0(PREFIX_IGL, "chain")

  if (igh_key_col %in% names(processed_df) && igl_key_col %in% names(processed_df)) {
    final_df <- processed_df %>%
      dplyr::filter(!is.na(.data[[igh_key_col]]) & !is.na(.data[[igl_key_col]]))

    if(nrow(final_df) == 0 && nrow(processed_df) > 0) {
      warning("Filtering removed all rows. No barcodes with both IGH and IGK/IGL chains found after merging.")
    }

  } else {
    missing_cols <- c()
    if (!igh_key_col %in% names(processed_df)) missing_cols <- c(missing_cols, igh_key_col)
    if (!igl_key_col %in% names(processed_df)) missing_cols <- c(missing_cols, igl_key_col)
    warning(paste("Key columns for filtering (", paste(missing_cols, collapse=", "), ") not found in the merged data frame. Skipping the final filtering step."))
    final_df <- processed_df
  }


  # --- 9. 結果を返す ---
  return(final_df)
}

# --- 関数の使用例 ---
# # 解析対象のCSVファイルのパスを指定
# csv_file_bcr <- "path/to/your/bcr_data.csv" # <--- 実際のパスに置き換えてください
#
# # ファイルが存在するか確認してから実行
# if (file.exists(csv_file_bcr)) {
#     # メイン関数を実行してBCRデータを処理
#     processed_bcr_data_paired <- bcr_csv_to_dataframe(csv_file_bcr)
#
#     # 結果の確認
#     print("Processed BCR data (paired chains only):")
#     print(head(processed_bcr_data_paired))
#     print(paste("Dimensions:", paste(dim(processed_bcr_data_paired), collapse=" x ")))
#
# } else {
#     warning("Example BCR CSV file not found: ", csv_file_bcr)
# }