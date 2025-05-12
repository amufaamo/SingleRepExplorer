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
#' @importFrom scRepertoire combineTCR
#' @importFrom stats setNames
#' @importFrom rlang sym .data

# --- 定数定義 ---
PREFIX_TCR_PAIR  <- "TCR_pair_"
PREFIX_TCR_ALPHA <- "TCR_TRA_" # Alpha鎖 (TRA)
PREFIX_TCR_BETA  <- "TCR_TRB_" # Beta鎖 (TRB)

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


# --- 関数の使用例 ---
# # 解析対象のCSVファイルのパスを指定
# csv_file_tcr <- "path/to/your/tcr_data.csv" # <--- 実際のパスに置き換えてください
#
# # ファイルが存在するか確認してから実行
# if (file.exists(csv_file_tcr)) {
#     # メイン関数を実行してTCRデータを処理
#     processed_tcr_data_paired <- tcr_csv_to_dataframe(csv_file_tcr)
#
#     # 結果の確認
#     print("Processed TCR data (paired chains only):")
#     print(head(processed_tcr_data_paired))
#     print(paste("Dimensions:", paste(dim(processed_tcr_data_paired), collapse=" x ")))
#
#     # 必要であれば、フィルタリング前の行数と比較
#     # merged_df_before_filter <- merge_tcr(csv_to_tcr_pair_dataframe(csv_file_tcr), csv_to_tcr_tra_dataframe(csv_file_tcr), csv_to_tcr_trb_dataframe(csv_file_tcr))
#     # print(paste("Rows before final filter:", nrow(merged_df_before_filter)))
#     # print(paste("Rows after final filter:", nrow(processed_tcr_data_paired)))
#
# } else {
#     warning("Example TCR CSV file not found: ", csv_file_tcr)
# }