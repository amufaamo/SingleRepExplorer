# UI部分
clonotypeInformationUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      vdjType(ns),
    ),
    mainPanel(
      DTOutput(ns("table")),
      downloadButton(ns("download_table"), "Download table (.csv)") # 追加
    )
  )
}

# Server部分

clonotypeInformationServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    table <- reactive({
      # --- DEBUG START (4_1_clonotypeInformation.R) ---
      print(paste0("DEBUG: clonotypeInformationServer - vdj_type: ", input$vdj_type))
      # browser() # ここで実行を一時停止して myReactives$tcr_df の内容を確認できます
      # --- DEBUG END ---
      if (input$vdj_type == "tcr") {
        # TCRデータフレームが存在することを確認
        req(myReactives$tcr_df)

        tcr_result_df <- myReactives$tcr_df %>%
          dplyr::select(
            dplyr::any_of(c(
            "raw_clonotype_id",
            "exact_subclonotype_id",
            "TCR_pair_CTgene",
            "TCR_pair_CTnt",
            "TCR_pair_CTaa",
            "TCR_pair_CTstrict",
            "TCR_TRA_v_gene",
            "TCR_TRA_j_gene",
            "TCR_TRA_c_gene",
            "TCR_TRA_fwr1",
            "TCR_TRA_fwr1_nt",
            "TCR_TRA_cdr1",
            "TCR_TRA_cdr1_nt",
            "TCR_TRA_fwr2",
            "TCR_TRA_fwr2_nt",
            "TCR_TRA_cdr2",
            "TCR_TRA_cdr2_nt",
            "TCR_TRA_fwr3",
            "TCR_TRA_fwr3_nt",
            "TCR_TRA_cdr3",
            "TCR_TRA_cdr3_nt",
            "TCR_TRA_fwr4",
            "TCR_TRA_fwr4_nt",
            "TCR_TRB_v_gene",
            "TCR_TRB_d_gene",
            "TCR_TRB_j_gene",
            "TCR_TRB_c_gene",
            "TCR_TRB_fwr1",
            "TCR_TRB_fwr1_nt",
            "TCR_TRB_cdr1",
            "TCR_TRB_cdr1_nt",
            "TCR_TRB_fwr2",
            "TCR_TRB_fwr2_nt",
            "TCR_TRB_cdr2",
            "TCR_TRB_cdr2_nt",
            "TCR_TRB_fwr3",
            "TCR_TRB_fwr3_nt",
            "TCR_TRB_cdr3",
            "TCR_TRB_cdr3_nt",
            "TCR_TRB_fwr4",
            "TCR_TRB_fwr4_nt",
            "TCR_TRA_full_length_nt",
            "TCR_TRB_full_length_nt",
            "TCR_TRA_full_length_aa",
            "TCR_TRB_full_length_aa"
          ))
           ) %>%
          distinct() %>%
          # clonotype_id が NA でない行のみを対象にする
          filter(!is.na(raw_clonotype_id)) %>%
          mutate(
            # parse_number を使って数値部分を抽出し、一時的な列 clonotype_num を作成
            clonotype_num = readr::parse_number(raw_clonotype_id),
            # subclonotype_id が NA の場合の処理を追加
            exact_subclonotype_id = tidyr::replace_na(exact_subclonotype_id, "")
          ) %>%
          # ID と subclonotype ID に基づいて並び替え
          arrange(
            clonotype_num, # まず clonotype 番号 (昇順)
            exact_subclonotype_id # 次に exact_subclonotype_id (昇順)
          ) %>%
          # 並び替えに使用した一時的な列を削除
          dplyr::select(-clonotype_num)
        # --- DEBUG START (4_1_clonotypeInformation.R) ---
        print(paste0("DEBUG: TCR table() reactive - final nrow: ", nrow(tcr_result_df) ))
        # browser() # ここで実行を一時停止して最終的なテーブルの内容を確認できます
        # --- DEBUG END ---
        return(tcr_result_df) # 処理されたデータフレームを返す
      } else if (input$vdj_type == "bcr") {
        # myReactives$bcr_df が存在するか確認 (任意だが推奨)
        req(myReactives$bcr_df)

        bcr_result_df <- myReactives$bcr_df %>% # BCR側も同様に変数に格納
          # any_of を使用して、万が一列が存在しなくてもエラーにならないようにする
          dplyr::select(
            dplyr::any_of(c(
              # 代表ID
              "raw_clonotype_id",
              "raw_consensus_id",
              "exact_subclonotype_id",
              # ペア情報
              "BCR_pair_CTgene",
              "BCR_pair_CTnt",
              "BCR_pair_CTaa",
              "BCR_pair_CTstrict",
              # 重鎖(IGH)情報
              "BCR_IGH_v_gene",
              "BCR_IGH_d_gene",
              "BCR_IGH_j_gene",
              "BCR_IGH_c_gene",
              "BCR_IGH_fwr1", "BCR_IGH_fwr1_nt",
              "BCR_IGH_cdr1", "BCR_IGH_cdr1_nt",
              "BCR_IGH_fwr2", "BCR_IGH_fwr2_nt",
              "BCR_IGH_cdr2", "BCR_IGH_cdr2_nt",
              "BCR_IGH_fwr3", "BCR_IGH_fwr3_nt",
              "BCR_IGH_cdr3", "BCR_IGH_cdr3_nt",
              "BCR_IGH_fwr4", "BCR_IGH_fwr4_nt",
              # 軽鎖(IGL/IGK)情報 (プレフィックスは BCR_IGL_ で統一されている想定)
              "BCR_IGL_v_gene",
              "BCR_IGL_j_gene",
              "BCR_IGL_c_gene",
              "BCR_IGL_fwr1", "BCR_IGL_fwr1_nt",
              "BCR_IGL_cdr1", "BCR_IGL_cdr1_nt",
              "BCR_IGL_fwr2", "BCR_IGL_fwr2_nt",
              "BCR_IGL_cdr2", "BCR_IGL_cdr2_nt",
              "BCR_IGL_fwr3", "BCR_IGL_fwr3_nt",
              "BCR_IGL_cdr3", "BCR_IGL_cdr3_nt",
              "BCR_IGL_fwr4", "BCR_IGL_fwr4_nt",
              # 全長配列
              "BCR_IGH_full_length_nt",
              "BCR_IGL_full_length_nt",
              "BCR_IGH_full_length_aa",
              "BCR_IGL_full_length_aa"
            ))
          ) %>%
          distinct() %>%
          # clonotype_id が NA でない行のみを対象にする
          filter(!is.na(raw_clonotype_id)) %>%
          mutate(
            clonotype_num = readr::parse_number(raw_clonotype_id)
          ) %>%
          # subclonotype_id が NA の場合の処理を追加
          mutate(exact_subclonotype_id = tidyr::replace_na(exact_subclonotype_id, "")) %>%
          arrange(
            clonotype_num,
            exact_subclonotype_id
          ) %>%
          dplyr::select(-clonotype_num)
        return(bcr_result_df) # 処理されたデータフレームを返す
      } else {
        # vdj_type が "tcr" でも "bcr" でもない場合の処理 (例: NULLを返す)
        NULL
      }
    })

    # req(input$vdj_type)

    # table <- reactive({
    #   if (vdj == 'tcr'){
    #     df <- myReactives$tcr_df$S1
    #   } else if (vdj == 'bcr'){
    #     df <- myReactives$bcr_df$S1
    #   }
    #   req(df)
    #   df <- df %>%
    #     mutate(
    #       exact_subclonotype_id = paste0(raw_clonotype_id, "_", exact_subclonotype_id),
    #       clonotype_num = as.numeric(gsub("clonotype(\\d+)_.*", "\\1", exact_subclonotype_id)), # clonotype の後の数字
    #       subclonotype_num = as.numeric(gsub(".*_(\\d+)", "\\1", exact_subclonotype_id))      # _ の後の数字
    #     ) %>%
    #     arrange(clonotype_num, subclonotype_num) %>% # 数値でソート
    #     select(raw_clonotype_id, exact_subclonotype_id, CTstrict, CTaa, CTnt) %>%
    #     distinct()
    #   return(df)
    # })

    # DTテーブルの表示
    output$table <- renderDT({
      datatable(
        table(),
        filter = "top",
        options = list(pageLength = 20)
      )
    })

    # CSVダウンロードハンドラー
    output$download_table <- downloadHandler(
      filename = function() {
        "table.csv"
      },
      content = function(file) {
        write.csv(table(), file) # table()で生成したデータをCSVとして保存
      }
    )
  })
}
