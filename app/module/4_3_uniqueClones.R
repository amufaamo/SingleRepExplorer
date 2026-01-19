#source("../utils.R")
source("utils.R")
# UI 関数 (変更なしと仮定)
uniqueClonesUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      # データタイプ選択 (vdjType が input$vdj_type を生成すると仮定)
      vdjType(ns),
      # グループ化変数選択 (groupByInput が input$group_by を生成すると仮定)
      groupByInput(ns),
      valueType(ns),
      # プロットタイプ選択 (必要であれば追加)
      commonPlotOptions(ns),
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      h3("Table"),
      downloadButton(ns("download_table"), "Downloadtable (.csv)"),
      DTOutput(ns("table"))

    ),
  )
}


# Server 関数修正版
uniqueClonesServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    reactive_df_raw <- reactive({
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df <- myReactives$bcr_df
      }
      return(df)
    })

     summary_data <- reactive({
      df <- reactive_df_raw()
      # req チェックはそのまま
      req(!is.null(df), input$group_by)
      req(is.character(input$group_by), nzchar(input$group_by))
      req(input$group_by %in% names(df))
      req("raw_clonotype_id" %in% names(df))

      grouping_column_name <- input$group_by
      message("Grouping column name: ", grouping_column_name) # デバッグ用

      # summarise までの処理
      summary_result <- df %>%
        group_by(!!sym(grouping_column_name)) %>%
        summarise(
          unique_clonotype_count = n_distinct(raw_clonotype_id),
          total_clonotype_count = n(),
          proportion_unique = ifelse(total_clonotype_count > 0,
                                     unique_clonotype_count / total_clonotype_count,
                                     NA_real_),
          .groups = "drop"
        )
        # ungroup() は .groups = "drop" があれば通常不要

      # ★ rename の代替: names() を使用 ★
      message("Column names before rename: ", paste(names(summary_result), collapse=", ")) # 変更前の列名確認

      # grouping_column_name と一致する列名を "Group" に変更
      # which() を使うとより安全かもしれません
      col_index <- which(names(summary_result) == grouping_column_name)
      if (length(col_index) == 1) { # 確実に1つだけ見つかった場合に変更
          names(summary_result)[col_index] <- "Group"
          message("Column names after rename: ", paste(names(summary_result), collapse=", ")) # 変更後の列名確認
      } else {
          warning("Could not find unique column to rename: ", grouping_column_name)
          # エラー処理: 見つからない場合や複数見つかる場合どうするか？
          # ここでは何もしない（元の列名のまま）か、エラーを発生させる等
      }

      # デバッグ用に結果を表示
      print(head(summary_result))

      return(summary_result)
    })

    output$table <- renderDT({
      summary_data()
    })
    plot <- reactive({
      plot_df <- summary_data()
      req(plot_df)
      req(all(c("Group", "unique_clonotype_count", "proportion_unique") %in% names(plot_df)))
      req(input$value_type)

      p_base <- ggplot(plot_df, aes(x = Group, fill = Group))

      # y軸のマッピングを決定
      y_mapping <- if (input$value_type == "count") {
                     aes(y = unique_clonotype_count)
                   } else if (input$value_type == "percentage") {
                     aes(y = proportion_unique)
                   } else {
                     # エラーハンドリング
                     warning("Invalid input$value_type for plot: ", input$value_type)
                     return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Invalid 'value_type'") + theme_void())
                   }

      # y軸のラベルを決定
      y_label <- if (input$value_type == "count") {
                   "Unique Clonotype Count"
                 } else { # percentage の場合
                   "Proportion of Unique Clonotypes (%)"
                 }

      # y軸のスケールラベル関数を決定
      y_scale_labels <- if (input$value_type == "count") {
                          waiver() # デフォルトのラベルを使用
                        } else { # percentage の場合
                          scales::percent_format(accuracy = 1)
                        }

      # プロットを構築
      p <- p_base +
        y_mapping + # y軸マッピングを追加
        geom_bar(stat = "identity") +
        # --- scale_y_continuous をここで一元管理 ---
        scale_y_continuous(
          labels = y_scale_labels, # value_type に応じたラベル関数
          # expand = expansion(mult = c(0, 0.05)) # 下0、上5%余裕を持たせる場合
          # または常に expand = c(0, 0) にする場合：
          expand = c(0, 0)
        ) +
        # -----------------------------------------
        labs(y = y_label) + # y軸ラベルを設定
        theme_classic() +
        theme(legend.position = input$legend)

      return(p)
    })

    
    output$plot <- renderPlot(
      {
        plot()
      },
      width = reactive(input$plot_width),
      height = reactive(input$plot_height),
    )

    output$download_plot <- downloadHandler(
      filename = function() {
        "plot.pdf"
      }, # ファイル名を修正
      content = function(file) {
        ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300) # renderCustomPlot() を plot() に修正
      }
    )

  output$download_table <- downloadHandler(
      filename = function() {
        "table.csv"
      },
      content = function(file) {
        write.csv(summary_data(), file, row.names = FALSE)
      }
    )

  })
}
