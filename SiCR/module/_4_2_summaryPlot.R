# --- UI ---
summaryPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      # データタイプ選択
      vdjType(ns),
      # グループ化変数選択
      groupByInput(ns),
      # プロットタイプ選択
      selectInput(ns("plot_type"), "Plot Type", choices = c("Normal" = "normal", "Stacked" = "stacked", "Dodged" = "dodged"), selected = "dodged"),
      # 表示値選択 (Count or Percentage)
      selectInput(ns("display_value"), "Display Value",
        choices = c("Count" = "count", "Percentage" = "percentage"),
        selected = "count"
      ),
      # 表示する上位N件
      numericInput(ns("top_n"), "Top N Clonotypes:", value = 10, min = 1, step = 1),
      # グループ要素フィルタリング用UI
      uiOutput(ns("filter_groups_ui")),
      conditionalPanel(
        condition = "input.plot_type == 'dodged'",
        ns = ns,
        selectInput(ns("sort_by_group"), "Sort X-axis by group:", choices = NULL) # 選択肢は動的に更新
      ),
      # 横向きグラフにするか
      checkboxInput(ns("horizontal"), "Horizontal Plot", value = FALSE),
      commonPlotOptions(ns),
    ),
    mainPanel(
      # プロット出力
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      h3("Table"),
      downloadButton(ns("download_table"), "Downloadtable (.csv)"),
      DTOutput(ns("table"))
    )
  )
}


# --- サーバー ---
summaryPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    # Seuratオブジェクトが変更されたらGroup byの選択肢を更新
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })
    # 1. 元データのリアクティブ (NAと空文字列を除去)
    reactive_df_raw <- reactive({
      req(input$vdj_type, input$group_by) # vdj_typeとgroup_byが選択されていることを要求
      df <- NULL
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df <- myReactives$bcr_df
      } else {
        warning(paste("Unexpected vdj_type:", input$vdj_type))
        return(NULL)
      }
      # 必要な列が存在するか確認
      req(df, "raw_clonotype_id" %in% names(df), nrow(df) > 0)
      validate(need(input$group_by %in% names(df), paste("Column", shQuote(input$group_by), "not found in the data.")))

      # group_by 列の NA と空文字列を除外
      df_filtered_na_empty <- df %>%
        dplyr::filter(!is.na(.data[[input$group_by]]) & .data[[input$group_by]] != "")
      # 必要であれば他の不要な値も除外 (例: スペースのみ)
      # dplyr::filter(!grepl("^\\s*$", .data[[input$group_by]]))

      # 除外後にデータが残っているか確認
      if (nrow(df_filtered_na_empty) == 0) {
        showNotification(paste("No data remaining after removing NA/empty values from column:", input$group_by), type = "warning", duration = 5)
        return(NULL)
      }

      return(df_filtered_na_empty)
    })

    # 2. フィルタリング前の合計数を計算 (group_col を character に)
    original_totals <- reactive({
      df <- reactive_df_raw()
      req(df, input$group_by)
      grouping_var <- input$group_by

      # グループごとの合計数
      group_tots <- df %>%
        dplyr::count(.data[[grouping_var]], name = "original_group_total") %>%
        dplyr::rename(group_col = .data[[grouping_var]]) %>%
        # ★ 結合のため character 型に統一
        dplyr::mutate(group_col = as.character(group_col))

      # 全体の合計数
      overall_tot <- nrow(df)

      # グループ合計が計算できたか確認 (0行でもOKとする)
      # req(nrow(group_tots) > 0)

      return(list(group = group_tots, overall = overall_tot))
    })

    # 3. グループフィルタリングUIの生成
    output$filter_groups_ui <- renderUI({
      df <- reactive_df_raw()
      req(df, input$group_by) # dfがNULLでないこと、group_byが指定されていること
      grouping_var <- input$group_by

      # NA/空文字は除外済み
      available_groups <- sort(unique(df[[grouping_var]]))

      # 利用可能なグループがない場合はUIを表示しない
      validate(
        need(length(available_groups) > 0, paste("No available groups found in column:", shQuote(grouping_var), "after removing NA/empty values."))
      )

      checkboxGroupInput(session$ns("filter_groups"),
        label = paste("Filter", tools::toTitleCase(gsub("_", " ", grouping_var)), ":"),
        choices = available_groups,
        selected = available_groups, # 初期状態は全て選択
        inline = TRUE
      )
    })

    # 4. フィルタリングされたデータフレーム
    filtered_df <- reactive({
      df <- reactive_df_raw()
      # UIが生成され、選択が存在することを要求
      req(df, input$group_by, input$filter_groups)

      grouping_var <- input$group_by
      selected_groups <- input$filter_groups

      # 選択されたグループでフィルタリング (NA/空文字はないはず)
      df_filtered <- df %>%
        dplyr::filter(.data[[grouping_var]] %in% selected_groups)

      # フィルタリング後にデータが残っているか確認
      if (nrow(df_filtered) == 0) {
        showNotification("No data remaining after filtering by selected groups.", type = "warning", duration = 5)
        return(NULL)
      }
      return(df_filtered)
    })

    # ★ 追加: 並び替え基準グループの選択肢を更新
    observe({
      # Dodgedプロットが選ばれ、フィルターグループが存在する場合のみ更新
      req(input$plot_type == "dodged", input$filter_groups)
      # 利用可能な選択肢（フィルターされたグループ + Overall）
      choices_with_overall <- c("Overall (All selected groups)" = "overall", input$filter_groups)
      # 現在の選択を維持しようと試みる、なければ"overall"を選択
      current_selection <- isolate(input$sort_by_group)
      selected_value <- if (!is.null(current_selection) && current_selection %in% choices_with_overall) {
        current_selection
      } else {
        "overall"
      }
      updateSelectInput(session, "sort_by_group",
        choices = choices_with_overall,
        selected = selected_value
      )
    })


    # 5. プロット/テーブル用データの計算 (★修正: Dodgedの並び替え)
    table_for_plot <- reactive({
      df_filtered <- filtered_df()
      req(df_filtered, input$plot_type, input$top_n, input$top_n >= 1, input$group_by, input$filter_groups)
      og_totals <- original_totals()
      req(og_totals)

      plot_data <- NULL
      n_top <- as.integer(input$top_n)
      grouping_var <- input$group_by


      # --- Normal Plot ---
      if (input$plot_type == "normal") {
        plot_data_raw <- df_filtered %>%
          dplyr::count(raw_clonotype_id, sort = TRUE, name = "n")

        total_cells_original <- og_totals$overall
        req(total_cells_original) # 元の合計数が存在すること

        plot_data <- plot_data_raw %>%
          dplyr::slice_head(n = n_top) %>%
          dplyr::mutate(percentage = ifelse(total_cells_original == 0, 0, n / total_cells_original)) %>%
          # raw_clonotype_id を factor 化 (プロットの順番のため)
          dplyr::mutate(raw_clonotype_id = factor(raw_clonotype_id, levels = .$raw_clonotype_id))

        # --- Stacked/Dodged Plot (★修正箇所あり) ---
      } else if (input$plot_type %in% c("stacked", "dodged")) {
        # フィルター後のデータでカウント
        counts_by_group_raw <- df_filtered %>%
          dplyr::count(raw_clonotype_id, .data[[grouping_var]], name = "n") %>%
          dplyr::rename(group_col = .data[[grouping_var]])

        # --- ★ 並び替え基準の決定 (Dodgedの場合のみinputを使用) ---
        clonotype_order_info <- NULL
        if (input$plot_type == "dodged") {
          # Dodgedの場合、並び替え基準の入力を要求
          req(input$sort_by_group)
          sort_criterion_group <- input$sort_by_group
        } else {
          # Stackedの場合は常にOverallでソート
          sort_criterion_group <- "overall"
        }

        # 全体のクロノタイプリスト（基準グループにない場合も考慮するため）
        all_clonotypes <- distinct(counts_by_group_raw, raw_clonotype_id)

        if (sort_criterion_group == "overall") {
          # Overall: フィルターされた全グループでの合計頻度でソート
          clonotype_order_info <- counts_by_group_raw %>%
            dplyr::group_by(raw_clonotype_id) %>%
            dplyr::summarise(sort_value = sum(n), .groups = "drop") %>%
            dplyr::right_join(all_clonotypes, by = "raw_clonotype_id") %>% # 全クロノタイプを保持
            tidyr::replace_na(list(sort_value = 0)) %>% # 合計がないものは0
            dplyr::arrange(dplyr::desc(sort_value)) # 降順ソート
        } else {
          # 特定グループ: そのグループでの頻度でソート
          clonotype_order_info <- counts_by_group_raw %>%
            dplyr::filter(group_col == sort_criterion_group) %>% # 基準グループでフィルタ
            dplyr::select(raw_clonotype_id, sort_value = n) %>% # nをsort_valueとする
            dplyr::right_join(all_clonotypes, by = "raw_clonotype_id") %>% # 全クロノタイプを保持
            tidyr::replace_na(list(sort_value = 0)) %>% # 基準グループにないものは0
            dplyr::arrange(dplyr::desc(sort_value)) # 降順ソート
        }
        # --- ★ 並び替え基準決定 終了 ---

        validate(need(!is.null(clonotype_order_info) && nrow(clonotype_order_info) > 0, "Could not determine clonotype order."))

        # 上位N件を決定 (新しい並び順で)
        top_clonotypes_sorted <- clonotype_order_info %>%
          dplyr::slice_head(n = n_top)

        validate(need(nrow(top_clonotypes_sorted) > 0, "No top clonotypes found based on the sorting criterion."))

        # ★ 上位N件のIDリスト (新しい並び順)
        top_clonotype_ids_new_order <- top_clonotypes_sorted$raw_clonotype_id

        # 実際にデータに存在する有効なグループ（フィルター後）
        valid_group_levels <- intersect(input$filter_groups, unique(counts_by_group_raw$group_col))
        validate(need(length(valid_group_levels) > 0, "No valid groups found after filtering."))

        # 上位Nと有効なグループでフィルタリングし、group_colをFactor化
        counts_by_group_filtered <- counts_by_group_raw %>%
          # ★ 新しい上位Nリストでフィルタ
          dplyr::filter(raw_clonotype_id %in% top_clonotype_ids_new_order) %>%
          dplyr::mutate(group_col = factor(as.character(group_col), levels = valid_group_levels))

        # complete で使用するレベル
        group_levels_to_complete <- factor(valid_group_levels, levels = valid_group_levels)
        # ★ 新しい並び順の Factor レベル
        clonotype_levels_to_complete <- factor(top_clonotype_ids_new_order, levels = top_clonotype_ids_new_order)

        # 全組み合わせを生成、結合前にgroup_colをcharacterに戻す
        plot_data_completed <- counts_by_group_filtered %>%
          tidyr::complete(
            raw_clonotype_id = clonotype_levels_to_complete,
            group_col = group_levels_to_complete, fill = list(n = 0)
          ) %>%
          tidyr::drop_na(raw_clonotype_id) %>%
          dplyr::mutate(group_col = as.character(group_col)) # ★結合のため

        # 元のグループ合計数を取得
        original_group_totals_df <- og_totals$group
        req(original_group_totals_df)

        # 元の合計数と結合しPercentage計算
        plot_data_with_perc <- plot_data_completed %>%
          dplyr::left_join(original_group_totals_df, by = "group_col") %>%
          dplyr::mutate(percentage = ifelse(is.na(original_group_total) | original_group_total == 0, 0, n / original_group_total)) %>%
          dplyr::select(-original_group_total)

        validate(need(nrow(plot_data_with_perc) > 0, "Data processing empty."))

        # フィルター後の合計頻度(total_n_filtered)を結合し、プロット用にFactor化
        plot_data <- plot_data_with_perc %>%
          dplyr::left_join(
            # 参考用にフィルター後の合計頻度を計算
            counts_by_group_raw %>%
              dplyr::group_by(raw_clonotype_id) %>%
              dplyr::summarise(total_n_filtered = sum(n), .groups = "drop"),
            by = "raw_clonotype_id"
          ) %>%
          dplyr::mutate(
            # ★ 新しい並び順で Factor 化
            raw_clonotype_id = factor(raw_clonotype_id, levels = top_clonotype_ids_new_order),
            group_col = factor(group_col, levels = valid_group_levels) # group_col はそのまま
          ) %>%
          dplyr::mutate(total_n_filtered = coalesce(total_n_filtered, 0)) # NAを0に
      } else { # plot_type が予期せぬ値の場合
        warning(paste("Unexpected plot_type:", input$plot_type))
        return(NULL)
      }

      validate(need(!is.null(plot_data) && nrow(plot_data) > 0, "No data for plot/table."))
      return(plot_data)
    })

    # 6. プロットオブジェクトの生成
    plot_obj <- reactive({
      plot_data <- table_for_plot()
      req(plot_data) # plot_data が NULL でないことを要求
      grouping_var <- input$group_by

      # 表示値に応じて y軸変数、ラベル、スケールを設定
      if (input$display_value == "count") {
        y_var <- "n"
        y_label <- "Frequency (Number of Cells)"
        y_scale <- ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
      } else { # percentage
        y_var <- "percentage"
        if (input$plot_type == "normal") {
          y_label <- "Frequency (% of Original Total)"
        } else {
          y_label <- paste("Frequency (% within Original", tools::toTitleCase(gsub("_", " ", grouping_var)), ")")
        }
        y_scale <- ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.05)))
      }

      # 基本プロット
      p <- ggplot2::ggplot(plot_data)
      fill_label <- if (input$plot_type %in% c("stacked", "dodged")) {
        tools::toTitleCase(gsub("_", " ", grouping_var))
      } else {
        NULL
      }

      # プロットタイプ別の geom と aes
      if (input$plot_type == "normal") {
        p <- p + ggplot2::aes(x = raw_clonotype_id, y = .data[[y_var]]) +
          ggplot2::geom_bar(stat = "identity", fill = "steelblue")
      } else if (input$plot_type %in% c("stacked", "dodged")) {
        # group_col は既に Factor になっているはず
        valid_levels <- levels(plot_data$group_col)
        validate(need(length(valid_levels) > 0, "No valid group levels for plotting."))

        p <- p + ggplot2::aes(x = raw_clonotype_id, y = .data[[y_var]], fill = group_col)

        if (input$plot_type == "stacked") {
          p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")
        } else { # dodged
          p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(preserve = "single"))
        }
      }

      # プロットタイトル（フィルター状況を反映）
      plot_title <- paste("Top", input$top_n, "Clonotype Frequencies")
      all_groups_original <- tryCatch(unique(reactive_df_raw()[[input$group_by]]), error = function(e) NULL)
      if (!is.null(all_groups_original) && length(input$filter_groups) < length(all_groups_original)) {
        plot_title <- paste(plot_title, "(Groups Filtered)")
      }

      # テーマ、ラベル、スケールの適用
      p <- p +
        ggplot2::scale_x_discrete(limits = levels(plot_data$raw_clonotype_id)) + # x軸の順番を制御
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          axis.title = ggplot2::element_text(size = 12),
          plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
        ) +
        theme(legend.position = input$legend) +
        y_scale + # y軸スケール適用
        ggplot2::labs(
          title = plot_title,
          x = "Raw Clonotype ID",
          y = y_label,
          fill = fill_label
        )

      # 横向きプロットの適用
      if (input$horizontal) {
        p <- p + ggplot2::coord_flip() +
          ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = 10)) # 横向き用にy軸テキスト調整
      }

      # 色の一貫性を保つためのスケール適用
      if (input$plot_type %in% c("stacked", "dodged")) {
        # 元データの全グループ取得 (NA/空除外後)
        all_groups <- tryCatch(sort(unique(reactive_df_raw()[[grouping_var]])), error = function(e) NULL)
        if (!is.null(all_groups) && length(all_groups) > 0) {
          color_palette <- scales::hue_pal()(length(all_groups))
          names(color_palette) <- all_groups
          # プロットデータに含まれる有効なレベルの色のみを適用
          valid_colors <- color_palette[levels(plot_data$group_col)]
          # 有効な色が存在する場合のみスケールを適用
          if (length(valid_colors) > 0 && !all(is.na(valid_colors))) {
            p <- p + ggplot2::scale_fill_manual(values = valid_colors, name = fill_label, drop = FALSE)
          }
        }
      }
      return(p)
    })

    # プロットのレンダリング
    # 7. プロット出力
    output$plot <- renderPlot(
      {
           plot_obj()
      },
      width = reactive(input$plot_width),
      height = reactive(input$plot_height),
    )

      # plot_obj()がNULLやエラーを返す可能性を考慮
      # plot_result <- tryCatch(
      #   {
#          plot_obj()
        # },
      # width = reactive(input$plot_width),
      # height = reactive(input$plot_height),
      #   error = function(e) {
      #     showNotification(paste("Error generating plot:", e$message), type = "error", duration = 10)
      #     # エラー時に表示する空のプロット
      #     ggplot() +
      #       theme_void() +
      #       geom_text(aes(0, 0, label = "Error generating plot. Check console/logs."), size = 5)
      #   }
#      )

    #   # validateの結果でNULLが返る場合も考慮
    #   if (is.null(plot_result)) {
    #     # 何も表示しないか、メッセージを表示
    #     plot.new()
    #     text(0.5, 0.5, "No data to display for current selections.", cex = 1.2)
    #   } else {
    #     print(plot_result)
    #   }
    # })

        # PDFダウンロードハンドラー
    output$download_plot <- downloadHandler(
      filename = function() {
        "plot.pdf"
      }, # ファイル名を修正
      content = function(file) {
        ggsave(file, plot = plot_obj(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300) # renderCustomPlot() を plot() に修正
      }
    )


    # 8. テーブル出力
    output$table <- renderDT({
      dt_data_full <- tryCatch(
        {
          table_for_plot()
        },
        error = function(e) {
          showNotification(paste("Error generating table data:", e$message), type = "error", duration = 10)
          return(NULL) # エラー時はNULLを返す
        }
      )

      # table_for_plot()がNULLを返す場合やエラーの場合
      validate(need(!is.null(dt_data_full), "No data available for the table with current selections."))

      # テーブル表示用にデータを整形
      dt_data <- dt_data_full %>%
        # total_n 列は通常不要なので削除 (もし存在すれば)
        dplyr::select(-any_of("total_n"))

      # グループ化列名を元に戻す
      if ("group_col" %in% names(dt_data)) {
        dt_data <- dt_data %>% dplyr::rename(!!input$group_by := group_col)
      }

      # 表示値に応じて列を選択・フォーマット
      if (input$display_value == "count") {
        dt_data <- dt_data %>% dplyr::select(-percentage)
      } else { # percentage
        dt_data <- dt_data %>%
          dplyr::mutate(percentage = scales::percent(percentage, accuracy = 0.1)) %>%
          dplyr::select(-n)
      }

      datatable(dt_data,
        options = list(scrollX = TRUE, pageLength = 5, searching = FALSE, lengthChange = FALSE),
        rownames = FALSE
      )
    })

  output$download_table <- downloadHandler(
      filename = function() {
        "table.csv"
      },
      content = function(file) {
        write.csv(table_for_plot(), file, row.names = FALSE)
      }
    )

  }) # moduleServer終了
} # summaryPlotServer終了
