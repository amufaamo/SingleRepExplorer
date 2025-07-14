geneUsageUI <- function(id) { # 関数名変更
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      # データタイプ選択
      vdjType(ns),
      # ★ 集計対象の遺伝子列を選択 ★
      selectInput(ns("target_gene_column"), "Target Gene Column",
        choices = NULL, # Server側で設定するため NULL に
        selected = NULL
      ),
      # グループ化変数選択
      groupByInput(ns),
      # プロットタイプ選択
      selectInput(ns("plot_type"), "Plot Type",
        choices = c(
          "Barplot (Dodged)" = "dodged", # グループ比較、横並び
          "Barplot (Stacked)" = "stacked", # グループ比較、積み上げ
          "Barplot (Overall)" = "normal", # 全体での頻度
          "Heatmap" = "heatmap" # グループ比較、ヒートマップ
        ),
        selected = "dodged"
      ),
      selectInput(ns("display_value"), "Display Value",
        choices = c("Count" = "count", "Percentage" = "percentage"),
        selected = "count"
      ),
      uiOutput(ns("filter_groups_ui")),
      # ★ Horizontal オプションは Barplot系のみ表示 ★
      conditionalPanel(
        condition = "input.plot_type == 'dodged' || input.plot_type == 'stacked' || input.plot_type == 'normal'",
        ns = ns, # conditionalPanel 内で ns を指定
        checkboxInput(ns("horizontal"), "Horizontal Barplot", value = FALSE) # ラベル変更
      ),
      # (Heatmap用オプションは任意で追加)
      commonPlotOptions(ns),
    ),
    mainPanel(
      # プロット出力
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      h3("Table"),
      downloadButton(ns("download_table"), "Download table (.csv)"),
      DTOutput(ns("table"))
    )
  )
}
geneUsageServer <- function(id, myReactives) { # 関数名変更
  moduleServer(id, function(input, output, session) {
    # Seuratオブジェクト更新時の処理 (変更なし)
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
      # TODO: target_gene_column の選択肢更新 observeEvent
    })



    # vdjType 変更時に Target Gene Column の選択肢を更新
    observeEvent(input$vdj_type,
      {
        req(input$vdj_type, input$vdj_type != "")

        target_choices <- list()
        selected_choice <- NULL
        # ★ 新しいV-Jペアの列名を定義 ★
        vj_pair_names <- list(
          tcr_tra = "TCR_TRA_v_j_pair",
          tcr_trb = "TCR_TRB_v_j_pair",
          bcr_igh = "BCR_IGH_v_j_pair",
          bcr_igl = "BCR_IGL_v_j_pair"
        )

        if (input$vdj_type == "tcr") {
          target_choices <- c(
            # 単一遺伝子
            "TRA V Gene" = "TCR_TRA_v_gene",
            "TRA J Gene" = "TCR_TRA_j_gene",
            "TRB V Gene" = "TCR_TRB_v_gene",
            "TRB D Gene" = "TCR_TRB_d_gene",
            "TRB J Gene" = "TCR_TRB_j_gene",
            # V-J ペア
            "TRA V-J Pair" = vj_pair_names$tcr_tra, # 定義した列名を使用
            "TRB V-J Pair" = vj_pair_names$tcr_trb # 定義した列名を使用
          )
          selected_choice <- "TCR_TRB_v_gene"
        } else if (input$vdj_type == "bcr") {
          target_choices <- c(
            # 単一遺伝子
            "IGH V Gene" = "BCR_IGH_v_gene",
            "IGH D Gene" = "BCR_IGH_d_gene",
            "IGH J Gene" = "BCR_IGH_j_gene",
            "IGH C Gene" = "BCR_IGH_c_gene",
            "IGK V Gene" = "BCR_IGK_v_gene",
            "IGK J Gene" = "BCR_IGK_j_gene",
            "IGL V Gene" = "BCR_IGL_v_gene",
            "IGL J Gene" = "BCR_IGL_j_gene",
            # V-J ペア
            "IGH V-J Pair" = vj_pair_names$bcr_igh,
            "IGK V-J Pair" = vj_pair_names$bcr_igk, # IGK追加
            "IGL V-J Pair" = vj_pair_names$bcr_igl
          )
          selected_choice <- "BCR_IGH_v_gene"
        } else {
          target_choices <- list("Select VDJ Type First" = "")
          selected_choice <- ""
        }

        # updateSelectInput で選択肢を更新
        updateSelectInput(session, "target_gene_column",
          choices = target_choices,
          selected = selected_choice
        )
      },
      ignoreNULL = FALSE,
      ignoreInit = FALSE
    )

    # 1. 元データのリアクティブ (★ V-J ペア列生成を追加 ★)
    reactive_df_raw <- reactive({
      req(input$vdj_type, input$group_by, input$target_gene_column)
      validate(need(input$target_gene_column != "", "Please select a Target Gene Column."))

      df_orig <- NULL # 元のDFを保持
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df_orig <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df_orig <- myReactives$bcr_df
      } else {
        return(NULL)
      }
      validate(need(nrow(df_orig) > 0, "Input data empty."))

      # --- V-Jペア列の生成 ---
      df_processed <- df_orig # コピーして処理

      # V-Jペア生成関数 (NA処理込み、列存在チェック強化)
      create_vj_pair <- function(df, v_col, j_col, new_col_name) {
        # 両方の列が存在するかチェック
        if (all(c(v_col, j_col) %in% names(df))) {
          df %>% dplyr::mutate(
            !!new_col_name := dplyr::if_else(
              # NAだけでなく空文字もチェック
              !is.na(.data[[v_col]]) & .data[[v_col]] != "" &
                !is.na(.data[[j_col]]) & .data[[j_col]] != "",
              paste(.data[[v_col]], .data[[j_col]], sep = "-"),
              NA_character_ # VかJがNA/空ならペアもNA
            )
          )
        } else {
          # 必要な列がない場合、警告を出してペア列は追加しない（または全てNAの列を追加）
          warning(paste(
            "Cannot create V-J pair '", new_col_name, "'. Missing required columns:",
            paste(setdiff(c(v_col, j_col), names(df)), collapse = ", ")
          ))
          df # 元のdfを返す
          # もし列がない場合でも列自体は追加したいなら:
          # df %>% dplyr::mutate(!!new_col_name := NA_character_)
        }
      }

      # observeEvent で定義した V-J ペア列名を再利用
      vj_pair_names <- list(
        tcr_tra = "TCR_TRA_v_j_pair", tcr_trb = "TCR_TRB_v_j_pair",
        bcr_igh = "BCR_IGH_v_j_pair", bcr_igk = "BCR_IGK_v_j_pair", bcr_igl = "BCR_IGL_v_j_pair"
      )

      # vdjType に応じてペア列を生成
      if (input$vdj_type == "tcr") {
        df_processed <- df_processed %>%
          create_vj_pair("TCR_TRA_v_gene", "TCR_TRA_j_gene", vj_pair_names$tcr_tra) %>%
          create_vj_pair("TCR_TRB_v_gene", "TCR_TRB_j_gene", vj_pair_names$tcr_trb)
      } else if (input$vdj_type == "bcr") {
        df_processed <- df_processed %>%
          create_vj_pair("BCR_IGH_v_gene", "BCR_IGH_j_gene", vj_pair_names$bcr_igh) %>%
          create_vj_pair("BCR_IGK_v_gene", "BCR_IGK_j_gene", vj_pair_names$bcr_igk) %>% # IGK追加
          create_vj_pair("BCR_IGL_v_gene", "BCR_IGL_j_gene", vj_pair_names$bcr_igl)
      }
      # --- V-Jペア列の生成 終了 ---

      # 選択された列 (group_by, target_gene_column) が最終的なDFに存在するか確認
      required_cols <- c(input$group_by, input$target_gene_column)
      actual_cols <- names(df_processed)
      missing_cols <- setdiff(required_cols, actual_cols)
      validate(need(
        length(missing_cols) == 0,
        paste(
          "Selected column(s) not found after processing:",
          paste(shQuote(missing_cols), collapse = ", "),
          "\nAvailable columns:", paste(actual_cols, collapse = ", ")
        )
      )) # 利用可能な列も表示

      # NA/空文字を除外 (group_by と target_gene_column)
      df_filtered_na_empty <- df_processed %>%
        dplyr::filter(
          !is.na(.data[[input$group_by]]) & .data[[input$group_by]] != "",
          !is.na(.data[[input$target_gene_column]]) & .data[[input$target_gene_column]] != ""
        )
      validate(need(
        nrow(df_filtered_na_empty) > 0,
        paste(
          "No data remaining after removing NA/empty values from columns:",
          shQuote(input$group_by), "and", shQuote(input$target_gene_column)
        )
      ))
      return(df_filtered_na_empty)
    })

    # 2. フィルタリング前の合計数 (変更なし)
    original_totals <- reactive({
      # ... (変更なし) ...
      df <- reactive_df_raw()
      req(df, input$group_by)
      grouping_var <- input$group_by
      group_tots <- df %>%
        dplyr::count(.data[[grouping_var]], name = "original_group_total") %>%
        dplyr::rename(group_col = .data[[grouping_var]]) %>%
        dplyr::mutate(group_col = as.character(group_col))
      overall_tot <- nrow(df)
      return(list(group = group_tots, overall = overall_tot))
    })

    # 3. グループフィルタリングUI (変更なし)
    output$filter_groups_ui <- renderUI({
      # ... (変更なし) ...
      df <- reactive_df_raw()
      req(df, input$group_by)
      grouping_var <- input$group_by
      available_groups <- sort(unique(df[[grouping_var]]))
      validate(need(length(available_groups) > 0, "No available groups."))
      checkboxGroupInput(session$ns("filter_groups"), paste("Filter", tools::toTitleCase(gsub("_", " ", grouping_var)), ":"),
        choices = available_groups, selected = available_groups, inline = TRUE
      )
    })

    # 4. フィルタリングされたデータフレーム (変更なし)
    filtered_df <- reactive({
      # ... (変更なし) ...
      df <- reactive_df_raw()
      req(df, input$group_by, input$filter_groups)
      grouping_var <- input$group_by
      selected_groups <- input$filter_groups
      df_filtered <- df %>% dplyr::filter(.data[[grouping_var]] %in% selected_groups)
      if (nrow(df_filtered) == 0) {
        return(NULL)
      } # エラー処理省略
      return(df_filtered)
    })

    # ★ 並び替え基準グループの observe ブロックを削除 ★
    # observe({ ... })
    # --- 5. プロット/テーブル用データの計算 (Percentage計算部分を修正) ---
    table_for_plot <- reactive({
        df_filtered <- filtered_df()
        req(df_filtered, input$plot_type %in% c("normal", "stacked", "dodged", "heatmap"),
            input$group_by, input$target_gene_column, input$target_gene_column != "",
            input$filter_groups)

        og_totals <- original_totals(); req(og_totals)
        plot_data <- NULL
        grouping_var <- input$group_by; target_col <- input$target_gene_column

        # 全てのターゲット遺伝子を自然順でソート
        all_target_genes_sorted <- df_filtered %>%
            dplyr::distinct(.data[[target_col]]) %>%
            dplyr::pull(.data[[target_col]]) %>%
            stringr::str_sort(numeric = TRUE)
        validate(need(length(all_target_genes_sorted) > 0, "No target genes found after filtering."))
        gene_levels_natural_order <- factor(all_target_genes_sorted, levels = all_target_genes_sorted)

        # 全データでカウント (group_by も含める)
        counts_by_group_raw <- df_filtered %>%
            dplyr::count(.data[[target_col]], .data[[grouping_var]], name = "n") %>%
            dplyr::rename(group_col = .data[[grouping_var]])

        # 有効なグループ (フィルター後)
        # input$group_by が選択されているか確認 (Normal Plot以外で必要)
        # Normal Plotの場合、group_col には単一の値 (e.g., "All") が入るように工夫が必要かも
        # → ここではinput$filter_groupsを使うことで、Normalでもエラーなく進む想定
        valid_group_levels <- intersect(input$filter_groups, unique(counts_by_group_raw$group_col))
        validate(need(length(valid_group_levels) > 0, "No valid groups found after filtering."))
        group_levels_to_complete <- factor(valid_group_levels, levels = valid_group_levels) # Use sorted unique groups

        # 全組み合わせを生成 (全遺伝子 x 全有効グループ)
        plot_data_completed <- counts_by_group_raw %>%
            dplyr::mutate(group_col = as.character(group_col)) %>% # Ensure character for complete
            tidyr::complete(
                !!target_col := gene_levels_natural_order,
                group_col = levels(group_levels_to_complete),
                fill = list(n = 0)
            ) %>%
            tidyr::drop_na(!!target_col)

        # 元の合計数を取得
        original_group_totals_df <- og_totals$group; req(original_group_totals_df)
        overall_total_original <- og_totals$overall; req(overall_total_original)

        # ★ Percentage計算: グループ内比率と全体比率の両方を計算 ★
        plot_data_with_perc <- plot_data_completed %>%
            dplyr::left_join(original_group_totals_df, by = "group_col") %>%
            dplyr::mutate(
                percentage_within_group = ifelse(is.na(original_group_total) | original_group_total == 0, 0, n / original_group_total),
                percentage_of_total = ifelse(overall_total_original == 0, 0, n / overall_total_original)
            ) %>%
            dplyr::select(-original_group_total)

        validate(need(nrow(plot_data_with_perc) > 0, "Data processing empty."))

        # Factor化
        plot_data <- plot_data_with_perc %>%
            dplyr::mutate(
                !!target_col := factor(.data[[target_col]], levels = levels(gene_levels_natural_order)),
                group_col = factor(group_col, levels = levels(group_levels_to_complete))
            ) %>%
            dplyr::arrange(.data[[target_col]], group_col)

        # display_value に応じて使用する percentage 列を決定 (plot_obj側で処理するため、ここでは両方保持)
        # ただし、後のテーブル表示用に 'percentage' 列はグループ内比率にしておく
        if("percentage_within_group" %in% names(plot_data)) {
             plot_data <- plot_data %>% dplyr::rename(percentage = percentage_within_group)
        }


        validate(need(!is.null(plot_data) && nrow(plot_data) > 0, "No data for plot/table."))
        return(plot_data)
    })

    # --- 6. プロットオブジェクトの生成 (★ Heatmap 対応追加 ★) ---
    plot_obj <- reactive({
        plot_data <- table_for_plot()
        req(plot_data)
        grouping_var <- input$group_by
        target_col <- input$target_gene_column

        # --- Heatmap Plot ---
        if (input$plot_type == "heatmap") {
            # Heatmap はグループ化データが前提
            validate(need("group_col" %in% names(plot_data) && length(levels(plot_data$group_col)) > 0,
                          "Heatmap requires grouped data (select a grouping variable)."))
             # Normal plot が選択された場合は Heatmap を表示しない方が親切かもしれない
             # validate(need(input$group_by != "none" || length(levels(plot_data$group_col)) > 1), "Select a grouping variable for Heatmap")


            # 表示値 (fill aesthetic)
            fill_var <- if(input$display_value == "count") "n" else "percentage" # 'percentage' はグループ内比率
            fill_label <- if(input$display_value == "count") "Count" else "Percentage (within Group)"

            # Y軸のレベル (遺伝子/ペア、自然順)
            y_levels <- levels(plot_data[[target_col]])
            # X軸のレベル (グループ)
            x_levels <- levels(plot_data$group_col)
            validate(need(length(y_levels)>0 && length(x_levels)>0, "Cannot determine axes for heatmap."))

            p <- ggplot2::ggplot(plot_data, aes(x = group_col, y = .data[[target_col]], fill = .data[[fill_var]])) +
                ggplot2::geom_tile(color = "grey85", linewidth = 0.2) + # タイル境界線調整
                ggplot2::scale_y_discrete(limits = rev(y_levels)) + # ヒートマップは通常Y軸逆順
                ggplot2::scale_x_discrete(limits = x_levels, position = "bottom") +
                # ★ 色スケール: 0が白、値が大きいほど濃い青になるように ★
                ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = fill_label,
                                             labels = if(fill_var == "percentage") scales::percent else waiver(),
                                             na.value = "grey80",
                                             guide = guide_colorbar(barwidth = 0.8, barheight = 10)) + # カラーバー調整
                ggplot2::labs(
                    title = paste("Gene Usage Frequencies Heatmap"),
                    x = tools::toTitleCase(gsub("_", " ", grouping_var)),
                    y = gsub("_", " ", tools::toTitleCase(target_col)),
                    fill = fill_label
                ) +
                ggplot2::theme_minimal(base_size = 11) + # ベースフォントサイズ調整
                ggplot2::theme(
                    axis.text.x = element_text(angle = 45, hjust = 1, size=9), # X軸ラベル調整
                    axis.text.y = element_text(size = 8),                    # Y軸ラベル調整
                    axis.ticks = element_blank(),
                    panel.grid = element_blank(),
                    legend.position = if (isTRUE(input$legend %||% TRUE)) "right" else "none",
                    plot.title = element_text(hjust = 0.5) # タイトル中央揃え
                )

             # 値表示オプション（もし追加する場合）
             # if (isTRUE(input$heatmap_show_values)) { ... }

            return(p)

        # --- Barplot (Dodged, Stacked, Normal) ---
        } else if (input$plot_type %in% c("dodged", "stacked", "normal")) {

            # y軸関連
            if (input$display_value == "count") {
                y_var <- "n"; y_label <- "Frequency (Number of Cells)"
                y_scale <- ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
            } else { # percentage
                # Normal Plot の場合は全体比率を使用、それ以外はグループ内比率
                if(input$plot_type == "normal" && "percentage_of_total" %in% names(plot_data)) {
                    y_var <- "percentage_of_total"
                    y_label <- "Frequency (% of Total Cells)"
                } else {
                    y_var <- "percentage" # グループ内比率 (rename済み)
                    y_label <- paste("Frequency (% within Original", tools::toTitleCase(gsub("_", " ", grouping_var)), ")")
                }
                # percentage 列が存在するか確認
                validate(need(y_var %in% names(plot_data), paste("Required percentage column", shQuote(y_var), "not found.")))
                y_scale <- ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.05)))
            }

            # 基本プロット
            p <- ggplot2::ggplot(plot_data)
            fill_label <- if (input$plot_type %in% c("stacked", "dodged")) { tools::toTitleCase(gsub("_", " ", grouping_var)) } else { NULL }

            # geom と aes
            if (input$plot_type == "normal") {
                 # y_var が percentage_of_total か n になる
                p <- p + ggplot2::aes(x = .data[[target_col]], y = .data[[y_var]]) +
                    ggplot2::geom_bar(stat = "identity", fill = "steelblue")
            } else { # stacked or dodged
                 # y_var が percentage (within group) か n になる
                validate(need("group_col" %in% names(plot_data), "Stacked/Dodged plot requires a grouping variable."))
                req(levels(plot_data$group_col))
                p <- p + ggplot2::aes(x = .data[[target_col]], y = .data[[y_var]], fill = group_col)
                if (input$plot_type == "stacked") {
                    p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")
                } else { # dodged
                    p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(preserve = "single"))
                }
            }

            # タイトル
            plot_title <- "Gene Usage Frequencies"
            # Check if grouping variable exists and is not 'none' or similar placeholder if you add one
            is_grouped <- "group_col" %in% names(plot_data) && length(levels(plot_data$group_col)) > 1
            all_groups_original <- tryCatch(unique(reactive_df_raw()[[input$group_by]]), error = function(e) NULL)
            # Only add "(Groups Filtered)" if it's actually grouped and filtered
            if (is_grouped && !is.null(all_groups_original) && length(input$filter_groups) < length(all_groups_original)) {
                plot_title <- paste(plot_title, "(Groups Filtered)")
            }

            # テーマ、ラベル、スケール
            p <- p +
                ggplot2::theme_classic(base_size = 14) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                    axis.title = ggplot2::element_text(size = 12),
                    plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
                    legend.position = if (isTRUE(input$legend %||% TRUE)) "right" else "none"
                ) +
                y_scale +
                ggplot2::labs( title = plot_title, x = gsub("_", " ", tools::toTitleCase(target_col)),
                               y = y_label, fill = fill_label)

            # 横向き (Barplotのみ)
            if (isTRUE(input$horizontal)) {
                 p <- p + ggplot2::coord_flip() + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = 10))
            }

            # 色 (Stacked/Dodgedのみ)
            if (input$plot_type %in% c("stacked", "dodged")) {
                # ...(色設定コードは変更なし)...
                 all_groups <- tryCatch(sort(unique(reactive_df_raw()[[grouping_var]])), error = function(e) NULL)
                 if (!is.null(all_groups) && length(all_groups) > 0) {
                     color_palette <- scales::hue_pal()(length(all_groups))
                     names(color_palette) <- all_groups
                     req(plot_data$group_col)
                     valid_colors <- color_palette[levels(plot_data$group_col)]
                     if (length(valid_colors) > 0 && !all(is.na(valid_colors))) {
                         p <- p + ggplot2::scale_fill_manual(values = valid_colors, name = fill_label, drop = FALSE)
                     }
                 }
            }
            return(p)
        } else {
            # input$plot_type が予期せぬ値の場合
            warning("Invalid plot_type selected: ", input$plot_type)
            return(ggplot() + theme_void() + labs(title = "Invalid Plot Type Selected"))
        }
    })

    # 7. プロット出力 (renderPlot, downloadHandler) (変更なし、エラーハンドリング強化)
    output$plot <- renderPlot(
      {
        # table_for_plot() が NULL を返す可能性を考慮
        plot_data_check <- tryCatch(table_for_plot(), error = function(e) NULL)
        if (is.null(plot_data_check) || nrow(plot_data_check) == 0) {
          plot.new()
          text(0.5, 0.5, "No data to display for current selections.", cex = 1.2)
          return() # プロット生成をスキップ
        }
        # プロット生成
        plot_result <- tryCatch(plot_obj(), error = function(e) {
          showNotification(paste("Error generating plot:", e$message), type = "error", duration = 10)
          ggplot() +
            theme_void() +
            geom_text(aes(0, 0, label = "Error generating plot."), size = 5)
        })
        if (!is.null(plot_result)) print(plot_result)
      },
      # commonPlotOptions からの入力を取得 (%||% でデフォルト値設定)
      width = reactive(input$plot_width %||% 600),
      height = reactive(input$plot_height %||% 500)
    )

    output$download_plot <- downloadHandler(
      filename = function() {
        # ファイル名 (変更なし)
        paste0("gene_usage_", input$target_gene_column, "_", input$plot_type, ".pdf")
      },
      content = function(file) {
        p <- plot_obj()
        req(p)
        # プロットサイズ (変更なし)
        plot_width_inch <- (input$plot_width %||% 600) / 72
        plot_height_inch <- (input$plot_height %||% 500) / 72
        ggsave(file, plot = p, width = plot_width_inch, height = plot_height_inch, device = "pdf", dpi = 300)
      }
    )

    # 8. テーブル出力 (★ DT オプション変更の可能性 ★)
    output$table <- renderDT({
      dt_data_full <- tryCatch(table_for_plot(), error = function(e) {
        showNotification(paste("Error generating table data:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      validate(need(!is.null(dt_data_full) && nrow(dt_data_full) > 0, "No data available for the table.")) # nrow > 0 を追加

      # テーブル表示用に整形 (変更なし)
      dt_data <- dt_data_full %>%
        dplyr::select(-any_of(c("total_n", "total_n_filtered", "sort_value")))
      target_col <- input$target_gene_column
      grouping_var <- input$group_by
      if ("group_col" %in% names(dt_data)) {
        dt_data <- dt_data %>% dplyr::rename(!!grouping_var := group_col)
      }
      if (input$display_value == "count") {
        dt_data <- dt_data %>% dplyr::select(-any_of("percentage"))
      } else { # percentage
        dt_data <- dt_data %>%
          dplyr::mutate(percentage = scales::percent(percentage, accuracy = 0.1)) %>%
          dplyr::select(-any_of("n"))
      }
      ordered_cols <- c(
        target_col, if (input$plot_type != "normal") grouping_var else NULL,
        if (input$display_value == "count") "n" else "percentage"
      )
      ordered_cols_exist <- intersect(ordered_cols, names(dt_data))
      dt_data <- dt_data %>% dplyr::select(all_of(ordered_cols_exist))

      # ★ DT オプション: 全要素表示のため pageLength を調整、検索有効化 ★
      datatable(dt_data,
        options = list(
          scrollX = TRUE,
          pageLength = 25, # または -1 で全件表示 (多いと重い)
          lengthMenu = c(10, 25, 50, 100, -1), # -1 = All
          searching = TRUE, # 検索を有効に
          lengthChange = TRUE # 表示件数変更を有効に
        ),
        rownames = FALSE
      )
    })

    # テーブルダウンロード (変更なし)
    output$download_table <- downloadHandler(
      filename = function() {
        # ファイル名 (変更なし)
        paste0("gene_usage_summary_", input$target_gene_column, "_", input$group_by, ".csv")
      },
      content = function(file) {
        data_to_write <- table_for_plot()
        # 列名整形 (変更なし)
        if ("group_col" %in% names(data_to_write) && input$plot_type != "normal") {
          data_to_write <- data_to_write %>% dplyr::rename(!!input$group_by := group_col)
        }
        data_to_write <- data_to_write %>% dplyr::select(-any_of(c("total_n_filtered", "sort_value")))
        write.csv(data_to_write, file, row.names = FALSE)
      }
    )
  }) # moduleServer終了
} # geneUsageServer終了


# --- UI ---

analysisPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      # --- 基本設定 ---
      # 1. 解析ターゲットの選択（統合の核となる部分）
      selectInput(ns("analysis_target"), "Analysis Target",
                  choices = c("Clonotype Frequency" = "clonotype", "Gene Usage" = "gene_usage"),
                  selected = "clonotype"),

      # 2. 共通の入力項目
      vdjType(ns),
      groupByInput(ns),

      # --- 解析ターゲットに応じたオプション ---
      # 3a. Gene Usage選択時にのみ表示
      conditionalPanel(
        condition = "input.analysis_target == 'gene_usage'",
        ns = ns,
        selectInput(ns("target_gene_column"), "Target Gene Column", choices = NULL) # サーバー側で選択肢を更新
      ),
      
      # 3b. Clonotype Frequency選択時にのみ表示
      conditionalPanel(
        condition = "input.analysis_target == 'clonotype'",
        ns = ns,
        numericInput(ns("top_n"), "Top N Clonotypes:", value = 10, min = 1, step = 1)
      ),

      # --- プロット関連の設定 ---
      # 4. プロットタイプの選択
      # Note: サーバー側で`analysis_target`に応じて選択肢を更新することを推奨
      selectInput(ns("plot_type"), "Plot Type",
                  choices = c(
                    "Barplot (Dodged)" = "dodged",
                    "Barplot (Stacked)" = "stacked",
                    "Barplot (Normal)" = "normal",
                    "Heatmap" = "heatmap" # Gene Usageでのみ利用
                  ),
                  selected = "dodged"),
                  
      # 5. 表示形式の選択 (共通)
      selectInput(ns("display_value"), "Display Value",
                  choices = c("Count" = "count", "Percentage" = "percentage"),
                  selected = "count"),
      
      # 6. グループフィルタリング (共通)
      uiOutput(ns("filter_groups_ui")),

      # --- プロット表示オプション ---
      # 7a. Dodged Barplot選択時にX軸のソート順を指定
      conditionalPanel(
        condition = "input.plot_type == 'dodged'",
        ns = ns,
        selectInput(ns("sort_by_group"), "Sort X-axis by group:", choices = NULL) # サーバー側で選択肢を更新
      ),
      
      # 7b. Barplot系選択時にグラフの向きを指定
      conditionalPanel(
        condition = "input.plot_type == 'dodged' || input.plot_type == 'stacked' || input.plot_type == 'normal'",
        ns = ns,
        checkboxInput(ns("horizontal"), "Horizontal Barplot", value = FALSE)
      ),

      # 8. 共通のプロットオプション (共通)
      commonPlotOptions(ns),
    ),
    mainPanel(
      # --- 出力部分 (共通) ---
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      h3("Table"),
      downloadButton(ns("download_table"), "Download table (.csv)"),
      DTOutput(ns("table"))
    )
  )
}

# --- サーバー ---
analysisPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {

    # 'analysis_target' の選択に応じて処理を分岐
    observeEvent(input$analysis_target, {
      # "Clonotype Frequency" が選択された場合のみ、以下の処理を実行
      if (input$analysis_target == "clonotype") {

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
            dplyr::mutate(group_col = as.character(group_col))

          # 全体の合計数
          overall_tot <- nrow(df)
          
          return(list(group = group_tots, overall = overall_tot))
        })

        # 3. グループフィルタリングUIの生成
        output$filter_groups_ui <- renderUI({
          df <- reactive_df_raw()
          req(df, input$group_by) # dfがNULLでないこと、group_byが指定されていること
          grouping_var <- input$group_by
          
          available_groups <- sort(unique(df[[grouping_var]]))

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
          req(df, input$group_by, input$filter_groups)

          grouping_var <- input$group_by
          selected_groups <- input$filter_groups
          
          df_filtered <- df %>%
            dplyr::filter(.data[[grouping_var]] %in% selected_groups)
          
          if (nrow(df_filtered) == 0) {
            showNotification("No data remaining after filtering by selected groups.", type = "warning", duration = 5)
            return(NULL)
          }
          return(df_filtered)
        })

        # ★ 追加: 並び替え基準グループの選択肢を更新
        observe({
          req(input$plot_type == "dodged", input$filter_groups)
          choices_with_overall <- c("Overall (All selected groups)" = "overall", input$filter_groups)
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

        # 5. プロット/テーブル用データの計算
        table_for_plot <- reactive({
          df_filtered <- filtered_df()
          req(df_filtered, input$plot_type, input$top_n, input$top_n >= 1, input$group_by, input$filter_groups)
          og_totals <- original_totals()
          req(og_totals)

          plot_data <- NULL
          n_top <- as.integer(input$top_n)
          grouping_var <- input$group_by
          
          # (元のコードをそのまま挿入...)
          # --- Normal Plot ---
          if (input$plot_type == "normal") {
              plot_data_raw <- df_filtered %>%
                dplyr::count(raw_clonotype_id, sort = TRUE, name = "n")

              total_cells_original <- og_totals$overall
              req(total_cells_original) # 元の合計数が存在すること

              plot_data <- plot_data_raw %>%
                dplyr::slice_head(n = n_top) %>%
                dplyr::mutate(percentage = ifelse(total_cells_original == 0, 0, n / total_cells_original)) %>%
                dplyr::mutate(raw_clonotype_id = factor(raw_clonotype_id, levels = .$raw_clonotype_id))
          # --- Stacked/Dodged Plot ---
          } else if (input$plot_type %in% c("stacked", "dodged")) {
              counts_by_group_raw <- df_filtered %>%
                dplyr::count(raw_clonotype_id, .data[[grouping_var]], name = "n") %>%
                dplyr::rename(group_col = .data[[grouping_var]])
              
              clonotype_order_info <- NULL
              if (input$plot_type == "dodged") {
                req(input$sort_by_group)
                sort_criterion_group <- input$sort_by_group
              } else {
                sort_criterion_group <- "overall"
              }
              
              all_clonotypes <- distinct(counts_by_group_raw, raw_clonotype_id)

              if (sort_criterion_group == "overall") {
                clonotype_order_info <- counts_by_group_raw %>%
                  dplyr::group_by(raw_clonotype_id) %>%
                  dplyr::summarise(sort_value = sum(n), .groups = "drop") %>%
                  dplyr::right_join(all_clonotypes, by = "raw_clonotype_id") %>% 
                  tidyr::replace_na(list(sort_value = 0)) %>% 
                  dplyr::arrange(dplyr::desc(sort_value))
              } else {
                clonotype_order_info <- counts_by_group_raw %>%
                  dplyr::filter(group_col == sort_criterion_group) %>%
                  dplyr::select(raw_clonotype_id, sort_value = n) %>%
                  dplyr::right_join(all_clonotypes, by = "raw_clonotype_id") %>%
                  tidyr::replace_na(list(sort_value = 0)) %>%
                  dplyr::arrange(dplyr::desc(sort_value))
              }
              
              validate(need(!is.null(clonotype_order_info) && nrow(clonotype_order_info) > 0, "Could not determine clonotype order."))
              
              top_clonotypes_sorted <- clonotype_order_info %>%
                dplyr::slice_head(n = n_top)
              
              validate(need(nrow(top_clonotypes_sorted) > 0, "No top clonotypes found based on the sorting criterion."))
              
              top_clonotype_ids_new_order <- top_clonotypes_sorted$raw_clonotype_id
              
              valid_group_levels <- intersect(input$filter_groups, unique(counts_by_group_raw$group_col))
              validate(need(length(valid_group_levels) > 0, "No valid groups found after filtering."))

              counts_by_group_filtered <- counts_by_group_raw %>%
                dplyr::filter(raw_clonotype_id %in% top_clonotype_ids_new_order) %>%
                dplyr::mutate(group_col = factor(as.character(group_col), levels = valid_group_levels))

              group_levels_to_complete <- factor(valid_group_levels, levels = valid_group_levels)
              clonotype_levels_to_complete <- factor(top_clonotype_ids_new_order, levels = top_clonotype_ids_new_order)

              plot_data_completed <- counts_by_group_filtered %>%
                tidyr::complete(
                  raw_clonotype_id = clonotype_levels_to_complete,
                  group_col = group_levels_to_complete, fill = list(n = 0)
                ) %>%
                tidyr::drop_na(raw_clonotype_id) %>%
                dplyr::mutate(group_col = as.character(group_col))
              
              original_group_totals_df <- og_totals$group
              req(original_group_totals_df)

              plot_data_with_perc <- plot_data_completed %>%
                dplyr::left_join(original_group_totals_df, by = "group_col") %>%
                dplyr::mutate(percentage = ifelse(is.na(original_group_total) | original_group_total == 0, 0, n / original_group_total)) %>%
                dplyr::select(-original_group_total)
              
              validate(need(nrow(plot_data_with_perc) > 0, "Data processing empty."))
              
              plot_data <- plot_data_with_perc %>%
                dplyr::left_join(
                  counts_by_group_raw %>%
                    dplyr::group_by(raw_clonotype_id) %>%
                    dplyr::summarise(total_n_filtered = sum(n), .groups = "drop"),
                  by = "raw_clonotype_id"
                ) %>%
                dplyr::mutate(
                  raw_clonotype_id = factor(raw_clonotype_id, levels = top_clonotype_ids_new_order),
                  group_col = factor(group_col, levels = valid_group_levels)
                ) %>%
                dplyr::mutate(total_n_filtered = coalesce(total_n_filtered, 0))
          } else { 
              warning(paste("Unexpected plot_type:", input$plot_type))
              return(NULL)
          }

          validate(need(!is.null(plot_data) && nrow(plot_data) > 0, "No data for plot/table."))
          return(plot_data)
        })

        # 6. プロットオブジェクトの生成
        plot_obj <- reactive({
          plot_data <- table_for_plot()
          req(plot_data)
          grouping_var <- input$group_by
          # (元のコードをそのまま挿入...)
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
          
          p <- ggplot2::ggplot(plot_data)
          fill_label <- if (input$plot_type %in% c("stacked", "dodged")) {
              tools::toTitleCase(gsub("_", " ", grouping_var))
          } else {
              NULL
          }
          
          if (input$plot_type == "normal") {
              p <- p + ggplot2::aes(x = raw_clonotype_id, y = .data[[y_var]]) +
              ggplot2::geom_bar(stat = "identity", fill = "steelblue")
          } else if (input$plot_type %in% c("stacked", "dodged")) {
              valid_levels <- levels(plot_data$group_col)
              validate(need(length(valid_levels) > 0, "No valid group levels for plotting."))
              
              p <- p + ggplot2::aes(x = raw_clonotype_id, y = .data[[y_var]], fill = group_col)
              
              if (input$plot_type == "stacked") {
                  p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")
              } else { # dodged
                  p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(preserve = "single"))
              }
          }
          
          plot_title <- paste("Top", input$top_n, "Clonotype Frequencies")
          all_groups_original <- tryCatch(unique(reactive_df_raw()[[input$group_by]]), error = function(e) NULL)
          if (!is.null(all_groups_original) && length(input$filter_groups) < length(all_groups_original)) {
              plot_title <- paste(plot_title, "(Groups Filtered)")
          }

          p <- p +
            ggplot2::scale_x_discrete(limits = levels(plot_data$raw_clonotype_id)) +
            ggplot2::theme_classic(base_size = 14) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                axis.title = ggplot2::element_text(size = 12),
                plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
            ) +
            theme(legend.position = input$legend) +
            y_scale +
            ggplot2::labs(
                title = plot_title,
                x = "Raw Clonotype ID",
                y = y_label,
                fill = fill_label
            )

          if (input$horizontal) {
              p <- p + ggplot2::coord_flip() +
              ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = 10))
          }
          
          if (input$plot_type %in% c("stacked", "dodged")) {
              all_groups <- tryCatch(sort(unique(reactive_df_raw()[[grouping_var]])), error = function(e) NULL)
              if (!is.null(all_groups) && length(all_groups) > 0) {
                  color_palette <- scales::hue_pal()(length(all_groups))
                  names(color_palette) <- all_groups
                  valid_colors <- color_palette[levels(plot_data$group_col)]
                  if (length(valid_colors) > 0 && !all(is.na(valid_colors))) {
                      p <- p + ggplot2::scale_fill_manual(values = valid_colors, name = fill_label, drop = FALSE)
                  }
              }
          }
          return(p)
        })

        # 7. プロット出力
        output$plot <- renderPlot({
          plot_obj()
        },
          width = reactive(input$plot_width),
          height = reactive(input$plot_height)
        )

        # PDFダウンロードハンドラー
        output$download_plot <- downloadHandler(
          filename = function() { "clonotype_plot.pdf" },
          content = function(file) {
            ggsave(file, plot = plot_obj(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
          }
        )

        # 8. テーブル出力
        output$table <- renderDT({
          dt_data_full <- tryCatch(
            { table_for_plot() },
            error = function(e) {
              showNotification(paste("Error generating table data:", e$message), type = "error", duration = 10)
              return(NULL)
            }
          )
          
          validate(need(!is.null(dt_data_full), "No data available for the table with current selections."))
          
          dt_data <- dt_data_full %>%
            dplyr::select(-any_of("total_n"))

          if ("group_col" %in% names(dt_data)) {
            dt_data <- dt_data %>% dplyr::rename(!!input$group_by := group_col)
          }

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

        # テーブルダウンロードハンドラー
        output$download_table <- downloadHandler(
          filename = function() { "clonotype_table.csv" },
          content = function(file) {
            write.csv(table_for_plot(), file, row.names = FALSE)
          }
        )

      } else if (input$analysis_target == "gene_usage") {
            observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
      # TODO: target_gene_column の選択肢更新 observeEvent
    })



    # vdjType 変更時に Target Gene Column の選択肢を更新
    observeEvent(input$vdj_type,
      {
        req(input$vdj_type, input$vdj_type != "")

        target_choices <- list()
        selected_choice <- NULL
        # ★ 新しいV-Jペアの列名を定義 ★
        vj_pair_names <- list(
          tcr_tra = "TCR_TRA_v_j_pair",
          tcr_trb = "TCR_TRB_v_j_pair",
          bcr_igh = "BCR_IGH_v_j_pair",
          bcr_igl = "BCR_IGL_v_j_pair"
        )

        if (input$vdj_type == "tcr") {
          target_choices <- c(
            # 単一遺伝子
            "TRA V Gene" = "TCR_TRA_v_gene",
            "TRA J Gene" = "TCR_TRA_j_gene",
            "TRB V Gene" = "TCR_TRB_v_gene",
            "TRB D Gene" = "TCR_TRB_d_gene",
            "TRB J Gene" = "TCR_TRB_j_gene",
            # V-J ペア
            "TRA V-J Pair" = vj_pair_names$tcr_tra, # 定義した列名を使用
            "TRB V-J Pair" = vj_pair_names$tcr_trb # 定義した列名を使用
          )
          selected_choice <- "TCR_TRB_v_gene"
        } else if (input$vdj_type == "bcr") {
          target_choices <- c(
            # 単一遺伝子
            "IGH V Gene" = "BCR_IGH_v_gene",
            "IGH D Gene" = "BCR_IGH_d_gene",
            "IGH J Gene" = "BCR_IGH_j_gene",
            "IGH C Gene" = "BCR_IGH_c_gene",
            "IGK V Gene" = "BCR_IGK_v_gene",
            "IGK J Gene" = "BCR_IGK_j_gene",
            "IGL V Gene" = "BCR_IGL_v_gene",
            "IGL J Gene" = "BCR_IGL_j_gene",
            # V-J ペア
            "IGH V-J Pair" = vj_pair_names$bcr_igh,
            "IGK V-J Pair" = vj_pair_names$bcr_igk, # IGK追加
            "IGL V-J Pair" = vj_pair_names$bcr_igl
          )
          selected_choice <- "BCR_IGH_v_gene"
        } else {
          target_choices <- list("Select VDJ Type First" = "")
          selected_choice <- ""
        }

        # updateSelectInput で選択肢を更新
        updateSelectInput(session, "target_gene_column",
          choices = target_choices,
          selected = selected_choice
        )
      },
      ignoreNULL = FALSE,
      ignoreInit = FALSE
    )

    # 1. 元データのリアクティブ (★ V-J ペア列生成を追加 ★)
    reactive_df_raw <- reactive({
      req(input$vdj_type, input$group_by, input$target_gene_column)
      validate(need(input$target_gene_column != "", "Please select a Target Gene Column."))

      df_orig <- NULL # 元のDFを保持
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df)
        df_orig <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df)
        df_orig <- myReactives$bcr_df
      } else {
        return(NULL)
      }
      validate(need(nrow(df_orig) > 0, "Input data empty."))

      # --- V-Jペア列の生成 ---
      df_processed <- df_orig # コピーして処理

      # V-Jペア生成関数 (NA処理込み、列存在チェック強化)
      create_vj_pair <- function(df, v_col, j_col, new_col_name) {
        # 両方の列が存在するかチェック
        if (all(c(v_col, j_col) %in% names(df))) {
          df %>% dplyr::mutate(
            !!new_col_name := dplyr::if_else(
              # NAだけでなく空文字もチェック
              !is.na(.data[[v_col]]) & .data[[v_col]] != "" &
                !is.na(.data[[j_col]]) & .data[[j_col]] != "",
              paste(.data[[v_col]], .data[[j_col]], sep = "-"),
              NA_character_ # VかJがNA/空ならペアもNA
            )
          )
        } else {
          # 必要な列がない場合、警告を出してペア列は追加しない（または全てNAの列を追加）
          warning(paste(
            "Cannot create V-J pair '", new_col_name, "'. Missing required columns:",
            paste(setdiff(c(v_col, j_col), names(df)), collapse = ", ")
          ))
          df # 元のdfを返す
          # もし列がない場合でも列自体は追加したいなら:
          # df %>% dplyr::mutate(!!new_col_name := NA_character_)
        }
      }

      # observeEvent で定義した V-J ペア列名を再利用
      vj_pair_names <- list(
        tcr_tra = "TCR_TRA_v_j_pair", tcr_trb = "TCR_TRB_v_j_pair",
        bcr_igh = "BCR_IGH_v_j_pair", bcr_igk = "BCR_IGK_v_j_pair", bcr_igl = "BCR_IGL_v_j_pair"
      )

      # vdjType に応じてペア列を生成
      if (input$vdj_type == "tcr") {
        df_processed <- df_processed %>%
          create_vj_pair("TCR_TRA_v_gene", "TCR_TRA_j_gene", vj_pair_names$tcr_tra) %>%
          create_vj_pair("TCR_TRB_v_gene", "TCR_TRB_j_gene", vj_pair_names$tcr_trb)
      } else if (input$vdj_type == "bcr") {
        df_processed <- df_processed %>%
          create_vj_pair("BCR_IGH_v_gene", "BCR_IGH_j_gene", vj_pair_names$bcr_igh) %>%
          create_vj_pair("BCR_IGK_v_gene", "BCR_IGK_j_gene", vj_pair_names$bcr_igk) %>% # IGK追加
          create_vj_pair("BCR_IGL_v_gene", "BCR_IGL_j_gene", vj_pair_names$bcr_igl)
      }
      # --- V-Jペア列の生成 終了 ---

      # 選択された列 (group_by, target_gene_column) が最終的なDFに存在するか確認
      required_cols <- c(input$group_by, input$target_gene_column)
      actual_cols <- names(df_processed)
      missing_cols <- setdiff(required_cols, actual_cols)
      validate(need(
        length(missing_cols) == 0,
        paste(
          "Selected column(s) not found after processing:",
          paste(shQuote(missing_cols), collapse = ", "),
          "\nAvailable columns:", paste(actual_cols, collapse = ", ")
        )
      )) # 利用可能な列も表示

      # NA/空文字を除外 (group_by と target_gene_column)
      df_filtered_na_empty <- df_processed %>%
        dplyr::filter(
          !is.na(.data[[input$group_by]]) & .data[[input$group_by]] != "",
          !is.na(.data[[input$target_gene_column]]) & .data[[input$target_gene_column]] != ""
        )
      validate(need(
        nrow(df_filtered_na_empty) > 0,
        paste(
          "No data remaining after removing NA/empty values from columns:",
          shQuote(input$group_by), "and", shQuote(input$target_gene_column)
        )
      ))
      return(df_filtered_na_empty)
    })

    # 2. フィルタリング前の合計数 (変更なし)
    original_totals <- reactive({
      # ... (変更なし) ...
      df <- reactive_df_raw()
      req(df, input$group_by)
      grouping_var <- input$group_by
      group_tots <- df %>%
        dplyr::count(.data[[grouping_var]], name = "original_group_total") %>%
        dplyr::rename(group_col = .data[[grouping_var]]) %>%
        dplyr::mutate(group_col = as.character(group_col))
      overall_tot <- nrow(df)
      return(list(group = group_tots, overall = overall_tot))
    })

    # 3. グループフィルタリングUI (変更なし)
    output$filter_groups_ui <- renderUI({
      # ... (変更なし) ...
      df <- reactive_df_raw()
      req(df, input$group_by)
      grouping_var <- input$group_by
      available_groups <- sort(unique(df[[grouping_var]]))
      validate(need(length(available_groups) > 0, "No available groups."))
      checkboxGroupInput(session$ns("filter_groups"), paste("Filter", tools::toTitleCase(gsub("_", " ", grouping_var)), ":"),
        choices = available_groups, selected = available_groups, inline = TRUE
      )
    })

    # 4. フィルタリングされたデータフレーム (変更なし)
    filtered_df <- reactive({
      # ... (変更なし) ...
      df <- reactive_df_raw()
      req(df, input$group_by, input$filter_groups)
      grouping_var <- input$group_by
      selected_groups <- input$filter_groups
      df_filtered <- df %>% dplyr::filter(.data[[grouping_var]] %in% selected_groups)
      if (nrow(df_filtered) == 0) {
        return(NULL)
      } # エラー処理省略
      return(df_filtered)
    })

    # ★ 並び替え基準グループの observe ブロックを削除 ★
    # observe({ ... })
    # --- 5. プロット/テーブル用データの計算 (Percentage計算部分を修正) ---
    table_for_plot <- reactive({
        df_filtered <- filtered_df()
        req(df_filtered, input$plot_type %in% c("normal", "stacked", "dodged", "heatmap"),
            input$group_by, input$target_gene_column, input$target_gene_column != "",
            input$filter_groups)

        og_totals <- original_totals(); req(og_totals)
        plot_data <- NULL
        grouping_var <- input$group_by; target_col <- input$target_gene_column

        # 全てのターゲット遺伝子を自然順でソート
        all_target_genes_sorted <- df_filtered %>%
            dplyr::distinct(.data[[target_col]]) %>%
            dplyr::pull(.data[[target_col]]) %>%
            stringr::str_sort(numeric = TRUE)
        validate(need(length(all_target_genes_sorted) > 0, "No target genes found after filtering."))
        gene_levels_natural_order <- factor(all_target_genes_sorted, levels = all_target_genes_sorted)

        # 全データでカウント (group_by も含める)
        counts_by_group_raw <- df_filtered %>%
            dplyr::count(.data[[target_col]], .data[[grouping_var]], name = "n") %>%
            dplyr::rename(group_col = .data[[grouping_var]])

        # 有効なグループ (フィルター後)
        # input$group_by が選択されているか確認 (Normal Plot以外で必要)
        # Normal Plotの場合、group_col には単一の値 (e.g., "All") が入るように工夫が必要かも
        # → ここではinput$filter_groupsを使うことで、Normalでもエラーなく進む想定
        valid_group_levels <- intersect(input$filter_groups, unique(counts_by_group_raw$group_col))
        validate(need(length(valid_group_levels) > 0, "No valid groups found after filtering."))
        group_levels_to_complete <- factor(valid_group_levels, levels = valid_group_levels) # Use sorted unique groups

        # 全組み合わせを生成 (全遺伝子 x 全有効グループ)
        plot_data_completed <- counts_by_group_raw %>%
            dplyr::mutate(group_col = as.character(group_col)) %>% # Ensure character for complete
            tidyr::complete(
                !!target_col := gene_levels_natural_order,
                group_col = levels(group_levels_to_complete),
                fill = list(n = 0)
            ) %>%
            tidyr::drop_na(!!target_col)

        # 元の合計数を取得
        original_group_totals_df <- og_totals$group; req(original_group_totals_df)
        overall_total_original <- og_totals$overall; req(overall_total_original)

        # ★ Percentage計算: グループ内比率と全体比率の両方を計算 ★
        plot_data_with_perc <- plot_data_completed %>%
            dplyr::left_join(original_group_totals_df, by = "group_col") %>%
            dplyr::mutate(
                percentage_within_group = ifelse(is.na(original_group_total) | original_group_total == 0, 0, n / original_group_total),
                percentage_of_total = ifelse(overall_total_original == 0, 0, n / overall_total_original)
            ) %>%
            dplyr::select(-original_group_total)

        validate(need(nrow(plot_data_with_perc) > 0, "Data processing empty."))

        # Factor化
        plot_data <- plot_data_with_perc %>%
            dplyr::mutate(
                !!target_col := factor(.data[[target_col]], levels = levels(gene_levels_natural_order)),
                group_col = factor(group_col, levels = levels(group_levels_to_complete))
            ) %>%
            dplyr::arrange(.data[[target_col]], group_col)

        # display_value に応じて使用する percentage 列を決定 (plot_obj側で処理するため、ここでは両方保持)
        # ただし、後のテーブル表示用に 'percentage' 列はグループ内比率にしておく
        if("percentage_within_group" %in% names(plot_data)) {
             plot_data <- plot_data %>% dplyr::rename(percentage = percentage_within_group)
        }


        validate(need(!is.null(plot_data) && nrow(plot_data) > 0, "No data for plot/table."))
        return(plot_data)
    })

    # --- 6. プロットオブジェクトの生成 (★ Heatmap 対応追加 ★) ---
    plot_obj <- reactive({
        plot_data <- table_for_plot()
        req(plot_data)
        grouping_var <- input$group_by
        target_col <- input$target_gene_column

        # --- Heatmap Plot ---
        if (input$plot_type == "heatmap") {
            # Heatmap はグループ化データが前提
            validate(need("group_col" %in% names(plot_data) && length(levels(plot_data$group_col)) > 0,
                          "Heatmap requires grouped data (select a grouping variable)."))
             # Normal plot が選択された場合は Heatmap を表示しない方が親切かもしれない
             # validate(need(input$group_by != "none" || length(levels(plot_data$group_col)) > 1), "Select a grouping variable for Heatmap")


            # 表示値 (fill aesthetic)
            fill_var <- if(input$display_value == "count") "n" else "percentage" # 'percentage' はグループ内比率
            fill_label <- if(input$display_value == "count") "Count" else "Percentage (within Group)"

            # Y軸のレベル (遺伝子/ペア、自然順)
            y_levels <- levels(plot_data[[target_col]])
            # X軸のレベル (グループ)
            x_levels <- levels(plot_data$group_col)
            validate(need(length(y_levels)>0 && length(x_levels)>0, "Cannot determine axes for heatmap."))

            p <- ggplot2::ggplot(plot_data, aes(x = group_col, y = .data[[target_col]], fill = .data[[fill_var]])) +
                ggplot2::geom_tile(color = "grey85", linewidth = 0.2) + # タイル境界線調整
                ggplot2::scale_y_discrete(limits = rev(y_levels)) + # ヒートマップは通常Y軸逆順
                ggplot2::scale_x_discrete(limits = x_levels, position = "bottom") +
                # ★ 色スケール: 0が白、値が大きいほど濃い青になるように ★
                ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = fill_label,
                                             labels = if(fill_var == "percentage") scales::percent else waiver(),
                                             na.value = "grey80",
                                             guide = guide_colorbar(barwidth = 0.8, barheight = 10)) + # カラーバー調整
                ggplot2::labs(
                    title = paste("Gene Usage Frequencies Heatmap"),
                    x = tools::toTitleCase(gsub("_", " ", grouping_var)),
                    y = gsub("_", " ", tools::toTitleCase(target_col)),
                    fill = fill_label
                ) +
                ggplot2::theme_minimal(base_size = 11) + # ベースフォントサイズ調整
                ggplot2::theme(
                    axis.text.x = element_text(angle = 45, hjust = 1, size=9), # X軸ラベル調整
                    axis.text.y = element_text(size = 8),                    # Y軸ラベル調整
                    axis.ticks = element_blank(),
                    panel.grid = element_blank(),
                    legend.position = if (isTRUE(input$legend %||% TRUE)) "right" else "none",
                    plot.title = element_text(hjust = 0.5) # タイトル中央揃え
                )

             # 値表示オプション（もし追加する場合）
             # if (isTRUE(input$heatmap_show_values)) { ... }

            return(p)

        # --- Barplot (Dodged, Stacked, Normal) ---
        } else if (input$plot_type %in% c("dodged", "stacked", "normal")) {

            # y軸関連
            if (input$display_value == "count") {
                y_var <- "n"; y_label <- "Frequency (Number of Cells)"
                y_scale <- ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
            } else { # percentage
                # Normal Plot の場合は全体比率を使用、それ以外はグループ内比率
                if(input$plot_type == "normal" && "percentage_of_total" %in% names(plot_data)) {
                    y_var <- "percentage_of_total"
                    y_label <- "Frequency (% of Total Cells)"
                } else {
                    y_var <- "percentage" # グループ内比率 (rename済み)
                    y_label <- paste("Frequency (% within Original", tools::toTitleCase(gsub("_", " ", grouping_var)), ")")
                }
                # percentage 列が存在するか確認
                validate(need(y_var %in% names(plot_data), paste("Required percentage column", shQuote(y_var), "not found.")))
                y_scale <- ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.05)))
            }

            # 基本プロット
            p <- ggplot2::ggplot(plot_data)
            fill_label <- if (input$plot_type %in% c("stacked", "dodged")) { tools::toTitleCase(gsub("_", " ", grouping_var)) } else { NULL }

            # geom と aes
            if (input$plot_type == "normal") {
                 # y_var が percentage_of_total か n になる
                p <- p + ggplot2::aes(x = .data[[target_col]], y = .data[[y_var]]) +
                    ggplot2::geom_bar(stat = "identity", fill = "steelblue")
            } else { # stacked or dodged
                 # y_var が percentage (within group) か n になる
                validate(need("group_col" %in% names(plot_data), "Stacked/Dodged plot requires a grouping variable."))
                req(levels(plot_data$group_col))
                p <- p + ggplot2::aes(x = .data[[target_col]], y = .data[[y_var]], fill = group_col)
                if (input$plot_type == "stacked") {
                    p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")
                } else { # dodged
                    p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(preserve = "single"))
                }
            }

            # タイトル
            plot_title <- "Gene Usage Frequencies"
            # Check if grouping variable exists and is not 'none' or similar placeholder if you add one
            is_grouped <- "group_col" %in% names(plot_data) && length(levels(plot_data$group_col)) > 1
            all_groups_original <- tryCatch(unique(reactive_df_raw()[[input$group_by]]), error = function(e) NULL)
            # Only add "(Groups Filtered)" if it's actually grouped and filtered
            if (is_grouped && !is.null(all_groups_original) && length(input$filter_groups) < length(all_groups_original)) {
                plot_title <- paste(plot_title, "(Groups Filtered)")
            }

            # テーマ、ラベル、スケール
            p <- p +
                ggplot2::theme_classic(base_size = 14) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                    axis.title = ggplot2::element_text(size = 12),
                    plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
                    legend.position = if (isTRUE(input$legend %||% TRUE)) "right" else "none"
                ) +
                y_scale +
                ggplot2::labs( title = plot_title, x = gsub("_", " ", tools::toTitleCase(target_col)),
                               y = y_label, fill = fill_label)

            # 横向き (Barplotのみ)
            if (isTRUE(input$horizontal)) {
                 p <- p + ggplot2::coord_flip() + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = 10))
            }

            # 色 (Stacked/Dodgedのみ)
            if (input$plot_type %in% c("stacked", "dodged")) {
                # ...(色設定コードは変更なし)...
                 all_groups <- tryCatch(sort(unique(reactive_df_raw()[[grouping_var]])), error = function(e) NULL)
                 if (!is.null(all_groups) && length(all_groups) > 0) {
                     color_palette <- scales::hue_pal()(length(all_groups))
                     names(color_palette) <- all_groups
                     req(plot_data$group_col)
                     valid_colors <- color_palette[levels(plot_data$group_col)]
                     if (length(valid_colors) > 0 && !all(is.na(valid_colors))) {
                         p <- p + ggplot2::scale_fill_manual(values = valid_colors, name = fill_label, drop = FALSE)
                     }
                 }
            }
            return(p)
        } else {
            # input$plot_type が予期せぬ値の場合
            warning("Invalid plot_type selected: ", input$plot_type)
            return(ggplot() + theme_void() + labs(title = "Invalid Plot Type Selected"))
        }
    })

    # 7. プロット出力 (renderPlot, downloadHandler) (変更なし、エラーハンドリング強化)
    output$plot <- renderPlot(
      {
        # table_for_plot() が NULL を返す可能性を考慮
        plot_data_check <- tryCatch(table_for_plot(), error = function(e) NULL)
        if (is.null(plot_data_check) || nrow(plot_data_check) == 0) {
          plot.new()
          text(0.5, 0.5, "No data to display for current selections.", cex = 1.2)
          return() # プロット生成をスキップ
        }
        # プロット生成
        plot_result <- tryCatch(plot_obj(), error = function(e) {
          showNotification(paste("Error generating plot:", e$message), type = "error", duration = 10)
          ggplot() +
            theme_void() +
            geom_text(aes(0, 0, label = "Error generating plot."), size = 5)
        })
        if (!is.null(plot_result)) print(plot_result)
      },
      # commonPlotOptions からの入力を取得 (%||% でデフォルト値設定)
      width = reactive(input$plot_width %||% 600),
      height = reactive(input$plot_height %||% 500)
    )

    output$download_plot <- downloadHandler(
      filename = function() {
        # ファイル名 (変更なし)
        paste0("gene_usage_", input$target_gene_column, "_", input$plot_type, ".pdf")
      },
      content = function(file) {
        p <- plot_obj()
        req(p)
        # プロットサイズ (変更なし)
        plot_width_inch <- (input$plot_width %||% 600) / 72
        plot_height_inch <- (input$plot_height %||% 500) / 72
        ggsave(file, plot = p, width = plot_width_inch, height = plot_height_inch, device = "pdf", dpi = 300)
      }
    )

    # 8. テーブル出力 (★ DT オプション変更の可能性 ★)
    output$table <- renderDT({
      dt_data_full <- tryCatch(table_for_plot(), error = function(e) {
        showNotification(paste("Error generating table data:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      validate(need(!is.null(dt_data_full) && nrow(dt_data_full) > 0, "No data available for the table.")) # nrow > 0 を追加

      # テーブル表示用に整形 (変更なし)
      dt_data <- dt_data_full %>%
        dplyr::select(-any_of(c("total_n", "total_n_filtered", "sort_value")))
      target_col <- input$target_gene_column
      grouping_var <- input$group_by
      if ("group_col" %in% names(dt_data)) {
        dt_data <- dt_data %>% dplyr::rename(!!grouping_var := group_col)
      }
      if (input$display_value == "count") {
        dt_data <- dt_data %>% dplyr::select(-any_of("percentage"))
      } else { # percentage
        dt_data <- dt_data %>%
          dplyr::mutate(percentage = scales::percent(percentage, accuracy = 0.1)) %>%
          dplyr::select(-any_of("n"))
      }
      ordered_cols <- c(
        target_col, if (input$plot_type != "normal") grouping_var else NULL,
        if (input$display_value == "count") "n" else "percentage"
      )
      ordered_cols_exist <- intersect(ordered_cols, names(dt_data))
      dt_data <- dt_data %>% dplyr::select(all_of(ordered_cols_exist))

      # ★ DT オプション: 全要素表示のため pageLength を調整、検索有効化 ★
      datatable(dt_data,
        options = list(
          scrollX = TRUE,
          pageLength = 25, # または -1 で全件表示 (多いと重い)
          lengthMenu = c(10, 25, 50, 100, -1), # -1 = All
          searching = TRUE, # 検索を有効に
          lengthChange = TRUE # 表示件数変更を有効に
        ),
        rownames = FALSE
      )
    })

    # テーブルダウンロード (変更なし)
    output$download_table <- downloadHandler(
      filename = function() {
        # ファイル名 (変更なし)
        paste0("gene_usage_summary_", input$target_gene_column, "_", input$group_by, ".csv")
      },
      content = function(file) {
        data_to_write <- table_for_plot()
        # 列名整形 (変更なし)
        if ("group_col" %in% names(data_to_write) && input$plot_type != "normal") {
          data_to_write <- data_to_write %>% dplyr::rename(!!input$group_by := group_col)
        }
        data_to_write <- data_to_write %>% dplyr::select(-any_of(c("total_n_filtered", "sort_value")))
        write.csv(data_to_write, file, row.names = FALSE)
      }
    )
      }
    }
    
    ) # observeEvent終了
  }) # moduleServer終了
} # AnalysisPlotServer終了