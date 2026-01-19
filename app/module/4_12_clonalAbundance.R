#source("../utils.R")
source("utils.R")
# --- UI ---
# UI部分は、軸スケールのラベルをより分かりやすく変更した以外は、元のままでOKです！
clonalAbundancePlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      # --- 流用するUI ---
      vdjType(ns),
      groupByInput(ns),

      # --- このモジュール特有のUI ---
      h4("Plot Settings"),
      # ★★★ Y軸の単位を「割合」か「細胞数」か選択するように変更 ★★★
      selectInput(ns("yaxis_unit"), "Y-Axis Unit",
                  choices = c("Percentage" = "proportion", "Count" = "count"),
                  selected = "proportion"),

      # グループ要素フィルタリング用UI
      uiOutput(ns("filter_groups_ui")),

      # --- 流用するUI ---
      commonPlotOptions(ns)
    ),
    mainPanel(
      # --- 流用するUI ---
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      h3("Table"),
      downloadButton(ns("download_table"), "Download table (.csv)"),
      DTOutput(ns("table"))
    )
  )
}

# --- Server ---
clonalAbundancePlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    # 1. 共通のリアクティブ要素 (この部分は変更ありません)
    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    reactive_df_raw <- reactive({
      req(input$vdj_type, input$group_by)
      df <- if (input$vdj_type == "tcr") myReactives$tcr_df else myReactives$bcr_df
      req(df, "raw_clonotype_id" %in% names(df), input$group_by %in% names(df))
      df_filtered <- df %>%
        dplyr::filter(!is.na(.data[[input$group_by]]) & .data[[input$group_by]] != "")
      validate(need(nrow(df_filtered) > 0, paste("No data after removing NA/empty from", input$group_by)))
      df_filtered
    })

    output$filter_groups_ui <- renderUI({
      df <- reactive_df_raw()
      req(df, input$group_by)
      available_groups <- sort(unique(df[[input$group_by]]))
      validate(need(length(available_groups) > 0, "No available groups found."))
      checkboxGroupInput(session$ns("filter_groups"),
                         label = paste("Filter", tools::toTitleCase(gsub("_", " ", input$group_by)), ":"),
                         choices = available_groups,
                         selected = available_groups,
                         inline = TRUE)
    })

    filtered_df <- reactive({
      df <- reactive_df_raw()
      req(df, input$filter_groups)
      df %>% dplyr::filter(.data[[input$group_by]] %in% input$filter_groups)
    })

    # 2. データ処理ロジック (この部分は変更ありません)
    # クローンの細胞数を数え、割合とランキングを計算します
    table_for_plot <- reactive({
      df <- filtered_df()
      req(df, nrow(df) > 0)
      grouping_var <- input$group_by

      # 1. グループごと、クローンIDごとに細胞数をカウント
      plot_data <- df %>%
        dplyr::count(.data[[grouping_var]], raw_clonotype_id, name = "count") %>%
        # 2. グループ内で、各クローンの割合(proportion)を計算
        dplyr::group_by(.data[[grouping_var]]) %>%
        dplyr::mutate(proportion = count / sum(count)) %>%
        # 3. 割合の大きい順に並び替え、ランキングを付与
        dplyr::arrange(desc(proportion), .by_group = TRUE) %>%
        dplyr::mutate(rank = row_number()) %>%
        dplyr::ungroup() %>%
        # 4. プロットしやすいようにグループ列の名前を変更
        dplyr::rename(group = .data[[grouping_var]])

      validate(need(nrow(plot_data) > 0, "No data to plot after processing."))
      return(plot_data)
    })

    # 3. ★★★ プロットオブジェクトを更新 ★★★
    # Y軸の単位をUIの選択に応じて動的に変更します
    plot_obj <- reactive({
      plot_data <- table_for_plot()
      req(plot_data, input$yaxis_unit)

      # Y軸のラベルを動的に設定
      y_label <- if (input$yaxis_unit == "proportion") "Abundance (Proportion)" else "Abundance (Count)"

      # プロットの基本設定
      p <- ggplot(plot_data, aes(x = rank, y = .data[[input$yaxis_unit]], color = group)) +
        geom_line(linewidth = 1) +
        # X軸は対数スケールで固定（ランキングプロットで一般的）
        scale_x_log10() +
        # ラベルとテーマ
        labs(
          title = "Clonal Abundance",
          x = "Rank",
          y = y_label,
          color = tools::toTitleCase(gsub("_", " ", input$group_by))
        ) +
        theme_classic(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = input$legend
        )

      # Y軸が割合(proportion)の時だけ、パーセント表示にフォーマット
      if (input$yaxis_unit == "proportion") {
        p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
      }

      return(p)
    })

    # 4. 出力 (この部分は変更ありません)
    output$plot <- renderPlot({
      plot_obj()
    },
      width = reactive(input$plot_width),
      height = reactive(input$plot_height)
    )

    output$download_plot <- downloadHandler(
        filename = function() {"clonal_abundance_rank_plot.pdf"},
        content = function(file) {
            ggsave(file, plot = plot_obj(), width = input$plot_width/72, height = input$plot_height/72, device = "pdf")
        }
    )

    output$table <- renderDT({
        table_display <- table_for_plot() %>%
            dplyr::select(
                !!input$group_by := group,
                raw_clonotype_id,
                rank,
                count,
                proportion
            )
        datatable(table_display,
                  options = list(scrollX = TRUE, pageLength = 10),
                  colnames = c(tools::toTitleCase(gsub("_", " ", input$group_by)), "Clonotype ID", "Rank", "Count", "Proportion")
        )
    })

    output$download_table <- downloadHandler(
        filename = function() {"clonal_abundance_rank_table.csv"},
        content = function(file) {
            table_display <- table_for_plot() %>%
                dplyr::rename(!!input$group_by := group)
            write.csv(table_display, file, row.names = FALSE)
        }
    )
  })
}
