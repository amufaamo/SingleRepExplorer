#source("../utils.R")
source("utils.R")
# このモジュールを動作させるには、以下のパッケージが必要です。
# install.packages("ggvenn")
# install.packages("DT")

# --- UI ---
publicClonotypeUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # --- 共通UI ---
      vdjType(ns),
      groupByInput(ns),

      # --- このモジュール特有のUI ---
      h4("Venn Diagram Settings"),
      # ベン図で比較するグループを選択するUI (2〜4グループ推奨)
      uiOutput(ns("group_selector_ui")),
      
      # テーブルに表示する単位を選択
      radioButtons(ns("table_unit"), "Table Value Unit",
                   choices = c("Count & Proportion" = "both", "Count" = "count", "Proportion" = "proportion"),
                   selected = "both"),

      # --- 共通UI ---
      commonPlotOptions(ns)
    ),
    mainPanel(
      width = 9,
      h3("Public Clonotypes Venn Diagram"),
      # ベン図のプロットエリア
      plotOutput(ns("venn_plot")),
      
      hr(),
      
      h3("Clonotype Details"),
      # 解析したい共通部分を選択するUI
      uiOutput(ns("intersection_selector_ui")),
      # 選択した共通部分のクローン情報を表示するテーブル
      DTOutput(ns("clonotype_table"))
    )
  )
}

# --- Server ---
publicClonotypeServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    # 1. 共通のリアクティブ要素 (既存コードから流用)
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

    # 2. ★★★ このモジュール特有のロジック ★★★

    # UI: ベン図で比較するグループを選択
    output$group_selector_ui <- renderUI({
      df <- reactive_df_raw()
      req(df, input$group_by)
      available_groups <- sort(unique(df[[input$group_by]]))
      
      # デフォルトで選択するグループ数を設定 (最大4つまで)
      default_selected <- head(available_groups, 4)

      selectizeInput(session$ns("selected_groups"),
                     label = "Select groups to compare (2-4 recommended):",
                     choices = available_groups,
                     selected = default_selected,
                     multiple = TRUE,
                     options = list(maxItems = 4)) # ベン図が見やすいように最大4つに制限
    })

    # Reactive: 選択されたグループのクローンIDリストを作成
    venn_input_list <- reactive({
      req(reactive_df_raw(), input$selected_groups, length(input$selected_groups) >= 2)
      
      df_venn <- reactive_df_raw() %>%
        dplyr::filter(.data[[input$group_by]] %in% input$selected_groups)
      
      # グループごとにユニークなクローンIDのリストを作成
      split(df_venn$raw_clonotype_id, df_venn[[input$group_by]]) %>%
        lapply(unique)
    })

    # Plot: ベン図を描画
    output$venn_plot <- renderPlot({
      req(venn_input_list())
      
      ggvenn::ggvenn(
        venn_input_list(),
        fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
        stroke_size = 0.5,
        set_name_size = 4,
        text_size = 5
      ) +
      labs(title = paste("Clonotype Overlap between", tools::toTitleCase(gsub("_", " ", input$group_by)))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })

    # Reactive: 全ての共通部分(intersection)を計算
    intersections <- reactive({
        req(venn_input_list())
        sets <- venn_input_list()
        set_names <- names(sets)
        
        # 組み合わせの全パターンを生成 (例: 1, 2, 1&2)
        combinations <- unlist(lapply(1:length(sets), function(i) combn(set_names, i, simplify = FALSE)), recursive = FALSE)
        
        intersection_list <- list()
        for (combo in combinations) {
            # 共通部分のクローンIDを計算
            intersect_clones <- Reduce(intersect, sets[combo])
            
            # 他のセットに含まれないクローンIDのみを抽出
            other_sets <- set_names[!set_names %in% combo]
            if (length(other_sets) > 0) {
                clones_in_others <- unlist(sets[other_sets])
                final_clones <- setdiff(intersect_clones, clones_in_others)
            } else {
                final_clones <- intersect_clones
            }
            
            # リストに追加
            combo_name <- paste(combo, collapse = " & ")
            if(length(final_clones) > 0) {
              intersection_list[[combo_name]] <- final_clones
            }
        }
        return(intersection_list)
    })

    # UI: 共通部分を選択するドロップダウンメニュー
    output$intersection_selector_ui <- renderUI({
      req(intersections())
      
      # 共通部分のクローン数も表示
      choices_with_counts <- purrr::map2_chr(names(intersections()), intersections(), ~ paste0(.x, " (", length(.y), " clones)"))
      
      selectInput(session$ns("selected_intersection"),
                  "Select intersection to inspect:",
                  choices = choices_with_counts)
    })
    
    # Table: 選択された共通部分のクローンの詳細を表示
    output$clonotype_table <- renderDT({
      req(reactive_df_raw(), input$selected_intersection, intersections())
      
      # (X clones)という部分を取り除く
      selected_name <- sub(" \\(.*\\)$", "", input$selected_intersection)
      
      # 選択された共通部分のクローンIDを取得
      target_clonotypes <- intersections()[[selected_name]]
      
      # 元データから、選択されたグループとクローンIDでフィルタリング
      table_data <- reactive_df_raw() %>%
        dplyr::filter(
          .data[[input$group_by]] %in% input$selected_groups,
          raw_clonotype_id %in% target_clonotypes
        ) %>%
        # グループごと、クローンIDごとに細胞数をカウント
        dplyr::count(.data[[input$group_by]], raw_clonotype_id, name = "count") %>%
        # グループごとの総細胞数で割り、割合を計算
        dplyr::group_by(.data[[input$group_by]]) %>%
        dplyr::mutate(proportion = count / sum(count)) %>%
        dplyr::ungroup() %>%
        # 見やすいように整形
        dplyr::rename(group = .data[[input$group_by]]) %>%
        tidyr::pivot_wider(
            names_from = group,
            values_from = c(count, proportion),
            values_fill = 0 # 存在しない場合は0で埋める
        )

      datatable(table_data, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    })
  })
}
