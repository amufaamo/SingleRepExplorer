# --- UI Definition ---
cdrLengthDistUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      vdjType(ns),
      selectInput(ns("target_cdr_column"), "Target Column (for Length Calculation)",
        choices = c("Select VDJ type first..." = ""), # Server側で設定
        selected = ""
      ),
      groupByInput(ns),
      uiOutput(ns("filter_groups_ui")),

      selectInput(ns("plot_type"), "Plot Type",
        choices = c(
          "Dodged Bars" = "dodge", # グループ比較、横並び
          "Stacked Bars" = "stack", # グループ比較、積み上げ
          "Boxplot" = "boxplot"
        ),
        selected = "dodge"
      ),
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] != 'boxplot'"),
        selectInput(ns("display_value"), "Display Value",
          choices = c("Count" = "count", "Percentage (%)" = "percentage"),
          selected = "count"
        )
      ),
      uiOutput(ns("filter_groups_ui")), # グループフィルターUI (動的に生成)
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] != 'boxplot'"),
        checkboxInput(ns("horizontal"), "Horizontal Barplot", value = FALSE)
      ),
      hr(),
      commonPlotOptions(ns) # 共通プロットオプション
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download Plot (.pdf)"),
      plotOutput(ns('plot')),
      hr(),
      h3("Table"),
      downloadButton(ns("download_table"), "Download Table (.csv)"),
      DTOutput(ns("table")),
    )
  )
}

# --- Server Definition ---
cdrLengthDistServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    observeEvent(myReactives$seurat_object, {
      update_group_by_select_input(session, myReactives)
    })

    observeEvent(input$vdj_type, {
      if (input$vdj_type == "tcr") {
        req(myReactives$tcr_df) # 対応するデータフレームが必要
        target_choices_all <- c(
          "TCR Paired CTnt" = "TCR_pair_CTnt", "TCR Paired CTaa" = "TCR_pair_CTaa",
          "TCR TRA CDR3 aa" = "TCR_TRA_cdr3", "TCR TRA CDR3 nt" = "TCR_TRA_cdr3_nt",
          "TCR TRB CDR3 aa" = "TCR_TRB_cdr3", "TCR TRB CDR3 nt" = "TCR_TRB_cdr3_nt"
        )
        selected_choice <- "TCR_TRB_cdr3"
      } else if (input$vdj_type == "bcr") {
        req(myReactives$bcr_df) # 対応するデータフレームが必要
        target_choices_all <- c(
          "BCR Paired CTnt" = "BCR_pair_CTnt", "BCR Paired CTaa" = "BCR_pair_CTaa",
          "BCR IGH CDR3 aa" = "BCR_IGH_cdr3", "BCR IGH CDR3 nt" = "BCR_IGH_cdr3_nt",
          "BCR IGK CDR3 aa" = "BCR_IGK_cdr3", "BCR IGK CDR3 nt" = "BCR_IGK_cdr3_nt", # IGKも追加しておく
          "BCR IGL CDR3 aa" = "BCR_IGL_cdr3", "BCR IGL CDR3 nt" = "BCR_IGL_cdr3_nt"
        )
        selected_choice <- "BCR_IGH_cdr3"
      }
      updateSelectInput(session, "target_cdr_column", choices = target_choices_all, selected = selected_choice)
    })

    # 2. Update Target Column Choices when VDJ Type changes
    reactive_df_raw <- reactive({
      if (input$vdj_type == "tcr" && !is.null(myReactives$tcr_df)) {
        df <- myReactives$tcr_df
      } else if (input$vdj_type == "bcr" && !is.null(myReactives$bcr_df)) {
        df <- myReactives$bcr_df
      }
      return(df)
    })

    table <- reactive({
      req(reactive_df_raw())
      target_col_sym <- sym(input$target_cdr_column)

      df <- reactive_df_raw() %>%
        mutate(
          length = nchar(!!target_col_sym)
        )
      df <- df %>% dplyr::select(input$group_by, length)
      df <- df %>%
        group_by(.data[[input$group_by]], length) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(.data[[input$group_by]]) %>%
        mutate(percentage = count / sum(count)) %>%
        ungroup()
      return(df)
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
    table_filtered <- reactive({
      df <- table()
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

    #boxplot用のデータフレーム
    boxplot_df <- reactive({
      req(reactive_df_raw())
      target_col_sym <- sym(input$target_cdr_column)

      df <- reactive_df_raw() %>%
        mutate(
          length = nchar(!!target_col_sym)
        )
      df <- df %>% dplyr::select(input$group_by, length)
      
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

    output$plot <- renderPlot({
      req(table_filtered())
      if (input$plot_type == "dodge" || input$plot_type == "stack") {
        if (input$display_value == "count") {
          p <-  table_filtered() %>% ggplot(aes(x = length, y = count, fill = .data[[input$group_by]]))
        } else if (input$display_value == "percentage"){
          p <-  table_filtered() %>% ggplot(aes(x = length, y = percentage, fill = .data[[input$group_by]]))
        }
        p <- p + geom_bar(stat = 'identity', position = input$plot_type) + 
          labs(x = "Length", y = input$display_value, fill = input$group_by) +
          theme_classic() + 
          scale_y_continuous(expand = c( 0, 0 ))

        # 横向きプロットの適用
        if (input$horizontal) {
          p <- p + ggplot2::coord_flip() +
            ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = 10)) # 横向き用にy軸テキスト調整
        }
      } else if (input$plot_type == "boxplot"){
        req(boxplot_df())
        p <- boxplot_df() %>% ggplot(aes(x = .data[[input$group_by]], y = length, fill = .data[[input$group_by]])) +
          geom_boxplot() +
          labs(x = input$group_by, y = "Length", fill = input$group_by) +
          theme_classic()
      }
      return(p)
    })

    output$table <- renderDT({
      table_filtered()
    })
  }) # moduleServer end
} # cdrLengthDistServer end
