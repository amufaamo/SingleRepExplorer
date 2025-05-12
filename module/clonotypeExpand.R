clonotypeExpandUI <- function(id, vdj = 'tcr'){
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns('group_by'),
                   label = 'Group',
                   choices = c('sample', 'seurat_clusters'),
                   selected = 'sample'),
      uiOutput(ns("group_value_selector")),
      radioButtons(ns('clone_call'),
                   label = "Clone call",
                   choices = list("Clonotype ID" = "raw_clonotype_id",
                                  "Sub clonotype ID" = "exact_subclonotype_id",
                                  "VDJC + CDR3 nucleotide" = "CTstrict",
                                  "CDR3 nucleotide" = "CTnt",
                                  "CDR3 amino acid" = "CTaa"),
                   selected = "raw_clonotype_id"),
      numericInput(ns('number'), label = "Number of Clonotypes", value = 20, min = 1, step = 1),
      radioButtons(ns("legend"), "Legend", choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
      sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
      sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
    ),
    mainPanel(
      plotOutput(ns("plot")),
      DTOutput(ns('table')),
      downloadButton(ns("download_data"), "Download data (.csv)"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)")
    )
  )
}


clonotypeExpandServer <- function(id, myReactives, vdj = 'tcr') {
  moduleServer(id, function(input, output, session) {

    observeEvent(myReactives$seurat_object,{
      req(myReactives$seurat_object)
      update_group_by(session, myReactives)
    })
    
    observeEvent(input$group_by, {
      req(input$group_by)
      selected_column <- input$group_by
      unique_values <- unique(data()[[selected_column]])
      
      # ユニークな値に基づいて radioButtons を動的に生成
      output$group_value_selector <- renderUI({
        radioButtons("group_value",
                     label = paste("Select", input$group_by),
                     choices = unique_values)
      })
    })
    
    # データテーブル作成部分
    table <- reactive({
      req(myReactives$seurat_object)
      group_by_column <- input$group_by  # 動的にgroup_byするカラムを指定

      df <- myReactives$seurat_object@meta.data %>%
        filter(!is.na(TCR_raw_clonotype_id)) %>%
        group_by(!!sym(group_by_column), TCR_raw_clonotype_id) %>%  # input$group_byを動的にgroup_by
        summarize(count = n(), .groups = "drop") %>%  # TCRごとのグループ内カウント
        ungroup()

      # TCRごとの合計カウントを計算し、その順序で並べる
      total_counts <- df %>%
        group_by(TCR_raw_clonotype_id) %>%
        summarize(total_count = sum(count)) %>%
        arrange(desc(total_count)) %>%
        head(input$number)

      # top20のTCR_raw_clonotype_idのみをフィルタリング
      df <- df %>% filter(TCR_raw_clonotype_id %in% total_counts$TCR_raw_clonotype_id)

      # 合計カウントに基づいてfactor化して並べ替え
      df$TCR_raw_clonotype_id <- factor(df$TCR_raw_clonotype_id, levels = total_counts$TCR_raw_clonotype_id)

      return(df)
    })

    # DTテーブルの表示
    output$table <- renderDT({
      datatable(
        table(),
        filter = 'top'
      )
    })


    # ggplot描画部分
    plot <- reactive({
      req(table())
      group_by_column <- input$group_by  # fillに使用するカラム
#      table()$TCR_raw_clonotype_id <- factor(table()$TCR_raw_clonotype_id)
      ggplot(table(), aes(x = TCR_raw_clonotype_id, y = count, fill = .data[[group_by_column]])) +
        geom_bar(stat = "identity", position = "dodge") +
#        labs(title = paste("Top 20 Clonotype Distribution by", group_by_column, "(raw_clonotype_id NA removed)"),
#             x = "TCR_raw_clonotype_id",
#             y = "Count",
#             fill = group_by_column) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = input$legend)
    })

    output$plot <- renderPlot({
      plot()
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))

  })
}