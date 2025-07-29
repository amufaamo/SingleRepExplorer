# 必要なライブラリがapp.Rで読み込まれていることを前提としています
# library(shiny)
# library(Seurat)
# library(ggplot2)
# library(dplyr)
# library(DT)
# library(shinyjs)
# library(shinyalert)


# UI部分
reductionPlotUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      reductionInput(ns),
      div(style = "display: flex; align-items: center; gap: 10px;",
        div(style = "flex-grow: 1;", groupByInput(ns)),
        div(style = "flex-shrink: 0;", actionButton(ns("rename_btn"), "Rename", icon = icon("pen"), style = "height: 100%;"))
      ),
      checkboxInput(ns("highlight_toggle"), "Highlight Cells", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s'] == true", ns("highlight_toggle")),
        checkboxGroupInput(ns("unique_group"), "Select Group", choices = NULL, inline = TRUE)
      ),
      pointSizeInput(ns),
      labelSizeInput(ns),
      commonPlotOptions(ns),
    ),
    mainPanel(
      h3("Plot"),
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      hr(),
      h3("Coordinates Table"),
      downloadButton(ns("download_table"), "Download Table (.csv)"),
      DT::dataTableOutput(ns("coordinates_table"))
    )
  )
}

# Server部分
reductionPlotServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      
      so <- myReactives$seurat_object
      
      # --- 修正点：TCR/BCR列の存在を確認し、なければ追加するロジックを改善 ---
      # このロジックはオブジェクトが最初にロードされたときなど、列が存在しない場合のみ実行
      if (!"TCR_BCR_status" %in% names(so@meta.data)) {
        # scRepertoireやimmunarch由来の一般的な列名をリストアップ
        clonotype_cols <- c("CTaa", "CTnt", "raw_clonotype_id", "clonotype_id", "clonotype") 
        
        if (any(clonotype_cols %in% names(so@meta.data))) {
          # メタデータに存在する最初のクロノタイプ列を特定
          col_to_use <- clonotype_cols[clonotype_cols %in% names(so@meta.data)][1]
          
          # 新しい列 'TCR_BCR_status' を作成
          so$TCR_BCR_status <- ifelse(
            is.na(so@meta.data[[col_to_use]]) | so@meta.data[[col_to_use]] == "None",
            "Absent",
            "Present"
          )
          # プロットの順序や色のためにfactor型に変換
          so$TCR_BCR_status <- factor(so$TCR_BCR_status, levels = c("Present", "Absent"))
          
          # 更新したSeuratオブジェクトをmyReactivesに戻す
          myReactives$seurat_object <- so
        }
      }
      # --- ここまでが修正点 ---

      # UIの選択肢を更新
      # myReactives$seurat_objectが更新された場合も、されなかった場合も、
      # 現在のオブジェクトの状態に基づいてUIを更新する
      update_reduction_choices(session, myReactives)
      update_group_by_choices(session, myReactives) # この中でデフォルト値が設定される
      
      # input$group_byの変更は下のobserveEventが検知して
      # update_unique_group_choicesを呼び出すので、ここでの呼び出しは不要
    }, ignoreNULL = TRUE) # NULLを無視するオプション

    observeEvent(input$group_by, {
      req(myReactives$seurat_object, input$group_by)
      update_unique_group_choices(session, myReactives, input$group_by)
      updateCheckboxInput(session, "highlight_toggle", value = FALSE)
    })
    
    observeEvent(input$rename_btn, {
      req(myReactives$seurat_object, input$group_by)
      so <- myReactives$seurat_object
      group_by_col <- input$group_by
      
      current_levels <- if (is.factor(so[[group_by_col, drop = TRUE]])) {
        levels(so[[group_by_col, drop = TRUE]])
      } else {
        sort(unique(so[[group_by_col, drop = TRUE]]))
      }
      
      showModal(modalDialog(
        title = "Rename/Regroup Clusters",
        textInput(ns("new_group_name"), "New Group Name", placeholder = "e.g., Celltype_manual"),
        hr(),
        lapply(current_levels, function(level) {
          textInput(ns(paste0("rename_", level)), label = paste("New name for cluster:", level), value = level)
        }),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("save_rename"), "Save")
        )
      ))
    })

    observeEvent(input$save_rename, {
      req(input$new_group_name, myReactives$seurat_object)
      
      if (grepl("\\s|[^A-Za-z0-9_.]", input$new_group_name)) {
          if(require(shinyalert)) shinyalert("Error", "Group name cannot contain spaces or special characters.", type = "error")
          return()
      }

      so <- myReactives$seurat_object
      group_by_col <- input$group_by

      current_levels <- if (is.factor(so[[group_by_col, drop = TRUE]])) {
        levels(so[[group_by_col, drop = TRUE]])
      } else {
        sort(unique(so[[group_by_col, drop = TRUE]]))
      }
      
      new_names_map <- setNames(
        sapply(current_levels, function(level) input[[paste0("rename_", level)]]),
        current_levels
      )
      
      original_vector <- so[[group_by_col, drop = TRUE]]
      new_group_vector <- factor(unname(new_names_map[as.character(original_vector)]))
      
      so@meta.data[[input$new_group_name]] <- new_group_vector

      myReactives$seurat_object <- so
      
      update_group_by_choices(session, myReactives, selected = input$new_group_name)
      
      removeModal()
    })

    plot <- reactive({
      req(myReactives$seurat_object, input$reduction, input$group_by)
      so <- myReactives$seurat_object
      
      if (input$highlight_toggle && !is.null(input$unique_group)) {
        p <- DimPlot(
          so,
          reduction = input$reduction,
          group.by = input$group_by,
          cells.highlight = WhichCells(so, expression = !!sym(input$group_by) %in% !!input$unique_group),
          cols.highlight = 'darkblue',
          cols = 'grey',
          pt.size = input$point_size,
          label.size = input$label_size,
          label = TRUE
        ) + theme(legend.position = input$legend)
      } else {
        p <- DimPlot(
          so,
          reduction = input$reduction,
          group.by = input$group_by,
          label = TRUE,
          pt.size = input$point_size,
          label.size = input$label_size
        ) + theme(legend.position = input$legend)
      }
      return(p)
    })

    output$plot <- renderPlot({
      plot()
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))

    coordinates_data <- reactive({
        req(myReactives$seurat_object, input$reduction, input$group_by)
        so <- myReactives$seurat_object
        reduction_data <- so@reductions[[input$reduction]]@cell.embeddings
        table_data <- cbind(so@meta.data, reduction_data)
        reduction_prefix <- paste0(input$reduction, "_")
        table_data %>% select(all_of(input$group_by), starts_with(reduction_prefix))
    })

    output$coordinates_table <- DT::renderDataTable({
      coordinates_data()
    })

    output$download_plot <- downloadHandler(
      filename = function() { paste0("DimPlot_", input$reduction, "_", input$group_by, ".pdf") },
      content = function(file) {
        ggsave(file, plot = plot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
      }
    )

    output$download_table <- downloadHandler(
      filename = function() {
        paste0("coordinates_", input$reduction, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(coordinates_data(), file, row.names = FALSE)
      }
    )
  })
}

# --- ヘルパー関数群（このモジュール内で使用） ---

update_reduction_choices <- function(session, myReactives) {
  req(myReactives$seurat_object)
  reduction_names <- names(myReactives$seurat_object@reductions)
  choices <- stats::setNames(reduction_names, toupper(reduction_names))
  default_selection <- if ("umap" %in% reduction_names) "umap" else reduction_names[1]
  updateSelectInput(session, "reduction", choices = choices, selected = default_selection)
}

update_group_by_choices <- function(session, myReactives, selected = NULL) {
  req(myReactives$seurat_object)
  so <- myReactives$seurat_object
  
  all_choices <- names(so@meta.data)[sapply(so@meta.data, function(x) is.factor(x) || is.character(x))]
  
  choices_to_exclude <- c("orig.ident")
  filtered_choices <- all_choices[!all_choices %in% choices_to_exclude]
  filtered_choices <- filtered_choices[!grepl("_snn_res\\.", filtered_choices)]
  
  if (is.null(selected) || !selected %in% filtered_choices) {
    if ("TCR_BCR_status" %in% filtered_choices) {
      selected <- "TCR_BCR_status"
    } else if ("seurat_clusters" %in% filtered_choices) {
      selected <- "seurat_clusters"
    } else {
      selected <- filtered_choices[1]
    }
  }
  
  updateSelectInput(session, "group_by", choices = filtered_choices, selected = selected)
}

update_unique_group_choices <- function(session, myReactives, group_by) {
  req(myReactives$seurat_object, group_by)
  so <- myReactives$seurat_object
  
  if (group_by %in% names(so@meta.data)) {
    unique_values <- if (is.factor(so@meta.data[[group_by]])) {
      levels(so@meta.data[[group_by]])
    } else {
      sort(unique(so@meta.data[[group_by]]))
    }
    updateCheckboxGroupInput(session, "unique_group", choices = unique_values, inline = TRUE)
  }
}

reductionInput <- function(ns) {
  selectInput(ns("reduction"), "Reduction", choices = NULL)
}

groupByInput <- function(ns) {
  selectInput(ns("group_by"), "Group by", choices = NULL)
}

pointSizeInput <- function(ns) {
  sliderInput(ns("point_size"), "Size of points", min = 0.01, max = 10, value = 0.1, step = 0.01)
}

labelSizeInput <- function(ns) {
  sliderInput(ns("label_size"), "Size of labels", min = 0, max = 20, value = 10, step = 1)
}

commonPlotOptions <- function(ns) {
  tagList(
    selectInput(ns("legend"), "Legend Position",
      choices = c("right", "left", "top", "bottom", "none"),
      selected = "right"
    ),
    sliderInput(ns("plot_width"), "Plot Width (px)", min = 200, max = 2000, value = 600),
    sliderInput(ns("plot_height"), "Plot Height (px)", min = 200, max = 2000, value = 500)
  )
}
