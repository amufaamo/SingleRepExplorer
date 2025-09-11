#source("../utils.R")
source("utils.R")
# # UI部分
# diversityUI <- function(id, vdj = 'tcr') {
#   ns <- NS(id)
#   sidebarLayout(
#     sidebarPanel(
#       selectInput(ns('group_by'), 
#                    label = 'Group', 
#                    choices = c('sample', 'seurat_clusters'),
#                    selected = 'sample'),
#       selectInput(ns('method'),
#                    label = 'Method (default: chao1)',
#                    choices = list('chao1', 'hill', 'div', 'gini.simp', 'inv.simp', 'gini', 'raref', 'd50', 'dxx'), selected = 'chao1'),
#       radioButtons(ns('clone_call'),
#                    label = "Clone call",
#                    choices = list("CDR3 nucleotide" = "nt",
#                                   "CDR3 amino acid" = "aa"),
#                    selected = "aa"),
#       radioButtons(ns("legend"), "Legend", choices = c("right", "left", "bottom", "top", "none"), selected = "right"),
#       sliderInput(ns("plot_width"), "Width", min = 100, max = 2000, value = 500, step = 100),
#       sliderInput(ns("plot_height"), "Height", min = 100, max = 2000, value = 500, step = 100),
#       downloadButton(ns("download_plot"), "Download plot (.pdf)"),
#       downloadButton(ns("download_table"), "Download table (.csv)") # 追加
#     ),
#     mainPanel(
#       plotOutput(ns("plot")),
#       DTOutput(ns('table'))
#     )
#   )
# }

# # Server部分

# diversityServer <- function(id, myReactives, vdj = 'tcr') {
#   moduleServer(id, function(input, output, session) {
    
#     observeEvent(myReactives$seurat_object,{
#       req(myReactives$seurat_object)
#       update_group_by_select_input(session, myReactives)
#     })
    
#     # immunarchを作る
#     immdata <- reactive({
#       imm <- NULL
#       imm <- immunarch_data(myReactives, input, vdj)
#       return(imm)
#     })
    

#     table <- reactive({
#       req(immdata())
#       repDiversity(immdata()$data,
#                    .method = input$method,
#                     .col = input$clone_call)
#     })

#     plot <- reactive({
#       req(table())
#       table() %>% vis
#     })

#     output$table <- renderDT({
#       table()
#     })


#     # プロットのレンダリング
#     output$plot <- renderPlot({
#       plot()
#     }, width = reactive(input$plot_width), height = reactive(input$plot_height))

#     # PDFダウンロードハンドラー
#     output$download_plot <- downloadHandler(
#       filename = function() { "UMAP_plot.pdf" },
#       content = function(file) {
#         ggsave(file, plot = renderCustomPlot(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
#       }
#     )


#     # CSVダウンロードハンドラー
#     output$download_table <- downloadHandler(
#       filename = function() { "table.csv" },
#       content = function(file) {
#         write.csv(table(), file)  # table()で生成したデータをCSVとして保存
#       }
#     )

    
#   })
# }



# immunarch_data <- function(myReactives, input, vdj){
#   unlink('immunarch', recursive = TRUE)
#   dir.create("immunarch", showWarnings = FALSE)
#   if(vdj == 'tcr'){
#     req(myReactives$tcr_path)
#     dt <- read.csv(myReactives$tcr_path)
#   } else if (vdj == 'bcr'){
#     req(myReactives$bcr_path)
#     dt <- read.csv(myReactives$bcr_path)
#   }
#   dt <- dt %>% mutate(sample = str_remove(barcode, "^.+-"))
#   split_column <- input$group_by
#   data_list <- split(dt, dt[[split_column]])
#   for (name in names(data_list)) {
#     write.csv(data_list[[name]], file = paste0("immunarch/", name, ".csv"), row.names = FALSE, quote = FALSE)
#   }
#   imm <- repLoad("immunarch/")
#   unlink('immunarch')
#   return(imm)
# }


# library(shiny)
# library(shinyjs)
# library(immunarch)
# library(dplyr)
# library(stringr)
# library(ggplot2)
# library(DT)
# library(Seurat) # update_group_by_select_inputのために必要かもしれません

# --- UI部分 ---
  ns <- NS(id)
diversityUI <- function(id) { # vdj引数を削除
  sidebarLayout(
    sidebarPanel(
      # VDJタイプを選択するselectInputを追加
      vdjType(ns),
      selectInput(ns('clone_call'),
                   label = "Clone call",
                   choices = list("CDR3 nucleotide" = "nt",
                                  "CDR3 amino acid" = "aa"),
                   selected = "aa"),
      selectInput(ns('group_by'),
                   label = 'Group',
                   choices = c('sample', 'seurat_clusters'), # 初期値。Server側で更新される
                   selected = 'sample'),
      selectInput(ns('method'),
                   label = 'Method (default: chao1)',
                   choices = list('chao1', 'hill', 'div', 'gini.simp', 'inv.simp', 'gini', 'raref', 'd50', 'dxx'), selected = 'chao1'),
                   commonPlotOptions(ns),
    ),
    mainPanel(
      downloadButton(ns("download_plot"), "Download plot (.pdf)"),
      plotOutput(ns("plot")),
      downloadButton(ns("download_table"), "Download table (.csv)"),
      DTOutput(ns('table')),
    )
  )
}

# --- Server部分 ---
diversityServer <- function(id, myReactives) { # vdj引数を削除
  moduleServer(id, function(input, output, session) {

    # update group_by
    observeEvent(myReactives$seurat_object,{
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    # immunarchデータを作成
    immdata <- reactive({
      # vdj_typeが選択されるまで待機
      req(input$vdj_type)
      # immunarch_dataに関数にinput$vdj_typeを渡す
      imm <- immunarch_data(myReactives, input, input$vdj_type)
      return(imm)
    })

    table <- reactive({
      req(immdata())
      # immunarch_dataがNULLでないことを確認
      validate(need(!is.null(immdata()$data), "immunarch data could not be loaded. Check input files and VDJ type selection."))
      repDiversity(immdata()$data,
                   .method = input$method,
                    .col = input$clone_call)
    })

    plot_obj <- reactive({ # plotオブジェクトを返すreactiveに変更
      req(table())
      # tableが空でないことを確認
      validate(need(nrow(table()) > 0, "Diversity table is empty."))
      vis(table()) + theme(legend.position = input$legend) # 凡例位置を適用
    })

    output$table <- renderDT({
      req(table())
      datatable(table(), options = list(scrollX = TRUE))
    })


    # プロットのレンダリング
    output$plot <- renderPlot({
      plot_obj() # plot_obj() を呼び出す
    }, width = reactive(input$plot_width), height = reactive(input$plot_height))

    # PDFダウンロードハンドラー
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("diversity_plot_", input$vdj_type, "_", input$group_by, "_", input$method, ".pdf")
      },
      content = function(file) {
        # plot_obj() を使用
        ggsave(file, plot = plot_obj(), width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300, device = "pdf")
      }
    )

    # CSVダウンロードハンドラー
    output$download_table <- downloadHandler(
      filename = function() {
        paste0("diversity_table_", input$vdj_type, "_", input$group_by, "_", input$method, ".csv")
      },
      content = function(file) {
        write.csv(table(), file, row.names = FALSE) # row.names=FALSE を推奨
      }
    )

  })
}


# --- データ読み込み・整形関数 ---
immunarch_data <- function(myReactives, input, vdj){
  # 一時ディレクトリのクリーンアップと作成
  unlink('immunarch_temp_dir', recursive = TRUE, force = TRUE)
  dir.create("immunarch_temp_dir", showWarnings = FALSE)

  dt <- NULL
  # 選択されたVDJタイプに基づいてデータを読み込む
  if(vdj == 'tcr'){
    req(myReactives$tcr_path) # 必要なファイルパスが存在するか確認
    # ファイルが存在するか確認
    validate(need(file.exists(myReactives$tcr_path), paste("TCR data file not found:", myReactives$tcr_path)))
    tryCatch({
        dt <- read.csv(myReactives$tcr_path)
    }, error = function(e) {
        showNotification(paste("Error reading TCR data:", e$message), type = "error")
        return(NULL) # エラー時はNULLを返す
    })
  } else if (vdj == 'bcr'){
    req(myReactives$bcr_path) # 必要なファイルパスが存在するか確認
    # ファイルが存在するか確認
    validate(need(file.exists(myReactives$bcr_path), paste("BCR data file not found:", myReactives$bcr_path)))
     tryCatch({
        dt <- read.csv(myReactives$bcr_path)
    }, error = function(e) {
        showNotification(paste("Error reading BCR data:", e$message), type = "error")
        return(NULL) # エラー時はNULLを返す
    })
  }

  # dtが正しく読み込めたか確認
  validate(need(!is.null(dt) && nrow(dt) > 0, "Input data is empty or could not be loaded."))

  # Group Byカラムが存在するか確認
  split_column <- input$group_by
  validate(need(split_column %in% colnames(dt), paste("Grouping column '", split_column, "' not found in the data.")))

  # sample列がない場合はbarcodeから作成 (エラーハンドリング強化)
  if (!"sample" %in% colnames(dt) && "barcode" %in% colnames(dt)) {
      if(any(grepl("^.+-", dt$barcode))){ # ハイフンが含まれるかチェック
         dt <- dt %>% mutate(sample = str_remove(barcode, "^.+-"))
      } else {
          # ハイフンがない場合の代替処理（例：barcode全体をsampleとする）
          # あるいはエラー通知
          showNotification("Barcode format unexpected for deriving 'sample'. Using full barcode.", type = "warning")
          dt <- dt %>% mutate(sample = barcode)
      }
  } else if (!"sample" %in% colnames(dt)) {
      validate(need(FALSE, "Column 'sample' or 'barcode' (for deriving sample) not found."))
  }


  imm <- NULL # immを初期化
  tryCatch({
    # データをグループごとに分割してファイルに書き出す
    data_list <- split(dt, dt[[split_column]])
    # ファイル名に使えない文字を置換する関数
    sanitize_filename <- function(name) {
      gsub("[^a-zA-Z0-9_.-]", "_", name)
    }
    for (name in names(data_list)) {
      sanitized_name <- sanitize_filename(name)
      write.csv(data_list[[name]], file = file.path("immunarch_temp_dir", paste0(sanitized_name, ".csv")), row.names = FALSE, quote = FALSE)
    }

    # immunarchで読み込む
    imm <- repLoad("immunarch_temp_dir/")

  }, error = function(e) {
      showNotification(paste("Error processing data for immunarch:", e$message), type = "error")
      imm <- NULL # エラーが発生したらimmをNULLにする
  }, finally = {
    # 一時ディレクトリを削除
    unlink('immunarch_temp_dir', recursive = TRUE, force = TRUE)
  })

  # immがNULLでないか、またはデータが含まれているかを確認
  if(is.null(imm) || length(imm$data) == 0){
      showNotification("Failed to load data into immunarch.", type = "error")
      return(NULL) # 失敗したらNULLを返す
  }

  return(imm)
}