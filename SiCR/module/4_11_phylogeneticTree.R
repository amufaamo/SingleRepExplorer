# 必要なパッケージをロード
library(shiny)
library(ggtree) # 系統樹プロット用だが、現在のコードでは未使用
library(dowser) # 系統樹構築用だが、現在のコードでは未使用
library(dplyr)
library(readr)
library(ggplot2)
library(tibble) # req()のため推奨

# --- Shiny UI 定義 ---
phylogeneticTreeUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      selectInput(ns("selected_clone"), "Select Clone:", # ラベルを少し変更
        choices = NULL, # Server側で設定
        selected = NULL
      ),
      # ボタンの役割に合わせてラベル変更も検討 (例: "Load & Display Data")
      actionButton(ns("run_phylotree"), "Display Clone Data", icon = icon("table"), class = "btn-primary"),
      # 必要であれば系統樹用の設定やボタンをここに追加
    ),
    mainPanel(
      # plotOutputは系統樹を表示する場合に使う。高さは固定値に変更。
      plotOutput(ns("plot"), width = "100%", height = "600px"),
      tags$h4("Selected Clone Data"), # テーブルのタイトル追加
      tableOutput(ns('table')) # ns() を使用
    )
  )
}

# --- Shiny Server 定義 ---
phylogeneticTreeServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {

    # クローン選択肢の更新ロジック (エラーハンドリング追加)
    observe({
      req(myReactives$bcr_df)
      df <- myReactives$bcr_df

      # 必要な列の存在チェック
      required_cols <- c("raw_clonotype_id", "exact_subclonotype_id")
      if (!all(required_cols %in% names(df))) {
        showNotification("Required columns for clone selection (raw_clonotype_id, exact_subclonotype_id) not found.", type = "error", duration = 10)
        updateSelectInput(session, "selected_clone", choices = c("Error: Missing columns" = ""), selected = "")
        return()
      }

      clones_summary <- tryCatch({
        df %>%
          dplyr::select(all_of(required_cols)) %>%
          dplyr::filter(!is.na(raw_clonotype_id) & raw_clonotype_id != "",
                        !is.na(exact_subclonotype_id) & exact_subclonotype_id != "") %>%
          dplyr::group_by(raw_clonotype_id) %>%
          dplyr::summarise(n_distinct_subclones = dplyr::n_distinct(exact_subclonotype_id), .groups = 'drop') %>%
          dplyr::filter(n_distinct_subclones > 1) %>%
          dplyr::pull(raw_clonotype_id) %>%
          sort()
      }, error = function(e){
          showNotification(paste("Error creating clone list:", e$message), type="error")
          return(character(0))
      })

      if (length(clones_summary) > 0) {
          updateSelectInput(session, "selected_clone", choices = clones_summary, selected = clones_summary[1])
      } else {
          updateSelectInput(session, "selected_clone", choices = c("No clones with >1 subclones found" = ""), selected = "")
          # Avoid notification flood if bcr_df is initially NULL
          if(!is.null(isolate(myReactives$bcr_df))) {
              showNotification("No suitable clones found.", type = "warning")
          }
      }
    })
 # 2. 選択されたクローンのデータを抽出 (ボタンクリックでトリガー)
    selected_clone_data <- eventReactive(input$run_phylotree, {
      req(myReactives$bcr_df, input$selected_clone, input$selected_clone != "")
      selected_clonotype_id <- input$selected_clone
      df <- myReactives$bcr_df

      # 系統樹構築に必要な基本列
      required_cols_tree <- c("barcode", "raw_clonotype_id", "BCR_IGH_full_length_nt")
      if (!all(required_cols_tree %in% names(df))) {
        showNotification(paste("Missing columns for tree building:", paste(setdiff(required_cols_tree, names(df)), collapse=", ")), type = "error")
        return(NULL) # 処理中断
      }

      clone_subset <- df %>%
        dplyr::filter(
          raw_clonotype_id == selected_clonotype_id,
          !is.na(BCR_IGH_full_length_nt) & BCR_IGH_full_length_nt != "", # 有効なシーケンス
          !is.na(barcode) & barcode != ""                        # 有効なID
        ) %>%
        dplyr::distinct(barcode, .keep_all = TRUE) # barcodeでユニークにする

      # シーケンス数が2未満の場合は系統樹を構築できない
      if (nrow(clone_subset) < 2) {
        showNotification(paste("Not enough valid sequences (< 2) for clone '", selected_clonotype_id, "' to build a tree."), type = "warning")
        # テーブル表示用にデータは返すかもしれないが、系統樹は作れない
        return(clone_subset) # テーブル表示は可能にする
      }

      showNotification(paste("Data loaded for clone:", selected_clonotype_id, "- Found", nrow(clone_subset), "unique sequences for tree."), type="message")
      return(clone_subset)
    }, ignoreNULL = FALSE) # input$selected_cloneが初期NULLでも反応開始（reqで止まる）

    # --- 系統樹構築ステップ ---

    # 3. dowser 用のデータフォーマット
    phylotree_input_data <- reactive({
      # まずデータが抽出され、有効であることを確認
      req(selected_clone_data())
      clone_subset <- selected_clone_data()
      # 系統樹構築可能なシーケンス数か確認
      req(nrow(clone_subset) >= 2)

      # メタデータ列の決定 (UIでの選択を反映)
      potential_meta_fields <- unique(c("barcode", "exact_subclonotype_id",
                                 input$color_by, input$label_by, # UIからの選択
                                 "sample", "seurat_clusters", "BCR_IGH_c_gene", "BCR_pair_CTaa")) # その他含めたい列
      text_fields_to_use <- intersect(potential_meta_fields, names(clone_subset))
      text_fields_to_use <- setdiff(text_fields_to_use, "None") # "None" は除外

      showNotification("Formatting data with formatClones...", id="format_msg", type = "message", duration = 4)
      formatted_data <- tryCatch({
        dowser::formatClones(
          clone_subset,                  # 入力データフレーム
          clone = "raw_clonotype_id",        # クローンID列名
          seq = "BCR_IGH_full_length_nt",
          id = 'barcode'
        )
      }, error = function(e) {
        removeNotification("format_msg")
        showNotification(paste("Error during formatClones:", e$message), type = "error", duration = 10)
        NULL # エラー時はNULLを返す
      })
      removeNotification("format_msg") # 成功時も消す

      # フォーマット結果が有効か確認
      req(formatted_data, inherits(formatted_data, "alakazamFormatClones"), nrow(formatted_data@data) > 0)
      print(formatted_data)
      return(formatted_data)
    })

    # # 4. 系統樹オブジェクトの生成
    # phylogenetic_tree <- reactive({
    #   # フォーマット済みデータが必要
    #   req(phylotree_input_data())
    #   formatted_clones <- phylotree_input_data()

    #   showNotification("Building phylogenetic tree with getTrees (using NJ)...", id="tree_build_msg", type = "message", duration = NULL) # 終わるまで表示
    #   tree_list <- tryCatch({
    #     # build="nj" (近隣結合法) は高速。エラーが出る場合や精度を上げたい場合は "pml" を試す（時間がかかる）
    #     dowser::getTrees(formatted_clones, build = "nj")
    #   }, warning = function(w) {
    #     showNotification(paste("Warning during getTrees:", w$message), type = "warning", duration = 8)
    #     # 警告が出ても処理を試みる
    #     suppressWarnings(trees <- dowser::getTrees(formatted_clones, build = "nj"))
    #     return(trees)
    #   }, error = function(e) {
    #     showNotification(paste("Error during getTrees:", e$message), type = "error", duration = 10)
    #     NULL # エラー時はNULLを返す
    #   })
    #   removeNotification("tree_build_msg") # 処理完了またはエラー時に通知を消す

    #   # 結果が有効か確認
    #   req(tree_list, is.list(tree_list), length(tree_list) > 0)
    #   showNotification("Tree generated successfully.", type = "message", duration = 4)
    #   return(tree_list[[1]]) # リストの最初のtreedataオブジェクトを返す
    # })

    # # 5. 系統樹プロットオブジェクトの生成
    # tree_plot_object <- reactive({
    #   # 系統樹オブジェクトが必要
    #   req(phylogenetic_tree())
    #   tree_data <- phylogenetic_tree() # treedataオブジェクト

    #   # UIからのプロット設定値を取得
    #   color_col <- input$color_by
    #   label_col <- input$label_by
    #   label_size_val <- input$label_size
    #   label_offset_val <- input$label_offset

    #   # Tip (葉ノード) のデータを取得し、ラベル列を確認
    #   tip_data <- as_tibble(tree_data) %>% dplyr::filter(isTip)
    #   if (!(label_col %in% names(tip_data))) {
    #     showNotification(paste("Label column '", label_col, "' not found. Using 'barcode'."), type = "warning")
    #     label_col <- "barcode"
    #     if (!("barcode" %in% names(tip_data))) { label_col <- "label" } # ggtreeのデフォルト
    #   }

    #   # --- ggtree プロット構築 ---
    #   p <- ggtree(tree_data) +
    #        geom_tree(linewidth = 0.6) + # 木の枝
    #        # スケールバーをプロット下部あたりに配置
    #        geom_treescale(x = 0.01, y = -Ntip(tree_data)*0.05, fontsize = 3, linesize = 0.5, offset = Ntip(tree_data)*0.02)

    #   # Tip ポイントの色分け
    #   if (color_col != "None" && color_col %in% names(tip_data)) {
    #      # 色分け対象列が数値やファクターでない場合、ファクターに変換 (凡例のため)
    #      if (!is.factor(tip_data[[color_col]]) && !is.numeric(tip_data[[color_col]])) {
    #          # tree_data@data を直接変更するのは推奨されない場合があるが、ここでは簡便のため行う
    #          tree_data@data <- tree_data@data %>% mutate(!!sym(color_col) := as.factor(!!sym(color_col)))
    #      }
    #      p <- p + geom_tippoint(aes(color = .data[[color_col]]), size = 2) # 色分け実行
    #   } else {
    #      p <- p + geom_tippoint(size = 2, color = "grey50") # 色分けなし
    #      if (color_col != "None") { # None 以外が選択されたのに列がなかった場合
    #          showNotification(paste("Color column '", color_col, "' not found. Tips are not colored."), type = "warning")
    #      }
    #   }

    #   # Tip ラベルの追加
    #   p <- p + geom_tiplab(
    #       aes(label = .data[[label_col]]), # ラベル列を指定
    #       align = TRUE,           # ラベルを右端に揃える
    #       linesize = 0.2,         # ラベルへの線の太さ
    #       size = label_size_val,  # ラベル文字サイズ (UIから)
    #       offset = label_offset_val # 木からの距離 (UIから)
    #   )

    #   # プロットテーマとマージン調整
    #   # 右側のマージンを推定して確保
    #   est_max_label_width_pt <- 0
    #   if(label_col %in% names(tip_data) && nrow(tip_data) > 0 && label_size_val > 0) {
    #        max_chars <- max(nchar(as.character(tip_data[[label_col]])), na.rm = TRUE)
    #        est_max_label_width_pt <- max_chars * label_size_val * 2.5 + label_offset_val * 1000 # 経験的係数
    #   }
    #   right_margin_pt <- max(50, est_max_label_width_pt) # 最低50ptのマージン

    #   p <- p +
    #     theme(
    #       legend.position = "right", # 凡例を右に配置
    #       # 上右下左のマージン (右を広めにとる)
    #       plot.margin = margin(t = 10, r = right_margin_pt, b = 20, l = 10, unit = "pt")
    #     ) +
    #     labs(color = gsub("_", " ", color_col)) # 凡例タイトルを設定 (例: "seurat clusters")

    #   # x軸の表示範囲を調整してラベル用スペースを確保
    #   p <- p + xlim_tree(max(p$data$x, na.rm=TRUE) * 1.1) # x軸を10%広げる

    #   return(p)
    # })

    # # --- Outputs ---

    # # 6. プロット出力
    # output$plot <- renderPlot({
    #   # プロットオブジェクトが生成されるのを待つ
    #   req(tree_plot_object())
    #   plot_obj <- tree_plot_object()

    #   if (!inherits(plot_obj, "ggplot")) {
    #     plot.new(); text(0.5, 0.5, "Failed to generate plot object.", cex = 1.2); return()
    #   }
    #   # プロットオブジェクトを表示
    #   print(plot_obj)
    # },
    # # 高さを動的に調整
    # height = function() {
    #   tree_obj <- phylogenetic_tree() # 依存関係
    #   if (!is.null(tree_obj) && inherits(tree_obj, "phylo")) {
    #     max(400, Ntip(tree_obj) * 18 + 50) # Tip数に応じて高さを調整 + 下部マージン分
    #   } else {
    #     800 # デフォルトの高さ
    #   }
#    })

    # テーブル表示
    output$table <- renderTable({
      # selected_clone_data() がNULLでないことを確認
      req(selected_clone_data())
      selected_clone_data()
    })
  })
}



# # 必要なパッケージをロード（インストールされていない場合は install.packages("package_name") でインストールしてください）
# library(shiny)
# library(ggtree)
# library(dowser)
# # library(alakazam) # 今回は直接は使用しませんが、関連パッケージとして
# library(dplyr)    # データ操作のため
# library(readr)    # ファイル読み込みのため
# library(ggplot2)  # プロット調整のため

# # --- Shiny UI 定義 ---
# phylogeneticTreeUI <- function(id) {
#   ns <- NS(id)
#   sidebarLayout(
#     sidebarPanel(
#       selectInput(ns("selected_clone"), "Select Clone",
#         choices = NULL, # Server側で設定するため NULL に
#         selected = NULL
#       ),
#       actionButton(ns("run_phylotree"), "Generate Phylogenetic Tree", icon = icon("tree"), class = "btn-primary"),

#     ),
#     mainPanel(
#       plotOutput("plot", width = "100%", height = "auto"),
#       tableOutput(ns('table'))
#        # 高さは自動調整 or sliderに連
#     )
#   )
# }

# phylogeneticTreeServer <- function(id, myReactives) {
#   moduleServer(id, function(input, output, session) {

#     observe({
#       req(myReactives$bcr_df)
#       df <- myReactives$bcr_df
#       clones_summary <- df %>%
#         # 必要な列を選択し、NA/空文字を除外
#         dplyr::select(raw_clonotype_id, exact_subclonotype_id) %>%
#         dplyr::filter(!is.na(raw_clonotype_id) & raw_clonotype_id != "",
#                       !is.na(exact_subclonotype_id) & exact_subclonotype_id != "") %>%
#         # raw_clonotype_id ごとにユニークな exact_subclonotype_id の数を数える
#         dplyr::group_by(raw_clonotype_id) %>%
#         dplyr::summarise(n_distinct_subclones = dplyr::n_distinct(exact_subclonotype_id), .groups = 'drop') %>%
#         # ユニークなサブクローンが2つ以上あるものをフィルタリング
#         dplyr::filter(n_distinct_subclones > 1) %>%
#         # raw_clonotype_id のリストを取得
#         dplyr::pull(raw_clonotype_id) %>%
#         # 必要であればソート
#         sort()
#       updateSelectInput(session, "selected_clone",
#         choices = clones_summary)
#     })

#     # 2. 選択されたクローンのデータを抽出 (ボタンクリックでトリガー)
#     selected_clone_data <- eventReactive(input$run_phylotree, {
#       req(myReactives$bcr_df, input$selected_clone)
#       selected_clonotype_id <- input$selected_clone
#       df <- myReactives$bcr_df

#       # データをフィルタリング
#       clone_subset <- df %>%
#         dplyr::filter(
#             raw_clonotype_id == selected_clonotype_id,
#             !is.na(BCR_IGH_full_length_nt) & BCR_IGH_full_length_nt != "",
#             !is.na(barcode) & barcode != ""
#         ) %>%
#         # barcode ごとに最初の行を選択（重複 barcode がある場合への対策）
#         dplyr::distinct(barcode, .keep_all = TRUE)

#       # シーケンス数が2未満の場合は系統樹を構築できない
#       if (nrow(clone_subset) < 2) {
#         showNotification(paste("Not enough valid sequences (< 2) found for clone '", selected_clonotype_id, "' to build a tree."), type = "warning", duration = 8)
#         return(NULL)
#       }
      
#       # データに含まれるメタデータ列を特定（色付け選択肢の更新などに使える）
#       # available_meta_cols <- intersect(names(clone_subset), c("sample", "seurat_clusters", ...))
#       # updateSelectInput(session, "color_by", choices = c("None", available_meta_cols))

#       return(clone_subset)
#     })

#     output$table <- renderTable({
#       selected_clone_data()
#     })



#   })
# }

#     # 2. 選択されたクローンのデータを抽出 (ボタンクリックでトリガー)
#     selected_clone_data <- eventReactive(input$run_phylotree, {
#       req(myReactives$bcr_df, input$selected_clone, input$selected_clone != "")
#       selected_clonotype_id <- input$selected_clone
#       df <- myReactives$bcr_df

#       # 系統樹構築に必要な列が存在するか確認
#       required_cols <- c("barcode", "raw_clonotype_id", "BCR_IGH_full_length_nt")
#       if (!all(required_cols %in% names(df))) {
#         showNotification("Required columns (barcode, raw_clonotype_id, BCR_IGH_full_length_nt) not found.", type = "error", duration = 10)
#         return(NULL)
#       }
      
#       # データをフィルタリング
#       clone_subset <- df %>%
#         dplyr::filter(
#             raw_clonotype_id == selected_clonotype_id,
#             !is.na(BCR_IGH_full_length_nt) & BCR_IGH_full_length_nt != "",
#             !is.na(barcode) & barcode != ""
#         ) %>%
#         # barcode ごとに最初の行を選択（重複 barcode がある場合への対策）
#         dplyr::distinct(barcode, .keep_all = TRUE)

#       # シーケンス数が2未満の場合は系統樹を構築できない
#       if (nrow(clone_subset) < 2) {
#         showNotification(paste("Not enough valid sequences (< 2) found for clone '", selected_clonotype_id, "' to build a tree."), type = "warning", duration = 8)
#         return(NULL)
#       }
      
#       # データに含まれるメタデータ列を特定（色付け選択肢の更新などに使える）
#       # available_meta_cols <- intersect(names(clone_subset), c("sample", "seurat_clusters", ...))
#       # updateSelectInput(session, "color_by", choices = c("None", available_meta_cols))

#       return(clone_subset)
#     })

#     # 3. dowser 用のデータフォーマット (選択データが変更されたら自動更新)
#     phylotree_input_data <- reactive({
#       req(selected_clone_data())
#       clone_subset <- selected_clone_data()
#       if(is.null(clone_subset)) return(NULL)

#       # メタデータとして保持する列を決定
#       # (UIで選択された列 + barcode + その他プロットで有用そうな列)
#       potential_meta_fields <- unique(c("barcode", "exact_subclonotype_id",
#                                  input$color_by, input$label_by,
#                                  "sample", "seurat_clusters", "BCR_IGH_c_gene", "BCR_pair_CTaa"))
#       # データフレームに存在する列のみを選択
#       text_fields_to_use <- intersect(potential_meta_fields, names(clone_subset))
#       # "None" は除外
#       text_fields_to_use <- setdiff(text_fields_to_use, "None")

#       showNotification("Formatting data for tree building...", id="format_msg", type = "message", duration = 4)

#       formatted_data <- tryCatch({
#         dowser::formatClones(
#           clone_subset,
#           clone = "raw_clonotype_id",       # クローンID列名
#           seq = "BCR_IGH_full_length_nt", # IGHのヌクレオチド配列列名
#           seq_id = "barcode",             # ユニークID列名
#           text_fields = text_fields_to_use, # 含めるメタデータ列名
#           num_fields = NULL               # 数値メタデータが必要な場合は指定
#         )
#       }, error = function(e) {
#         removeNotification("format_msg")
#         showNotification(paste("Error during formatClones:", e$message), type = "error", duration = 10)
#         return(NULL)
#       })

#       if (is.null(formatted_data) || nrow(formatted_data@data) == 0) {
#          removeNotification("format_msg")
#          showNotification("formatClones did not return valid data. Check input sequences and IDs.", type = "error", duration = 10)
#          return(NULL)
#       }
#       removeNotification("format_msg")
#       return(formatted_data)
#     })

#     # 4. 系統樹オブジェクトの生成 (フォーマット済みデータが変更されたら自動更新)
#     phylogenetic_tree <- reactive({
#       req(phylotree_input_data())
#       formatted_clones <- phylotree_input_data()
#       if (is.null(formatted_clones)) return(NULL)

#       showNotification("Building phylogenetic tree... (This may take time)", id="tree_build_msg", type = "message", duration = NULL) # 手動で閉じるまで表示

#       tree_list <- tryCatch({
#         # build="nj" (近隣結合法) は高速。build="pml" (最尤法) は時間がかかるが精度が高いとされる。
#         dowser::getTrees(formatted_clones, build = "nj")
#       }, warning = function(w) {
#         showNotification(paste("Warning during getTrees:", w$message), type = "warning", duration = 8)
#         # Warningが出ても処理を試みる
#         suppressWarnings(trees <- dowser::getTrees(formatted_clones, build = "nj"))
#         return(trees)
#       }, error = function(e) {
#         showNotification(paste("Error during getTrees:", e$message), type = "error", duration = 10)
#         return(NULL) # エラー時は NULL を返す
#       })

#       removeNotification("tree_build_msg") # 処理完了またはエラー時に通知を消す

#       if (is.null(tree_list)) {
#         showNotification("Tree generation failed (getTrees returned NULL).", type = "error")
#         return(NULL)
#       } else if (length(tree_list) == 0) {
#         showNotification("Could not generate a phylogenetic tree (getTrees returned empty list).", type = "warning")
#         return(NULL)
#       } else {
#         showNotification("Phylogenetic tree generated successfully.", type = "message", duration = 5)
#         # dowser はリストで返すことが多いが、通常は最初の木を使用
#         return(tree_list[[1]]) # treedata オブジェクトを返す
#       }
#     })

#     # 5. 系統樹プロットオブジェクトの生成 (系統樹オブジェクトが更新されたら自動更新)
#     tree_plot_object <- reactive({
#       req(phylogenetic_tree())
#       tree_data <- phylogenetic_tree() # treedataオブジェクト
#       if (is.null(tree_data)) return(NULL)

#       # プロット設定値を取得
#       color_col <- input$color_by
#       label_col <- input$label_by
#       label_size_val <- input$label_size
#       label_offset_val <- input$label_offset

#       # Tip (葉ノード) のデータを取得
#       tip_data <- as_tibble(tree_data) %>% dplyr::filter(isTip)

#       # ラベル列が存在するか確認
#       if (!(label_col %in% names(tip_data))) {
#         showNotification(paste("Label column '", label_col, "' not found in tree data. Using 'barcode' instead."), type = "warning")
#         label_col <- "barcode" # デフォルトに戻す
#         if (!("barcode" %in% names(tip_data))) {
#              label_col <- "label" # それもなければ ggtree デフォルトの tip ID
#         }
#       }
      
#       # ggtree を使ってプロットを構築
#       p <- ggtree(tree_data) +
#         geom_tree(linewidth = 0.6) + # 木の枝
#         # geom_treescale をプロットの下部中央あたりに配置 (y座標は負の値で指定)
#         geom_treescale(x = 0.01, y = -Ntip(tree_data)*0.05, fontsize = 3, linesize = 0.5, offset = Ntip(tree_data)*0.02)

#       # Tip ポイントの色分け (UIで "None" 以外が選択され、列が存在する場合)
#       if (color_col != "None" && color_col %in% names(tip_data)) {
#          # 色分け対象列のデータ型を確認し、必要ならFactorに変換 (凡例のため)
#          if (!is.factor(tip_data[[color_col]]) && !is.numeric(tip_data[[color_col]])) {
#              tree_data@data <- tree_data@data %>% mutate(!!sym(color_col) := as.factor(!!sym(color_col)))
#          }
#          p <- p + geom_tippoint(aes(color = .data[[color_col]]), size = 2) # .data[[]] で動的に列指定
#       } else {
#          # 色分けしない場合は灰色で点をプロット
#          p <- p + geom_tippoint(size = 2, color = "grey50")
#          if (color_col != "None") { # None 以外が選択されたのに列がなかった場合
#              showNotification(paste("Color column '", color_col, "' not found. Tips are not colored."), type = "warning")
#          }
#       }

#       # Tip ラベルの追加
#       p <- p + geom_tiplab(
#           aes(label = .data[[label_col]]), # ラベル列を動的に指定
#           align = TRUE,           # ラベルを右端に揃える
#           linesize = 0.2,         # ラベルへの線の太さ
#           size = label_size_val,  # ラベル文字サイズ (UIから)
#           offset = label_offset_val # 木からの距離 (UIから)
#       )

#       # プロットテーマとマージン調整
#       # 右側のマージンを確保してラベルがはみ出ないようにする
#       # 最大ラベル長の推定 (簡易的)
#       est_max_label_width_pt <- 0
#       if(label_col %in% names(tip_data) && nrow(tip_data) > 0 && label_size_val > 0) {
#            max_chars <- max(nchar(as.character(tip_data[[label_col]])), na.rm = TRUE)
#            # ポイント単位での大まかな幅を推定 (係数はフォント等に依存するため要調整)
#            est_max_label_width_pt <- max_chars * label_size_val * 2.5 + label_offset_val * 1000
#       }
#       right_margin_pt <- max(50, est_max_label_width_pt) # 最低50ptのマージン

#       p <- p +
#         theme(
#           legend.position = "right", # 凡例を右に配置
#           plot.margin = margin(t = 10, r = right_margin_pt, b = 20, l = 10, unit = "pt") # 上右下左のマージン (右を広めに)
#         ) +
#         # scale_color_viridis_d() # 必要なら離散値用のカラースケールを適用 (option="plasma" など)
#         labs(color = color_col) # 凡例タイトルを設定

#       # x軸の表示範囲を調整してラベル用スペースを確保 (ggtree > 3.0)
#       p <- p + xlim_tree(max(p$data$x, na.rm=TRUE) * 1.1) # x軸を10%広げる (右側)


#       return(p)
#     })

#     # --- Outputs ---

#     # 6. プロット情報 (タイトル代わり)
#     output$plot_info <- renderText({
#       # ボタンが押され、データが準備されたら情報を表示
#       req(selected_clone_data())
#       data_info <- selected_clone_data()
#       # 系統樹オブジェクトが生成されていればTip数も表示
#       tree_obj <- phylogenetic_tree() # NULLかもしれない

#       n_seqs <- if (!is.null(data_info)) nrow(data_info) else 0
#       n_tips <- if (!is.null(tree_obj) && inherits(tree_obj, "phylo")) Ntip(tree_obj) else 0 # phyloクラスか確認

#       paste("Phylogenetic Tree for Clone:", input$selected_clone,
#             "| Input sequences:", n_seqs,
#             "| Tree tips:", n_tips)
#     })

#     # 7. 系統樹プロットのレンダリング
#     output$plot <- renderPlot({
#       # プロットオブジェクトが生成されるのを待つ
#       plot_obj <- tree_plot_object()
#       req(plot_obj) # NULLでないことを確認

#       if (!inherits(plot_obj, "ggplot")) {
#         # ggplotオブジェクトでない場合のエラー表示
#         plot.new(); text(0.5, 0.5, "Failed to generate ggplot object.", cex = 1.2)
#         return()
#       }
#       # プロットオブジェクトを表示
#       print(plot_obj)
#     },
#     # 高さを動的に調整 (Tip数に基づいて)
#     height = function() {
#       tree_obj <- phylogenetic_tree() # 依存関係を設定
#       if (!is.null(tree_obj) && inherits(tree_obj, "phylo")) {
#         # Tip数に応じて高さを調整 (最小400px、Tipあたり約18pxを目安)
#         max(400, Ntip(tree_obj) * 18)
#       } else {
#         600 # デフォルトの高さ
#       }
#     })

#     # 8. プロットダウンロード機能
#     output$download_plot <- downloadHandler(
#       filename = function() {
#         # ファイル名に選択したクローンIDを含める
#         paste0("phylogenetic_tree_", input$selected_clone, "_", Sys.Date(), ".pdf")
#       },
#       content = function(file) {
#         plot_obj <- tree_plot_object() # ダウンロード時にもプロットオブジェクトを再取得
#         tree_obj <- phylogenetic_tree() # Tip数取得のため

#         if (is.null(plot_obj) || !inherits(plot_obj, "ggplot")) {
#           showNotification("Cannot download: Plot object is not available.", type = "error")
#           # 空のPDFを作成してエラーを示す
#           pdf(file, width=7, height=5)
#           plot.new(); text(0.5, 0.5, "Plot could not be generated.")
#           dev.off()
#           return()
#         }

#         # PDFのサイズを動的に決定 (renderPlotの高さ計算ロジックを参考にインチへ変換)
#         plot_height_in <- if (!is.null(tree_obj) && inherits(tree_obj, "phylo")) {
#                              max(6, Ntip(tree_obj) * 0.25) # インチ単位 (Tipあたり0.25 inch)
#                            } else {
#                              8 # デフォルト高さ (インチ)
#                            }
#         plot_width_in <- 8 # 固定幅 (インチ) or 必要なら調整

#         tryCatch({
#           # ggsave で PDF として保存
#           ggsave(file, plot = plot_obj,
#                  width = plot_width_in,
#                  height = plot_height_in,
#                  units = "in",
#                  device = "pdf",
#                  dpi = 300,
#                  limitsize = FALSE) # 大きなプロットの保存を許可
#           showNotification("Plot downloaded successfully as PDF.", type="message", duration=5)
#         }, error = function(e) {
#           showNotification(paste("Error saving PDF:", e$message), type="error", duration=10)
#           # エラー時も空PDFを作成
#           pdf(file, width=7, height=5)
#           plot.new(); text(0.5, 0.5, "Error saving plot to PDF.")
#           dev.off()
#         })
#       }
#     )

#   }) # moduleServer終了
# }

# --- アプリケーションの実行例 ---
# 実際にアプリとして実行するには、以下の部分のコメントを解除し、
# `myReactives$bcr_df` にデータをロードする処理を実装してください。

# ui <- fluidPage(
#   titlePanel("BCR Clonal Phylogenetic Tree Viewer"),
#   phylogeneticTreeUI("phyloTreeModule") # UIモジュールを呼び出し
# )
# 
# server <- function(input, output, session) {
# 
#   # リアクティブな値を保持するリスト
#   myReactives <- reactiveValues(bcr_df = NULL)
# 
#   # --- データ読み込み処理 ---
#   # ここで myReactives$bcr_df にデータをロードします。
#   # 例: ファイルアップロード機能を使う、あるいは固定パスから読み込むなど。
#   # ユーザー提供のデータ例 (bcr_data_example.csv) を読み込む場合:
#   observe({
#     data_file <- "bcr_data_example.csv" # ユーザー提供データを保存したファイル名
#     if (file.exists(data_file)) {
#       bcr_data <- tryCatch({
#         read_csv(data_file, show_col_types = FALSE)
#       }, error = function(e) {
#         showNotification(paste("Error reading data file:", e$message), type="error", duration=15)
#         NULL
#       })
# 
#       if (!is.null(bcr_data)) {
#         # モジュールが必要とする基本列 + UIで使う列の存在チェック
#         required_cols_for_app <- c("barcode", "raw_clonotype_id", "exact_subclonotype_id",
#                                    "BCR_IGH_full_length_nt", "sample", "seurat_clusters", "BCR_IGH_c_gene")
#         missing_cols <- setdiff(required_cols_for_app, names(bcr_data))
# 
#         if (length(missing_cols) > 0) {
#            showNotification(paste("Warning: Missing columns in loaded data, some features might not work:", paste(missing_cols, collapse=", ")), type="warning", duration=15)
#            # 基本列がなければエラーにする
#            base_required <- c("barcode", "raw_clonotype_id", "exact_subclonotype_id", "BCR_IGH_full_length_nt")
#            if(!all(base_required %in% names(bcr_data))){
#                showNotification("Error: Crucial columns missing. Module cannot run.", type="error", duration=20)
#                myReactives$bcr_df <- NULL # データ設定せず終了
#                return()
#            }
#         }
# 
#         # データ型を適切に設定 (例: クラスタリング結果をカテゴリカルに)
#         if ("seurat_clusters" %in% names(bcr_data)) {
#           bcr_data$seurat_clusters <- as.factor(bcr_data$seurat_clusters)
#         }
#         if ("sample" %in% names(bcr_data)) {
#           bcr_data$sample <- as.factor(bcr_data$sample)
#         }
# 
#         myReactives$bcr_df <- bcr_data
#         showNotification("BCR data loaded successfully.", type="message", duration=5)
# 
#       }
#     } else {
#       showNotification(paste("Data file not found:", data_file, ". Please place it in the app directory."), type="error", duration=15)
#       # アプリがクラッシュしないよう、最小限のデータ構造を持つ空のtibbleを設定
#       myReactives$bcr_df <- tibble(
#          barcode=character(), raw_clonotype_id=character(), exact_subclonotype_id=character(),
#          BCR_IGH_full_length_nt=character(), sample=character(), seurat_clusters=character(), BCR_IGH_c_gene=character()
#       )
#     }
#   })
# 
#   # --- モジュールサーバーの呼び出し ---
#   phylogeneticTreeServer("phyloTreeModule", myReactives)
# 
# }
# 
# shinyApp(ui = ui, server = server)