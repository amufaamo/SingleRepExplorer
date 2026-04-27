# Gene expression analysis module

# UI
geneExpressionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      selectizeInput(ns("gene"),
        label = "Enter or select gene names:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = 'e.g., CD3E, CD19',
          create = TRUE,
          plugins = list('remove_button')
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] != 'feature_plot'", ns("plot_type")),
        groupByInput(ns),
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'feature_plot'", ns("plot_type")),
        reductionInput(ns),
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'feature_plot' || input['%s'] == 'violin_plot'", ns("plot_type"), ns("plot_type")),
        pointSizeInput(ns),
      ),
      # --- Feature Plot specific ---
      conditionalPanel(
        condition = sprintf("input['%s'] == 'feature_plot'", ns("plot_type")),
        selectInput(ns("feature_col_scale"), "Color Scale",
                    choices = c("Viridis" = "viridis", "Reds" = "Reds",
                                "Blues" = "Blues", "Red-Blue (diverging)" = "RdBu"),
                    selected = "viridis")
      ),
      hr(),
      commonPlotOptions(ns),
    ),
    mainPanel(
      tabsetPanel(id = ns("plot_type"), type = "tabs",
        tabPanel("Feature plot", value = "feature_plot"),
        tabPanel("Violin plot", value = "violin_plot"),
        tabPanel("Dot plot", value = "dot_plot"),
        tabPanel("Heatmap", value = "heatmap")
      ),
      br(),
      downloadButton(ns("download_plot"), "Download Plot (.pptx)"),
      br(), br(),
      plotOutput(ns('plot'))
    )
  )
}

# Server
geneExpressionServer <- function(id, myReactives) {
  moduleServer(id, function(input, output, session) {

    observeEvent(myReactives$seurat_object, {
      req(myReactives$seurat_object)
      myReactives$available_genes <- rownames(myReactives$seurat_object[["RNA"]])

      updateSelectizeInput(
        session,
        "gene",
        choices = myReactives$available_genes,
        selected = c('CD3E', 'CD19'),
        server = TRUE
      )
      update_group_by_select_input(session, myReactives)
      update_reduction_choices(session, myReactives)
    })

    observeEvent(myReactives$grouping_updated, {
      req(myReactives$seurat_object)
      update_group_by_select_input(session, myReactives)
    })

    plot <- reactive({
      req(myReactives$seurat_object, input$gene, myReactives$available_genes)
      shiny::validate(
        shiny::need(length(input$gene) > 0, "Please enter at least one gene name.")
      )

      valid_genes <- intersect(input$gene, myReactives$available_genes)

      shiny::validate(
        shiny::need(length(valid_genes) > 0,
             "No valid genes found. Please check the entered gene names.")
      )

      req(
        if (input$plot_type != "feature_plot") input$group_by else TRUE,
        if (input$plot_type == "feature_plot") input$reduction else TRUE
      )

      so <- myReactives$seurat_object
      if (input$plot_type == "heatmap") {
         tryCatch({
            so <- ScaleData(so, features = valid_genes, verbose = FALSE)
         }, error = function(e) {})
      }

      leg_pos   <- input$legend %||% "right"
      angle     <- as.numeric(input$x_axis_angle %||% "45")
      hjust_val <- if (angle == 0) 0.5 else 1

      # Feature plot color scale
      feature_cols <- switch(input$feature_col_scale %||% "viridis",
        "viridis" = c("#440154", "#31688e", "#35b779", "#fde725"),
        "Reds"    = c("lightgrey", "red"),
        "Blues"   = c("lightgrey", "blue"),
        "RdBu"    = c("blue", "lightgrey", "red"),
        c("lightgrey", "blue")
      )

      dot_scale      <- input$dot_scale %||% 6
      heatmap_lbl_sz <- input$heatmap_label_size %||% 5.5

      tryCatch({
        switch(input$plot_type,
          "feature_plot" = FeaturePlot(myReactives$seurat_object, features = valid_genes,
                            reduction = input$reduction, pt.size = input$point_size,
                            cols = feature_cols) +
                           theme(legend.position = leg_pos),
          "violin_plot"  = VlnPlot(myReactives$seurat_object, features = valid_genes,
                            group.by = input$group_by, pt.size = input$point_size) +
                           theme(legend.position = leg_pos),
          "dot_plot"     = DotPlot(myReactives$seurat_object, features = valid_genes,
                            group.by = input$group_by, scale = dot_scale) +
                           theme(
                             legend.position = leg_pos,
                             axis.text.x     = element_text(angle = angle, hjust = hjust_val)
                           ),
          "heatmap"      = DoHeatmap(so, features = valid_genes, group.by = input$group_by,
                            size = heatmap_lbl_sz)
        )
      }, error = function(e) {
        shiny::validate(
          shiny::need(FALSE, paste("An error occurred while generating the plot:", e$message))
        )
      })
    })

    output$plot <- renderPlot({
      req(plot())
      plot()
    }, width = reactive(input$plot_width %||% 700), height = reactive(input$plot_height %||% 500))

    output$download_plot <- downloadHandler(
      filename = function() {
        message("[GeneExp Download] filename() called")
        plot_type <- isolate(input$plot_type)
        genes <- isolate(input$gene)

        req(myReactives$available_genes)
        valid_genes <- intersect(genes, myReactives$available_genes)

        gene_string <- if (length(valid_genes) > 3) {
          paste0(paste(head(valid_genes, 3), collapse="_"), "_etc")
        } else {
          paste(valid_genes, collapse="_")
        }
        if (gene_string == "") gene_string <- "no_valid_genes"

        fname <- sprintf("%s_%s.pptx", plot_type, gene_string)
        message("[GeneExp Download] filename = ", fname)
        fname
      },
      content = function(file) {
        message("[GeneExp Download] content() called, target file = ", file)
        
        tryCatch({
          message("[GeneExp Download] Step 1: Getting plot object...")
          p <- plot()
          req(p)
          message("[GeneExp Download] Step 1 OK: plot object class = ", paste(class(p), collapse=", "))
          
          plot_width  <- input$plot_width %||% 500
          plot_height <- input$plot_height %||% 500
          message("[GeneExp Download] Step 2: plot_width = ", plot_width, ", plot_height = ", plot_height)
          
          width_in  <- plot_width / 72
          height_in <- plot_height / 72
          
          # Step 3: Save as PNG first (this is the safest approach)
          message("[GeneExp Download] Step 3: Saving as temp PNG...")
          tmp_img <- tempfile(fileext = ".png")
          
          grDevices::png(tmp_img, width = plot_width * 3, height = plot_height * 3, res = 300)
          tryCatch({
            print(p)
          }, finally = {
            grDevices::dev.off()
          })
          
          message("[GeneExp Download] Step 3 OK: PNG saved, size = ", file.info(tmp_img)$size, " bytes")
          
          # Step 4: Create PPTX with embedded image
          message("[GeneExp Download] Step 4: Creating PPTX...")
          pptx <- officer::read_pptx()
          pptx <- officer::add_slide(pptx, layout = "Blank")
          pptx <- officer::ph_with(
            pptx,
            value = officer::external_img(src = tmp_img, width = width_in, height = height_in),
            location = officer::ph_location(
              left = 0.5, top = 0.5,
              width = width_in, height = height_in
            )
          )
          message("[GeneExp Download] Step 4 OK: PPTX object created")
          
          # Step 5: Write PPTX to file
          message("[GeneExp Download] Step 5: Writing PPTX to target file...")
          print(pptx, target = file)
          message("[GeneExp Download] Step 5 OK: PPTX written, size = ", file.info(file)$size, " bytes")
          
          # Cleanup
          unlink(tmp_img)
          message("[GeneExp Download] DONE - download should succeed")
          
        }, error = function(e) {
          message("[GeneExp Download] ERROR: ", conditionMessage(e))
          message("[GeneExp Download] TRACEBACK: ")
          message(paste(capture.output(traceback()), collapse = "\n"))
          
          # Emergency fallback: write a minimal valid pptx
          tryCatch({
            message("[GeneExp Download] Attempting emergency fallback...")
            pptx <- officer::read_pptx()
            pptx <- officer::add_slide(pptx, layout = "Blank")
            pptx <- officer::ph_with(
              pptx,
              value = paste("Error generating plot:", conditionMessage(e)),
              location = officer::ph_location(left = 1, top = 3, width = 8, height = 1)
            )
            print(pptx, target = file)
            message("[GeneExp Download] Emergency fallback PPTX written")
          }, error = function(e2) {
            message("[GeneExp Download] EMERGENCY FALLBACK ALSO FAILED: ", conditionMessage(e2))
          })
        })
      }
    )

    output$available_feature <- downloadHandler(
      filename = function() {"available_features.txt"},
      content = function(file) {
        req(myReactives$available_genes)
        write.table(sort(myReactives$available_genes), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
    )
  })
}
