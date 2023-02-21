#' hidecan_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS useShinyjs
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS tagList  
#' 
#' @noRd 
#'
mod_hidecan_view_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      verticalLayout(
        fluidRow(
          column(width = 12,
                 div(style = "position:absolute;right:1em;", 
                     div(
                       actionButton(ns("goMap"), "Go to Map", icon("arrow-circle-left", verify_fa = FALSE), class = "btn btn-primary"))
                 )
          ),
          tags$h2(tags$b("HIDECAN")), br(), hr(),
          column(12,
                 column(12,
                        box(
                          background = "light-blue",
                          "Required inputs (*)", br(),
                        )
                 ),
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", title = "Define thresholds *",
                            ## Input sliders for GWAS score threshold
                                sliderInput(ns("score_thr_gwas"), "Score threshold for GWAS results", value = 4, min = 0, max = 10, step = 0.1),br(),
                            ## Input sliders for DE score and log2FC threshold
                                sliderInput(ns("score_thr_de"), "Score threshold for DE results", value = 1.3, min = 0, max = 10, step = 0.1), br(),
                                sliderInput(ns("log2fc_thr_de"), "log2(fold-change) threshold for DE results", value = 1, min = 0, max = 10, step = 0.1)
                         
                            
                        )
                 ),
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", collapsible = TRUE,  collapsed = TRUE, title = "HIDECAN plot options",
                            textInput(ns("title"), "Title"),
                            textInput(ns("subtitle"), "Subtitle"),
                            fluidRow(
                              column(6,
                                     numericInput(ns("nrows"), "Number of rows", value = NULL, min = 1, max = Inf)),
                              column(6,
                                     numericInput(ns("ncols"), "Number of columns", value = 2, min = 1, max = Inf))
                            ),
                            selectInput(ns("legend_position"), "Legend position", c("bottom", "top", "left", "right", "none")),
                            fluidRow(
                              column(6,
                                     numericInput(ns("point_size"), "Point size", value = 3, min = 0, max = Inf))
                            ),
                            fluidRow(
                              column(6,
                                     numericInput(ns("label_size"), "Label size", value = 3.5, min = 0, max = Inf)),
                              column(6,
                                     numericInput(ns("label_padding"), "Label padding", value = 0.15, min = 0, max = Inf))
                            ),
                            checkboxInput(ns("colour_genes_by_score"), "Colour genes by score?", value = TRUE),
                            checkboxInput(ns("remove_empty_chrom"), "Remove empty chromosomes?", value = TRUE)
                        )
                 )
          ),
          column(12,
                 box(id = ns("box_hidecan"), width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = actionLink(inputId = ns("hidecanID"), label = "HIDECAN plot"),
                     column(12,
                            box(
                              background = "light-blue",
                              "* HIDECAN analysis files or viewpoly object or example dataset (check `Input data` tab)"
                            )
                     ), 
                     column(12,
                            column(3,
                                   useShinyjs(),
                                   tags$head(tags$style(".butt{background-color:#add8e6; border-color: #add8e6; color: #337ab7;}")),
                                   downloadButton(ns('bn_download'), "Download", class = "butt")
                            ),
                            column(3,
                                   radioButtons(ns("fformat"), "File type", choices=c("png","tiff","jpeg","pdf", "RData"), selected = "png", inline = T)
                            ),                     
                            column(2,
                                   numericInput(ns("width_hidecan"), "Width (mm)", value = 180),
                            ),
                            column(2,
                                   numericInput(ns("height_hidecan"), "Height (mm)", value = 120),
                            ),
                            column(2,
                                   numericInput(ns("dpi_hidecan"), "DPI", value = 300)
                            )), br(), 
                     column(12,
                            hr(),
                            uiOutput(ns("plot.ui"))
                     )
                 )
          )
        )
      )
    )
  )
}

# Inputs for tests
# input <- list()
# input$colour_genes_by_score <- TRUE
# input$remove_empty_chrom <- TRUE
# input$title <- NULL 
# input$subtitle <- NULL
# input$ncols <- 2
# input$legend_position <- "bottom"
# input$point_size <- 3
# input$label_size <- 3.5
# input$label_padding <- 0.15
# input$score_thr_gwas <- 4
# input$score_thr_de <- 1.3
# input$log2fc_thr_de <- 1

#' hidecan_view Server Functions
#'
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom dplyr `%>%`
#' @importFrom shinyjs js
#' @noRd 
mod_hidecan_view_server <- function(input, output, session, 
                                    loadHidecan,
                                    parent_session){
  ns <- session$ns
  
  observeEvent(input$goMap, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "map")
  })
  
  observe({
    if(!is.null(loadHidecan()$GWASpoly)){
      cat("alo!")
      shinyjs::disable("score_thr_gwas")
      shinyjs::disable("score_thr_de")
      shinyjs::disable("log2fc_thr_de")
    } else if(!is.null(loadHidecan()$GWAS)){
      updateSliderInput(inputId = ns("score_thr_gwas"), max = round(max(loadHidecan()$GWAS$score),2), step =  round(max(loadHidecan()$GWAS$score)/20,1))
    } else {
      shinyjs::disable("score_thr_gwas")
    }
    
    if(!is.null(loadHidecan()$DE)){
      if(!is.null(loadHidecan()$DE$score))
        updateSliderInput(inputId = ns("score_thr_de"), max = round(max(loadHidecan()$DE$score),2), step = round(max(loadHidecan()$DE$score)/20,1))
      updateSliderInput(inputId = ns("log2fc_thr_de"), max = round(max(abs(loadHidecan()$DE$log2FoldChange)),2), step = round(max(abs(loadHidecan()$DE$log2FoldChange))/10,1))
    } else {
      shinyjs::disable("score_thr_de")
      shinyjs::disable("log2fc_thr_de")
    }
    
  })
  
  observe({
    if (!is.null(hidecan_data()) & input$width_hidecan > 1 & input$height_hidecan > 1 & input$dpi_hidecan > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download")
    } else {
      shinyjs::disable("bn_download")
    }
  })
  
  plot_nrows <- reactive({
    res <- input$nrows
    if(missing(res) | is.na(res)) res <- NULL
    res
  })
  
  hidecan_data <- reactive({
    # temp <- load(system.file("ext/gwaspoly_thre.RData", package = "viewpoly"))
    # gwaspoly_res_thr <- get(temp)
    
    # loadHidecan() <- get_example_data()
    # loadHidecan() <- prepare_hidecan_examples()
    plot_list <- plot_list_thr <- list(NULL)
    if(!is.null(loadHidecan()[["GWAS"]])){
      plot_list[[1]] <- GWAS_data(loadHidecan()[["GWAS"]])
      plot_list_thr[[1]] <- apply_threshold(plot_list[[1]], 
                                            score_thr = input$score_thr_gwas)
    }
    
    if(!is.null(loadHidecan()[["DE"]])){
      plot_list[[2]] <- DE_data(loadHidecan()[["DE"]])
      plot_list_thr[[2]] <- apply_threshold(plot_list[[2]], 
                                            score_thr = input$score_thr_de,
                                            log2fc_thr = input$log2fc_thr_de)
    }
    
    if(!is.null(loadHidecan()[["CAN"]])){
      plot_list[[3]] <- CAN_data(loadHidecan()[["CAN"]])
      plot_list_thr[[3]] <- apply_threshold(plot_list[[3]])
    }
    
    idx <- sapply(plot_list, is.null)
    if(!all(idx)){
      if(length(which(idx)) > 0) plot_list <- plot_list[-which(idx)]
      chrom_length <- combine_chrom_length(plot_list)    
    }
    
    idx <- sapply(plot_list_thr, is.null)
    if(!all(idx)){
      if(length(which(idx)) > 0) plot_list_thr <- plot_list_thr[-which(idx)]
    }
    
    if(!is.null(loadHidecan()[["GWASpoly"]])){
      if(!inherits(loadHidecan()[["GWASpoly"]], "GWASpoly.thresh")) stop("It is not a GWASpoly.thresh object.")
      p <- hidecan_plot_from_gwaspoly(
        loadHidecan()[["GWASpoly"]],
        remove_empty_chrom = input$remove_empty_chrom,
        title = input$title,
        subtitle = input$subtitle,
        n_rows = plot_nrows(),
        n_cols = input$ncols,
        legend_position = input$legend_position,
        point_size = input$point_size,
        label_size = input$label_size,
        label_padding = input$label_padding  
      )
    } else {
      p <- create_hidecan_plot(plot_list_thr,
                               chrom_length,
                               colour_genes_by_score = input$colour_genes_by_score,
                               remove_empty_chrom = input$remove_empty_chrom,
                               title = input$title,
                               subtitle = input$subtitle,
                               n_rows = plot_nrows(),
                               n_cols = input$ncols,
                               legend_position = input$legend_position,
                               point_size = input$point_size,
                               label_size = input$label_size,
                               label_padding = input$label_padding)
    }
    p
  })
  
  output$plot_hidecan <- renderPlot({
    validate(
      need(!is.null(loadHidecan()), "Upload HIDECAN information in the upload session to access this feature.")
    )
    hidecan_data()
  })
  
  
  plotHeight <- reactive({
    
    validate(
      need(!all(c(is.null(loadHidecan()$GWAS),is.null(loadHidecan()$GWASpoly))), "Upload HIDECAN information in upload session to access this feature."),
    )
    
    if(!is.null(loadHidecan()$GWAS))
      n.chr <- length(unique(loadHidecan()$GWAS$chromosome))
    else if(!is.null(loadHidecan()$GWASpoly))
      n.chr <- length(unique(loadHidecan()$GWASpoly@map$Chrom))
    
    size <- (n.chr/input$ncols)*110
    
    size
  })
  
  output$plot.ui <- renderUI({
    plotOutput(ns("plot_hidecan"), height =  plotHeight())
  })
  
  # HIDECAN download
  fn_downloadname <- reactive({
    seed <- sample(1:1000,1)
    if(input$fformat=="png") filename <- paste0("hidecan","_",seed,".png")
    if(input$fformat=="tiff") filename <- paste0("hidecan","_",seed,".tiff")
    if(input$fformat=="jpeg") filename <- paste0("hidecan","_",seed,".jpg")
    if(input$fformat=="pdf") filename <- paste0("hidecan","_",seed,".pdf")
    if(input$fformat=="RData") filename <- paste0("hidecan","_",seed,".RData")
    return(filename)
  })
  
  # download profile 
  fn_download <- function()
  {
    if(input$fformat!="RData"){
      ggsave(hidecan_data(), filename = fn_downloadname(), 
             width = input$width_hidecan, height = input$height_hidecan, 
             units = "mm", dpi = input$dpi_hidecan)    
    } else save(hidecan_data(), file = fn_downloadname())
  }
  
  observe({
    if (!is.null(hidecan_data()) & input$width_hidecan > 1 & input$height_hidecan > 1 & input$dpi_hidecan > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download")
    } else {
      shinyjs::disable("bn_download")
    }
  })
  
  # download handler
  output$bn_download <- downloadHandler(
    filename = fn_downloadname,
    content = function(file) {
      fn_download()
      file.copy(fn_downloadname(), file, overwrite=T)
      file.remove(fn_downloadname())
    }
  )
}

## To be copied in the UI
# mod_hidecan_view_ui("hidecan_view_ui_1")

## To be copied in the server
# mod_hidecan_view_server("hidecan_view_ui_1")
