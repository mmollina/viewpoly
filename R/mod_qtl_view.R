#' qtl_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @import shinydashboard
#' @import shinyWidgets
#' 
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_qtl_view_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      verticalLayout(
        fluidRow(
          column(width = 12,
                 div(style = "position:absolute;right:1em;", 
                     actionButton(ns("server_off"), "Exit",icon("times-circle"), class = "btn btn-danger"),
                 )
          ),
          tags$h2(tags$b("View QTL")), br(), hr(),
          column(6,
                 column(6,
                        box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE, status="primary",
                            pickerInput(ns("group"),
                                        label = h6("Select linkage groups"),
                                        choices = "This will be updated",
                                        selected = "This will be updated",
                                        options = list(
                                          `actions-box` = TRUE, 
                                          size = 10,
                                          `selected-text-format` = "count > 3"
                                        ), 
                                        multiple = TRUE)
                        )
                 ),
                 column(6,
                        box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = FALSE, status="primary",
                            pickerInput(ns("phenotypes"),
                                        label = h6("Select phenotypes"),
                                        choices = "This will be updated",
                                        selected = "This will be updated",
                                        options = list(
                                          `actions-box` = TRUE, 
                                          size = 10,
                                          `selected-text-format` = "count > 3"
                                        ), 
                                        multiple = TRUE)
                        )
                 )
          ),
          column(12,
                 box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("LOD curve"),
                     column(2,
                            downloadBttn(ns('bn_download'), style = "gradient")
                     ),
                     column(10,
                            radioButtons(ns("fformat"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
                     ), br(), 
                     column(12,
                            hr(),
                            plotOutput(ns("plot_qtl"), 
                                       click=ns("plot_click"), brush = ns("plot_brush"))
                     ),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Effects"),
                         column(2,
                                downloadBttn(ns('bn_download_effects'), style = "gradient")
                         ),
                         column(10,
                                radioButtons(ns("fformat_effects"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
                         ), br(),
                         column(12,
                                hr(),
                                uiOutput(ns("plot.ui"))
                         )
                     ), br(),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("QTL info"),
                         DT::dataTableOutput(ns("info"))
                     ),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Breeding values"),
                         DT::dataTableOutput(ns("breeding_values"))
                     )
                 )
          )
        )
      )
    )
  )
}

#' qtl_view Server Functions
#'
#' @importFrom ggpubr ggarrange
#' @import shinydashboard
#' 
#' @noRd 
mod_qtl_view_server <- function(input, output, session, loadMap, loadJBrowse, loadQTL){
  ns <- session$ns
  
  observe({
    # Dynamic linkage group number
    group_choices <- as.list(1:length(loadMap()$dp))
    names(group_choices) <- 1:length(loadMap()$dp)
    
    updatePickerInput(session, "group",
                      label="Select linkage groups",
                      choices = group_choices,
                      selected= group_choices[[1]])
    
    
    # Dynamic QTL
    pheno_choices <- as.list(unique(loadQTL()$profile$pheno))
    names(pheno_choices) <- unique(loadQTL()$profile$pheno)
    
    updatePickerInput(session, "phenotypes",
                      label = "Select phenotypes",
                      choices = pheno_choices,
                      selected=unlist(pheno_choices)[1])
    
  })
  
  qtl.data <- reactive({
    idx <- which(unique(loadQTL()$profile$pheno) %in% input$phenotypes)
    pl <- plot_profile(profile = loadQTL()$profile, 
                       qtl_info = loadQTL()$qtl_info, 
                       selected_mks = loadQTL()$selected_mks,
                       pheno.col = idx,
                       lgs.id = as.numeric(input$group), 
                       by_range=F, plot = F)
    
    pl
  })
  
  output$plot_qtl <- renderPlot({
    only_plot_profile(pl.in = qtl.data())
  })
  
  output$effects <- renderPlot({
    if(!is.null(input$plot_brush)){
      df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      df <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "x", yvar = "y.dat")
    } else {
      stop("Select a point or region on LOD profile graphic.") 
    }
    plots <- plot_qtlpoly.effects(loadQTL()$qtl_info, loadQTL()$effects,
                                  pheno.col = as.character(df$Trait), 
                                  lgs = df$LG, position = df$`Position (cM)`)
    
    rows <- ceiling(length(plots)/4)
    if(rows == 0) rows <- 1
    
    ggarrange(plotlist = plots, ncol = 4, nrow = rows)
  })
  
  plotHeight <- reactive({
    if(!is.null(input$plot_brush)){
      dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      dframe <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "x", yvar = "y.dat")
    } else {
      stop("Select a point or region on LOD profile graphic.")
    }
    counts <- nrow(dframe)
    counts <- ceiling(counts/4)
    if(counts == 0) counts <- 1
    size <- counts*350
    size
  })
  
  output$plot.ui <- renderUI({
    plotOutput(ns("effects"), height =  plotHeight())
  })
  
  output$info <- DT::renderDataTable({
    if(!is.null(input$plot_brush)){
      dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      dframe <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "x", yvar = "y.dat")
    } else {
      stop("Select a point or region on graphic.")
    }
    dframe <- dframe[,-c(dim(dframe)[2]-1,dim(dframe)[2])]
    colnames(dframe)[c(2,4,5,6)] <- c("Linkage group", "Inferior interval", "Superior interval", "p-value")
    DT::datatable(dframe, extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtlp',
                    buttons = c('copy', 'csv', 'excel', 'pdf')
                  ),
                  class = "display")
  })
  
  # Breeding values
  output$breeding_values <- DT::renderDataTable({
    if(!is.null(input$plot_brush)){
      dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      dframe <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "x", yvar = "y.dat")
    } else {
      stop("Select a point or region on graphic.")
    }
    
    pos <- split(dframe$`Position (cM)`, dframe$Trait)
    dt <- breeding_values(loadQTL()$qtl_info, loadQTL()$probs, 
                          loadQTL()$selected_mks, loadQTL()$blups, 
                          loadQTL()$beta.hat, pos)
    DT::datatable(dt, extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtlp',
                    buttons = c('copy', 'csv', 'excel', 'pdf')
                  ),
                  class = "display")
  })
  
  # Download profile
  # create filename
  fn_downloadname <- reactive({
    
    seed <- sample(1:1000,1)
    if(input$fformat=="png") filename <- paste0("profile","_",seed,".png")
    if(input$fformat=="tiff") filename <- paste0("profile","_",seed,".tif")
    if(input$fformat=="jpeg") filename <- paste0("profile","_",seed,".jpg")
    if(input$fformat=="pdf") filename <- paste0("profile","_",seed,".pdf")
    return(filename)
  })
  
  # download profile 
  fn_download <- function()
  {
    p <- only_plot_profile(pl.in = qtl.data())
    ggsave(p, file = fn_downloadname())    
  }
  
  # download handler
  output$bn_download <- downloadHandler(
    filename = fn_downloadname,
    content = function(file) {
      fn_download()
      file.copy(fn_downloadname(), file, overwrite=T)
    }
  )
  
  # Download effects
  # create filename
  fn_downloadname_effects <- reactive({
    
    seed <- sample(1:1000,1)
    if(input$fformat_effects=="png") filename <- paste0("profile","_",seed,".png")
    if(input$fformat_effects=="tiff") filename <- paste0("profile","_",seed,".tif")
    if(input$fformat_effects=="jpeg") filename <- paste0("profile","_",seed,".jpg")
    if(input$fformat_effects=="pdf") filename <- paste0("profile","_",seed,".pdf")
    return(filename)
  })
  
  # download 
  fn_download_effects <- function()
  {
    if(!is.null(input$plot_brush)){
      df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      df <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "x", yvar = "y.dat")
    } else {
      stop("Select a point or region on LOD profile graphic.") 
    }
    plots <- plot_qtlpoly.effects(loadQTL()$qtl_info, loadQTL()$effects,
                                  pheno.col = as.character(df$Trait), 
                                  lgs = df$LG, position = df$`Position (cM)`)
    
    rows <- ceiling(length(plots)/4)
    if(rows == 0) rows <- 1
    
    p <- ggarrange(plotlist = plots, ncol = 4, nrow = rows)
    ggsave(p, file = fn_downloadname_effects(), height = plotHeight(), 
           width = 2.5*plotHeight(), units = "mm")    
  }
  
  # download handler
  output$bn_download_effects <- downloadHandler(
    filename = fn_downloadname_effects,
    content = function(file) {
      fn_download_effects()
      file.copy(fn_downloadname_effects(), file, overwrite=T)
    }
  )
}

# position = c(146.02,146.02,147.31,147.31,144.38,150.05,153.82,158.13)
# lgs = 12.00
# pheno.col <- which(unique(qtls[[1]]$pheno) %in% c("Starch", "Sucrose", "Bcar", "Bcar0", "Protein", "DM", "Maltose", "Ca"))

## To be copied in the UI
# mod_qtl_view_ui("qtl_view_ui_1")

## To be copied in the server
# mod_qtl_view_server("qtl_view_ui_1")
