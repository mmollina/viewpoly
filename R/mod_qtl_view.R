#' qtl_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
#' @importFrom RColorBrewer brewer.pal
#' @import ggpubr
#' @import shinydashboard
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
                        box(width = 6, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE, status="primary", title = h4("Linkage group"),
                            checkboxGroupInput(ns("group"),
                                               label = h6("Select linkage groups"),
                                               choices = "This will be updated",
                                               selected = "This will be updated"), 
                            actionLink(ns("selectall1"),"Select all groups"), br(),
                        )
                 ),
                 column(6,
                        box(width = 6, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Phenotypes"),
                            checkboxGroupInput(ns("phenotypes"),
                                               label = h6("Select phenotypes"),
                                               choices = "This will be updated",
                                               selected = "This will be updated"), 
                            actionLink(ns("selectall2"),"Select all phenotypes"), br(),
                        )
                 )
          ),
          column(12,
                 plotOutput(ns("plot_qtl"), 
                            click=ns("plot_click"), brush = ns("plot_brush")),
                 uiOutput(ns("plot.ui")),
                 tableOutput(ns("info"))
          )
        )
      )
    )
  )
}

#' qtl_view Server Functions
#'
#' @import JBrowseR
#' @import ggpubr
#' @importFrom shinyjs inlineCSS
#'
#' @noRd 
mod_qtl_view_server <- function(input, output, session, loadMap, loadJBrowse, loadQTL){
  ns <- session$ns
  
  observe({
    # Dynamic linkage group number
    group_choices <- as.list(1:length(loadMap()$dp))
    names(group_choices) <- 1:length(loadMap()$dp)
    
    if (input$selectall1%%2 == 0)
    {
      updateCheckboxGroupInput(session, "group",
                               label="Select linkage groups",
                               choices = group_choices,
                               selected= group_choices[[1]])
    }
    else
    {
      updateCheckboxGroupInput(session, "group",
                               label="Select linkage groups",
                               choices = group_choices,
                               selected= unlist(group_choices))
    }
    
    # Dynamic QTLs
    pheno_choices <- as.list(unique(loadQTL()[[1]]$pheno))
    names(pheno_choices) <- unique(loadQTL()[[1]]$pheno)
    
    if (input$selectall2%%2 == 0)
    {
      updateCheckboxGroupInput(session, "phenotypes",
                               label = "Select phenotypes",
                               choices = pheno_choices,
                               selected=unlist(pheno_choices)[1])
    }
    else
    {
      updateCheckboxGroupInput(session, "phenotypes",
                               label = "Select phenotypes",
                               choices = pheno_choices,
                               selected=unlist(pheno_choices))
    }
  })
  
  qtl.data <- reactive({
    idx <- which(names(loadQTL()[[3]]$results) %in% input$phenotypes)
    pl <- plot_profile(lgs.info = loadQTL()[[2]], model = loadQTL()[[3]],
                       pheno.col = idx,
                       lgs.id = input$group, by_range=F, plot = F)
    pl
  })
  
  output$plot_qtl <- renderPlot({
    only_plot_profile(pl.in = qtl.data())
  })
  
  output$effects <- renderPlot({
    if(!is.null(input$plot_brush)){
      df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "Position (cM)", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      df <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "Position (cM)", yvar = "y.dat")
    } else {
      stop("Select a point or region on graphic.") 
    }
    
    plots <- plot_qtlpoly.effects(x = loadQTL()[[4]], 
                                  pheno.col = which(unique(loadQTL()[[1]]$pheno) %in% as.character(df$Trait)), 
                                  df.info = loadQTL()[[1]], lgs = df$LG, position = df$`Position (cM)`)
    
    rows <- ceiling(length(plots)/4)
    if(rows == 0) rows <- 1

    ggarrange(plotlist = plots, ncol = 4, nrow = rows)
  })
  
  plotHeight <- reactive({
    if(!is.null(input$plot_brush)){
      dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "Position (cM)", yvar = "y.dat")
    } else if(!is.null(input$plot_click)){
      dframe <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "Position (cM)", yvar = "y.dat")
    } else {
      stop("Select a point or region on graphic.")
    }
    counts <- nrow(dframe)
    counts <- ceiling(counts/4)
    if(counts == 0) counts <- 1
    size <- counts*350
    size
  })

  output$plot.ui <- renderUI({
    plotOutput(ns("effects"), height = plotHeight())
  })
  
  output$info <- renderTable({
    req(input$plot_brush)
    brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "Position (cM)", yvar = "y.dat")
  })
}

# position = c(146.02,146.02,147.31,147.31,144.38,150.05,153.82,158.13)
# lgs = 12.00
# pheno.col <- which(unique(loadQTL()[[1]]$pheno) %in% c("Starch", "Sucrose", "Bcar", "Bcar0", "Protein", "DM", "Maltose", "Ca"))

## To be copied in the UI
# mod_qtl_view_ui("qtl_view_ui_1")

## To be copied in the server
# mod_qtl_view_server("qtl_view_ui_1")
