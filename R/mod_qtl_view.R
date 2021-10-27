#' qtl_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
#' @importFrom RColorBrewer brewer.pal
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
          
          column(4,
                 checkboxGroupInput(ns("phenotypes"),
                                    label = h4("Phenotypes"),
                                    choices = "This will be updated",
                                    selected = "This will be updated"), br(),
                 
          ),
          column(5,
                 checkboxGroupInput(ns("group"),
                                    label = h4("Linkage group"),
                                    choices = "This will be updated",
                                    selected = "This will be updated"), br(),
          ),
          column(12,
                 plotOutput(ns("plot_qtl"), 
                            click=ns("plot_click")),
                 tableOutput(ns("info")),
                 plotOutput(ns("effects"),  width = "100%")
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
    
    updateCheckboxGroupInput(session, "group",
                             label="Linkage group",
                             choices = group_choices,
                             selected= group_choices[[1]])
    
    # Dynamic QTLs
    pheno_choices <- as.list(unique(loadQTL()[[1]]$pheno))
    names(pheno_choices) <- unique(loadQTL()[[1]]$pheno)
    
    updateCheckboxGroupInput(session, "phenotypes",
                             label = "Phenotypes",
                             choices = pheno_choices,
                             selected=unlist(pheno_choices)[1])
  })
  
  qtl.data <- reactive({
    idx <- which(names(loadQTL()[[3]]$results) %in% input$phenotypes)
    pl <- plot_profile(lgs.info = loadQTL()[[2]], model = loadQTL()[[3]],
                       pheno.col = idx,
                       lgs.id = input$group, by_range=F, plot = F)
    pl
  })
  
  
  output$plot_qtl <- renderPlot({
    only_plot_profile(qtl.data()[[1]], 
                      qtl.data()[[2]], 
                      qtl.data()[[3]], 
                      qtl.data()[[4]], 
                      qtl.data()[[5]], 
                      qtl.data()[[6]],
                      qtl.data()[[7]],
                      qtl.data()[[8]],
                      qtl.data()[[9]])
  })
  
  output$info <- renderTable({
    req(input$plot_click)
    nearPoints(qtl.data()[[2]], input$plot_click, xvar = "Position (cM)", yvar = "y.dat")
  })
  
  output$effects <- renderPlot({
    req(input$plot_click)
    df <- nearPoints(qtl.data()[[2]], input$plot_click, xvar = "Position (cM)", yvar = "y.dat")
    plots <- plot_qtlpoly.effects(x = loadQTL()[[4]], 
                                  pheno.col = which(unique(loadQTL()[[1]]$pheno) %in% as.character(df$Trait)), 
                                  df.info = loadQTL()[[1]], lgs = df$LG, position = df$`Position (cM)`)
    ggarrange(plotlist = unlist(plots, recursive = F), ncol = 4)
  })
}

## To be copied in the UI
# mod_qtl_view_ui("qtl_view_ui_1")

## To be copied in the server
# mod_qtl_view_server("qtl_view_ui_1")
