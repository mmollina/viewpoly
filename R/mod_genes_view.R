#' genes_view UI Function
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
mod_genes_view_ui <- function(id){
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
                 selectInput(inputId = ns("group"), label = p("Linkage group"), choices = 1:15, selected = 1)
          )
        )
      )
    )
  )
}

#' genes_view Server Functions
#'
#' @import JBrowseR
#' @importFrom shinyjs inlineCSS
#'
#' @noRd 
mod_genes_view_server <- function(input, output, session, loadMap, loadJBrowse, loadQTL){
  ns <- session$ns
  
}

## To be copied in the UI
# mod_genes_view_ui("genes_view_ui_1")

## To be copied in the server
# mod_genes_view_server("genes_view_ui_1")
