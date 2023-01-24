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
          inlineCSS(".form-group {margin-bottom: 0;}
                                .irs-with-grid {bottom: 0px;}
                                .irs-grid {height: 13px;}
                                .irs-grid-text {height: 0px;}
                                "
          ),
          tags$h2(tags$b("Under development ;)")), br(), hr(),
          
        )
      )
    )
  )
}

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
  
  observe({
    str(loadHidecan())
  })
  
}

## To be copied in the UI
# mod_hidecan_view_ui("hidecan_view_ui_1")

## To be copied in the server
# mod_hidecan_view_server("hidecan_view_ui_1")
