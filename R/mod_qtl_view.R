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
      "Under development.."
    )
  )
}

#' qtl_view Server Functions
#'
#' @import JBrowseR
#' @importFrom shinyjs inlineCSS
#'
#' @noRd 
mod_qtl_view_server <- function(input, output, session, loadMap, loadJBrowse, loadQTL){
  ns <- session$ns
  
}

## To be copied in the UI
# mod_qtl_view_ui("qtl_view_ui_1")

## To be copied in the server
# mod_qtl_view_server("qtl_view_ui_1")
