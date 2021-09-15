#' download UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_download_ui <- function(id){
  ns <- NS(id)
  tagList(
    column(width = 6,
           fluidPage( "Add here options to download tables"
            )
    )
  )
}

#' download Server Functions
#'
#' @noRd 
mod_download_server <- function(input, output, session){
  ns <- session$ns
  
  # Format examples
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$downloadType, ".rds")
    },
    content = function(file) {
      if(input$downloadType == "dosages") {
        filetemp <- readRDS(system.file("ext/dosage.rds", package = "viewpoly"))
      } else if(input$downloadType == "phases") {
        filetemp <- readRDS(system.file("ext/phases.rds", package = "viewpoly"))
      } else if(input$downloadType == "genetic_map") {
        filetemp <- readRDS(system.file("ext/map.rds", package = "viewpoly"))
      } 
      saveRDS(filetemp, file = file)
    }
  )

}

## To be copied in the UI
# mod_download_ui("download_ui_1")

## To be copied in the server
# mod_download_server("download_ui_1")
