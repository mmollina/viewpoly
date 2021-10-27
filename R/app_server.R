#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic
  # Upload size
  options(shiny.maxRequestSize=200*1024^2)
  
  ## Start modules
  datas <- callModule(mod_upload_server,"upload_ui_1", parent_session=session)
  
  # Map view
  callModule(mod_map_view_server,
             "map_view_ui_1", 
             loadMap = datas$loadMap,
             loadJBrowse = datas$loadJBrowse,
             loadQTL = datas$loadQTL,
             parent_session=session)
  
  # QTL view
  callModule(mod_qtl_view_server,
             "qtl_view_ui_1", 
             loadMap = datas$loadMap,
             loadJBrowse = datas$loadJBrowse,
             loadQTL = datas$loadQTL)
  
  # Download
  callModule(mod_download_server,
             "download_ui_1")
}
