#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shiny NS tagList 
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic 
    navbarPage("VIEWpoly",id = "viewpoly",
               theme = shinythemes::shinytheme("flatly"),  # <--- Specify theme here
               tabPanel("About",
                        includeMarkdown(system.file("ext", "about.Rmd", package = "viewpoly"))
                        ),
               tabPanel("Upload data", value = "upload",
                        mod_upload_ui("upload_ui_1")
                        ),
               tabPanel("Map", value = "map",
                        mod_map_view_ui("map_view_ui_1")),
               tabPanel("Download infos", value = "download",
                        mod_download_ui("download_ui_1")),
               navbarMenu("More",
                          tabPanel("Summary", "Under development..."),
                          tabPanel("Table", "Under development...")
               )
    )
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
  
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'viewpoly'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

