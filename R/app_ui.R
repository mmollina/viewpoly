#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @importFrom shiny NS tagList 
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic 
    dashboardPage(
      dashboardHeader(disable = TRUE),
      dashboardSidebar(disable = TRUE),
      dashboardBody(
        fluidPage(
          navbarPage( 
            tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://www.polyploids.org/\"><img src=\"logo.png\" alt=\"alt\" style=\"float:right;width:150px;height:100px;padding-top:10px;padding-bottom:10px;\"> </a>`</div>');
    console.log(header)")
            ),
            title =  "ViewPoly", 
            id = "viewpoly",
            theme = shinythemes::shinytheme("flatly"),  # <--- Specify theme here
            tabPanel("About",
                     includeMarkdown(system.file("ext", "about.Rmd", package = "viewpoly"))
            ),
            tabPanel("Upload", value = "upload",
                     mod_upload_ui("upload_ui_1")),
            tabPanel("QTL", value = "qtl",
                     mod_qtl_view_ui("qtl_view_ui_1")),
            tabPanel("Genes", value = "genes",
                     mod_genes_view_ui("genes_view_ui_1")),
            tabPanel("Map", value = "map",
                     mod_map_view_ui("map_view_ui_1"))
            # tabPanel("Downloads", value = "download",
            #          mod_download_ui("download_ui_1"))
            # navbarMenu("More",
            #            tabPanel("Summary", "Under development..."),
            #            tabPanel("Table", "Under development...")
            #)
          )
        )
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

