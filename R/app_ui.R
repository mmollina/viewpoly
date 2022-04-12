#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import markdown
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
        # Lab colors
        tags$head(tags$style(HTML('
        .navbar-static-top {background-color: green;}
        
        .navbar-default .navbar-nav>.active>a {background-color: green;}
        
        .body {
            background-color: #22284c;
        }
                              
        .box.box-solid.box-primary>.box-header {
          color:#fff;
          background:#6c81c0
                            }
        
        .box.box-solid.box-primary{
        border-bottom-color:#6c81c0;
        border-left-color:#6c81c0;
        border-right-color:#6c81c0;
        border-top-color:#6c81c0;
        }
        
        .box.box-solid.box-info>.box-header {
        color:#fff;
        background:#22284c
        }
        
        .box.box-solid.box-info{
        border-bottom-color:#22284c;
        border-left-color:#22284c;
        border-right-color:#22284c;
        border-top-color:#22284c;
        }
        
        .box.box-solid.box-warning>.box-header {
        color:#fff;
        background:#a91021ff
        }
        
        .box.box-solid.box-warning{
        border-bottom-color:#a91021ff;
        border-left-color:#a91021ff;
        border-right-color:#a91021ff;
        border-top-color:#a91021ff;
        }
                              '))),
        tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://www.polyploids.org/\"><img src=\"www/logo_white.png\" alt=\"alt\" style=\"float:right;width:120px;height:80px;padding-top:10px;padding-bottom:10px;\"> </a>`</div>');
    console.log(header)")
        ),
        tags$head(tags$style(HTML('.navbar-static-top {background-color: #22284c;}',
                                  '.navbar-default .navbar-nav>.active>a {background-color: #22284c;}'))),
        fluidPage(
          navbarPage(
            title =  "VIEWpoly", 
            id = "viewpoly",
            theme = shinythemes::shinytheme("flatly"),  # <--- Specify theme here
            tabPanel("About",
                     includeMarkdown(system.file("ext", "about.Rmd", package = "viewpoly"))
            ),
            tabPanel("Input data", value = "upload",
                     mod_upload_ui("upload_ui_1")),
            tabPanel("QTL", value = "qtl",
                     mod_qtl_view_ui("qtl_view_ui_1")),
            tabPanel("Genome", value = "genes",
                     mod_genes_view_ui("genes_view_ui_1")),
            tabPanel("Map", value = "map",
                     mod_map_view_ui("map_view_ui_1"))
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

