#' upload UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_upload_ui <- function(id){
  ns <- NS(id)
  tagList(
    column(width = 12,
           tags$h3(tags$b("Upload data")), 
           div(style = "position:absolute;right:1em;", 
               actionButton(ns("goMap"), "Next",icon("arrow-circle-right", class = ))
           ),
           "Use this session to upload your data or select one of our examples. 
           Click `Next` to build graphics in the next sessions.", hr(),
    ), 
    column(width = 12,
           tags$h4(tags$b("View Example")), hr(),
           selectInput(ns("example_map"), width = 500, label = "If you still don't have your own dataset, you can view one of or examples:", 
                       choices = list("Sweetpotato genetic map - Beauregard x Tanzania (BT)" = "bt_map"),
                       selected = "bt_map"), hr(),
           
    ),
    column(width = 6,
           fluidPage(
             tags$h4(tags$b("View Map")), hr(),
             tags$h5(tags$b("Upload mappoly:::prepare_map list:")),
             fileInput(ns("mappoly_in"), label = h6("File: mappoly_map.RData"), multiple = F),
             tags$h5(tags$b("Or upload the markers dosages, genetic map and phases information.")), br(),
             fileInput(ns("dosages"), label = h6("File: dosages.csv"), multiple = F),
             fileInput(ns("genetic_map"), label = h6("File: genetic_map.csv"), multiple = F),
             fileInput(ns("phases"), label = h6("File: phases.csv"), multiple = F),
             "Check the input file formats with the example files:", br(),
             radioButtons(ns("downloadType"), "", 
                          choices = c("dosages.csv" = "dosages",
                                      "genetic_map.csv" = "genetic_map",
                                      "phases.csv" = "phases"),
                          inline = TRUE),
             downloadButton(ns("downloadData"), "Download"), hr()
           )
    ),
    column(width = 6,
           fluidPage(
             tags$h4(tags$b("View QTLs")), hr(),
             tags$h5(tags$b("Upload QTLpoly final object:")),
             fileInput(ns("qtlpoly_in"), label = h6("File: QTLpoly.RData"), multiple = F),
             "Or upload the markers dosages, genetic map and phases information.", br(),
             fileInput(ns("dosages"), label = h6("File: dosages.csv"), multiple = F),
             fileInput(ns("genetic_map"), label = h6("File: genetic_map.csv"), multiple = F),
             fileInput(ns("phases"), label = h6("File: phases.csv"), multiple = F),
             "Check the input format with the example file:", br(), br(),
             downloadButton("downloadQTL", "Download QTL input example file"), hr()
           ),
    ),
    column(width = 12,
           tags$div(class="alert alert-dismissible alert-warning",
                    tags$h4(class="alert-heading", "Warning!"),
                    tags$p(class="mb-0", "the uploaded .fasta and .gff genome version should be the same one used to build the genetic map")
           )
    ),
    column(width = 6,
           fluidPage(
             tags$h4(tags$b("View genomes")), hr(),
             tags$h5(tags$b("Upload .fasta file with genome information")),
             fileInput(ns("fasta"), label = h6("File: genome_v2.fasta"), multiple = F),
           )
    ),
    column(width = 6,
           fluidPage(
             tags$h4(tags$b("View genes")), hr(),
             tags$h5(tags$b("Upload .gff3 file with genes information")),
             fileInput(ns("gff3"), label = h6("File: genome_v2.gff3"), multiple = F),
           )
    )
  )
}

#' upload Server Functions
#'
#' @noRd 
mod_upload_server <- function(input, output, session, parent_session){
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
  
  observeEvent(input$goMap, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "map")
  })
  
  return(
    list(
      loadMap = reactive({prepare_Mapdata(input$mappoly_in,
                                          input$dosages,
                                          input$phases,
                                          input$genetic_map,
                                          input$example_map)}),
      
      loadJBrowse = reactive({list(fasta = input$fasta,
                                   gff3 = input$gff3,
                                   example = input$example_map)})
    )
  )
}

## To be copied in the UI
# mod_upload_ui("upload_ui_1")

## To be copied in the server
# mod_upload_server("upload_ui_1")
