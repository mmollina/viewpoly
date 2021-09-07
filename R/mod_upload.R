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
    column(width = 6,
           fluidPage(
             tags$h4(tags$b("View Map")), hr(),
             tags$h5(tags$b("Upload mappoly:::prepare_map list:")),
             fileInput(ns("mappoly_in"), label = h6("File: mappoly_map.RData"), multiple = F),
             tags$h5(tags$b("Or upload the markers dosages, genetic map and phases information.")), br(),
             fileInput(ns("dosages"), label = h6("File: dosages.csv"), multiple = F),
             fileInput(ns("genetic_map"), label = h6("File: genetic_map.csv"), multiple = F),
             fileInput(ns("phases"), label = h6("File: phases.csv"), multiple = F),
             "Check the proper input file formats with the example files:", br(),
             radioButtons(ns("downloadType"), "", 
                          choices = c("dosages.csv" = "dosages",
                                      "genetic_map.csv" = "genetic_map",
                                      "phases.csv" = "phases"),
                          inline = TRUE),
             downloadButton(ns("downloadData"), "Download"), hr(),
             
             "If you still don't have your own dataset, you can view one of or examples:",
             selectInput(ns("example_map"), label = h4(tags$b("Examples")), 
                         choices = list("Sweetpotato genetic map - Beauregard x Tanzania (BT)" = "bt_map",
                                        "None" = "none"),
                         selected = "none"),
           )
    ),
    # column(width = 6,
    #        fluidPage(
    #          tags$h4(tags$b("Upload QTLpoly final object:")),
    #          fileInput(ns("qtlpoly_in"), label = h6("File: QTLpoly.RData"), multiple = F),
    #          "Or upload the markers dosages, genetic map and phases information.", br(),
    #          fileInput(ns("dosages"), label = h6("File: dosages.csv"), multiple = F),
    #          fileInput(ns("genetic_map"), label = h6("File: genetic_map.csv"), multiple = F),
    #          fileInput(ns("phases"), label = h6("File: phases.csv"), multiple = F),
    #          "Check the proper input formats with the example files:",
    #          radioButtons("downloadType", "Download Type", 
    #                       choices = c("dosages.csv" = "dosages",
    #                                   "genetic_map.csv" = "genetic_map",
    #                                   "phases.csv" = "phases"),
    #                       inline = TRUE),
    #          downloadButton("downloadData", "Download")
    #        ),
    #        
    #        fluidPage(
    #          "If you still don't have your own dataset, you can view one of or examples:",
    #          selectInput(ns("example_qtl"), label = h4(tags$b("Examples")), 
    #                      choices = list("Sweetpotato genetic qtls - Beauregard x Tanzania (BT)" = "bt_qtl"),
    #                      selected = "bt_qtl"),
    #        )
    #        
    # ),
    column(width = 6,
           fluidPage(
             tags$h4(tags$b("View genes")), hr(),
             tags$h5(tags$b("Upload .gff3 file with genes information")),
             tags$div(class="alert alert-dismissible alert-warning",
                      tags$h4(class="alert-heading", "Warning!"),
                      tags$p(class="mb-0", "the .gff3 genome version should be the same one used to build the genetic map")
             ),
             fileInput(ns("gff3"), label = h6("File: genome_v2.gff3"), multiple = F),
             
             
             "If you still don't have your own dataset, you can view one of or examples:",
             selectInput(ns("example_gene"), label = h4(tags$b("Examples")), 
                         choices = list("Sweetpotato genetic qtls - Beauregard x Tanzania (BT)" = "bt_gene"),
                         selected = "bt_gene"), br(), br(),
             hr(),
             actionButton(ns("goMap"), "Next",icon("arrow-circle-right", class = ))
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
    
    # Upload map info
    # loadMap <- reactive({prepare_Mapdata(input$mappoly_in,
    #                                      input$dosages,
    #                                      input$phases,
    #                                      input$genetic_map,
    #                                      input$example_map)})
    # 
    # loadQtl <- reactive({prepare_Qtldata(input$qtlpoly_in,
    #                                      input$example_qtl)})
    # 
    # loadGene <- reactive({prepare_Genedata(input$gff3,
    #                                      input$example_gene)})
    
    # Next button
    
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
                                            input$example_map)})
      )
    )
}

## To be copied in the UI
# mod_upload_ui("upload_ui_1")

## To be copied in the server
# mod_upload_server("upload_ui_1")
