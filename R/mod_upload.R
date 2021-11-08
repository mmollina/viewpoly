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
    fluidRow(
      column(width = 12,
             div(style = "position:absolute;right:1em;", 
                 actionButton(ns("server_off"), "Exit",icon("times-circle"), class = "btn btn-danger"), br(), br(),
                 actionButton(ns("goMap"), "Next",icon("arrow-circle-right"), class = "btn btn-success")
             ),
             tags$h2(tags$b("Upload data")), br(),
             "Use this session to upload your data or select one of our examples. 
           Click `Next` to build graphics in the next sessions.", br(), br()
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = tags$h4(tags$b("View Example")),
                   selectInput(ns("example_map"), width = 500, label = "If you still don't have your own dataset, you can view one of our examples:", 
                               choices = list("Sweetpotato genetic map - Beauregard x Tanzania (BT)" = "bt_map"),
                               selected = "bt_map"),
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = tags$h4(tags$b("View Map")),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload MAPpoly output")),
                       fileInput(ns("mappoly_in"), label = h6("File: mappoly_map.RData"), multiple = F)
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE, title = tags$h5(tags$b("Upload map informations standard format (.tsv or .tsv.gz)")),
                       tags$div(class="alert alert-dismissible alert-warning",
                                tags$h4(class="alert-heading", "Warning!"),
                                tags$p(class="mb-0", "This is the most difficult way of upload your data. Good luck!")),
                       fileInput(ns("dosages"), label = h6("File: dosages.tsv"), multiple = F),
                       fileInput(ns("genetic_map"), label = h6("File: genetic_map.tsv"), multiple = F),
                       fileInput(ns("phases"), label = h6("File: phases.tsv"), multiple = F),
                       p("Upload here an RDS file with table with three columns: 1) marker ID; 2) genome position; 3) chromosome"),
                       fileInput(ns("mks.pos"), label = h6("File: marker information"), multiple = F),
                       "Check the input file formats with the example files:", br(),
                       radioButtons(ns("downloadType_map"), "", 
                                    choices = c("dosages.tsv" = "dosages",
                                                "genetic_map.tsv" = "genetic_map",
                                                "phases.tsv" = "phases"),
                                    inline = TRUE), br(), br(),
                       downloadButton(ns("downloadData_map"), "Download"), 
                   )
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title =   tags$h4(tags$b("View QTL")),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload QTLpoly output")),
                       fileInput(ns("qtlpoly_data"), label = h6("File: QTLpoly_data.RData"), multiple = F),
                       fileInput(ns("qtlpoly_remim.mod"), label = h6("File: QTLpoly_remim.mod.RData"), multiple = F),
                       fileInput(ns("qtlpoly_est.effects"), label = h6("File: QTLpoly_est.effects.RData"), multiple = F),
                       fileInput(ns("qtlpoly_fitted.mod"), label = h6("File: QTLpoly_fitted.mod.RData"), multiple = F)
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload QTL informations standard format (.tsv or .tsv.gz)")),
                       tags$div(class="alert alert-dismissible alert-warning",
                                tags$h4(class="alert-heading", "Warning!"),
                                tags$p(class="mb-0", "This is the most difficult way of upload your data. Good luck!")),
                       fileInput(ns("selected_mks"), label = h6("File: selected_mks.tsv"), multiple = F),
                       fileInput(ns("qtl_info"), label = h6("File: qtl_info.tsv"), multiple = F),
                       fileInput(ns("blups"), label = h6("File: blups.tsv"), multiple = F),
                       fileInput(ns("beta.hat"), label = h6("File: beta.hat.tsv"), multiple = F),
                       fileInput(ns("profile"), label = h6("File: profile.tsv"), multiple = F),
                       fileInput(ns("effects"), label = h6("File: effects.tsv"), multiple = F),
                       fileInput(ns("probs"), label = h6("File: probs.tsv"), multiple = F),
                       "Check the input format with the example file:", br(), br(),
                       radioButtons(ns("downloadType_qtl"), "", 
                                    choices = c("selected_mks.tsv" = "selected_mks",
                                                "qtl_info.tsv" = "qtl_info",
                                                "blups.tsv" = "blups",
                                                "beta.hat.tsv" = "beta.hat",
                                                "profile.tsv" = "profile",
                                                "effects.tsv" = "effects",
                                                "probs.tsv" = "probs"),
                                    inline = TRUE), br(), br(),
                       downloadButton(ns("downloadData_qtl"), "Download"), 
                   )
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = tags$h4(tags$b("View genome")),
                   tags$h5(tags$b("Upload genome information")),
                   p("Here you must upload the genome FASTA file compressed with bgzip, and the index files .fai and .gzi"),
                   tags$div(class="alert alert-dismissible alert-warning",
                            tags$h4(class="alert-heading", "Warning!"),
                            tags$p(class="mb-0", "the uploaded .fasta and .gff genome version should be the same one used to build the genetic map")),
                   fileInput(ns("fasta"), label = h6("File: genome_v2.fasta"), multiple = T),
               )
             )
      ),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = tags$h4(tags$b("View genes")),
                   tags$h5(tags$b("Upload .gff3 file with genes information")),
                   tags$div(class="alert alert-dismissible alert-warning",
                            tags$h4(class="alert-heading", "Warning!"),
                            tags$p(class="mb-0", "the uploaded .fasta and .gff genome version should be the same one used to build the genetic map")),
                   fileInput(ns("gff3"), label = h6("File: genome_v2.gff3"), multiple = T),
               )
             )
      )
    )
  )
}

#' upload Server Functions
#'
#' @import vroom
#' 
#' @noRd 
mod_upload_server <- function(input, output, session, parent_session){
  ns <- session$ns
  
  #  Trying to fix server issue
  observeEvent(input$server_off, {
    httpuv::stopAllServers()
  })
  
  # Format examples
  output$downloadData_map <- downloadHandler(
    filename = function() {
      paste0(input$downloadType_map, ".tsv")
    },
    content = function(file) {
      if(input$downloadType_map == "dosages") {
        filetemp <- vroom(system.file("ext/dosage.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_map == "phases") {
        filetemp <- vroom(system.file("ext/phases.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_map == "genetic_map") {
        filetemp <- vroom(system.file("ext/map.tsv.gz", package = "viewpoly"))
      } 
      vroom_write(filetemp, file = file)
    }
  )
  
  output$downloadData_qtl <- downloadHandler(
    filename = function() {
      paste0(input$downloadType_qtl, ".tsv")
    },
    content = function(file) {
      if(input$downloadType_qtl == "qtl_info") {
        filetemp <- vroom(system.file("ext/qtl_info.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_qtl == "blups") {
        filetemp <- vroom(system.file("ext/blups.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_qtl == "beta.hat") {
        filetemp <- vroom(system.file("ext/beta.hat.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_qtl == "profile.hat") {
        filetemp <- vroom(system.file("ext/profile.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_qtl == "effects.hat") {
        filetemp <- vroom(system.file("ext/effects.tsv.gz", package = "viewpoly"))
      } else if(input$downloadType_qtl == "probs") {
        filetemp <- vroom(system.file("ext/probs.tsv.gz", package = "viewpoly"))
      } 
      vroom_write(filetemp, file = file)
    }
  )
  
  observeEvent(input$goMap, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "map")
  })
  
  return(
    list(
      loadMap = reactive({read_Mapdata(input$mappoly_in,
                                       input$dosages,
                                       input$phases,
                                       input$genetic_map,
                                       input$example_map)}),
      
      loadJBrowse = reactive({list(fasta = input$fasta,
                                   gff3 = input$gff3,
                                   mks.pos = input$mks.pos,
                                   example = input$example_map)}),
      
      loadQTL = reactive({read_QTLdata(input$qtlpoly_data,
                                       input$qtlpoly_remim.mod,
                                       input$qtlpoly_est.effects,
                                       input$qtlpoly_fitted.mod,
                                       input$selected_mks,
                                       input$qtl_info,
                                       input$blups,
                                       input$beta.hat,
                                       input$profile,
                                       input$effects,
                                       input$probs,
                                       input$example_map)})
      
    )
  )
}

## To be copied in the UI
# mod_upload_ui("upload_ui_1")

## To be copied in the server
# mod_upload_server("upload_ui_1")
