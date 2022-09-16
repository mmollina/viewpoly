#' upload UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#' 
#' @importFrom shinyjs inlineCSS useShinyjs
#' @importFrom shiny NS tagList 
mod_upload_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 12,
             div(style = "position:absolute;right:1em;", 
                 div(
                   actionButton(ns("goAbout"), "Go to About",icon("arrow-circle-left", verify_fa = FALSE), class = "btn btn-primary"), 
                   actionButton(ns("goQTL"), label = div("Go to QTL", icon("arrow-circle-right", verify_fa = FALSE)), class = "btn btn-primary")), br(),
             ),
             tags$h2(tags$b("Input data")), br(),
             "Use this module to select an dataset.", br(), br()
      ), br(),
      column(width = 12,
             fluidPage(
               box(id= ns("box_example"), width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status="primary", title = actionLink(inputId = ns("exampleID"), label = tags$b("Available datasets")),
                   radioButtons(ns("example_map"), label = "Select dataset:", 
                                choices = "This will be updated",
                                selected = "This will be updated"), br(), br(), hr()
               )
             )
      )
    )
  )
}

#' upload Server Functions
#'
#' @import vroom
#' @importFrom shinyjs js
#' @importFrom utils packageVersion
#' 
#' @noRd 
mod_upload_server <- function(input, output, session, parent_session){
  ns <- session$ns
  
  #Collapse boxes
  observeEvent(input$exampleID, {
    js$collapse(ns("box_example"))
  })
  
  observeEvent(input$goQTL, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "qtl")
  })
  
  observeEvent(input$goAbout, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "about")
  })
  
  # load datas
  observe({
    files <- list.files(system.file("ext/my_viewpoly_objects/", package = "viewpoly"))
    ids <- read.csv(system.file("ext/info_data.csv", package = "viewpoly"), header = T)

    choices <- as.list(ids$ID)
    names(choices) <- ids$ID
    
    updateRadioButtons(session, "example_map", label = "Select dataset:", 
                       choices = choices,
                       selected = choices[[1]])
  })
  
  # Wait system for the uploads
  loadExample = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("Uploading example map data..."))
      example <- prepare_examples(input$example_map)
      return(example)
    })
  })
  
  loadJBrowse_fasta = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading fasta path..."))
      if(!is.null(loadExample()$fasta)){
        return(loadExample()$fasta)
      } else {
        return(NULL)
      }
    })
  })
  
  loadJBrowse_gff3 = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading gff3 path..."))
      if(!is.null(loadExample()$gff3)){
        return(loadExample()$gff3)
      } else 
        return(NULL)
    })
  })
  
  loadJBrowse_vcf = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading VCF path..."))
      if(!is.null(loadExample()$vcf_server)) {
        return(loadExample()$vcf_server)
      } else 
        return(NULL)
    })
  })
  
  loadJBrowse_align = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading BAM or CRAM alignment data path..."))
      if(!is.null(loadExample()$align_server)) {
        return(loadExample()$align_server)
      } else 
        return(NULL)
    })
  })
  
  loadJBrowse_wig = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading bigWig data path..."))
      if(!is.null(loadExample()$wig_server)) {
        return(loadExample()$wig_server)
      } else 
        return(NULL)
    })
  })
  
  loadMap = reactive({
    if(is.null(loadExample())){
      warning("Select one of the options in `upload` session")
      return(NULL)
    } else if(!is.null(loadExample())){
      return(loadExample()$map)
    }
  })
  
  loadQTL = reactive({
    if(is.null(loadExample())){
      warning("Select one of the options in `upload` session")
      return(NULL)
    } else if(!is.null(loadExample())){
      return(loadExample()$qtl)
    }
  })
  
  return(list(loadMap = reactive(loadMap()), 
              loadQTL = reactive(loadQTL()), 
              loadJBrowse_fasta = reactive(loadJBrowse_fasta()), 
              loadJBrowse_gff3 = reactive(loadJBrowse_gff3()), 
              loadJBrowse_vcf = reactive(loadJBrowse_vcf()),
              loadJBrowse_align = reactive(loadJBrowse_align()),
              loadJBrowse_wig = reactive(loadJBrowse_wig())))
}

## To be copied in the UI
# mod_upload_ui("upload_ui_1")

## To be copied in the server
# mod_upload_server("upload_ui_1")
