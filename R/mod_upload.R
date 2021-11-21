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
                 actionButton(ns("goQTL"), "Next",icon("arrow-circle-right"), class = "btn btn-success")
             ),
             tags$h2(tags$b("Upload data")), br(),
             "Use this session to select the available datasets or to upload your data. 
           Click `Next` to build graphics in the next sessions.", br(), br()
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status="primary", title = tags$h4(tags$b("Available datasets")),
                   radioButtons(ns("example_map"), width = 500, label = "Check one of the availables datasets:", 
                                choices = "tetra_map",
                                selected = "tetra_map"), br(),
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = tags$h4(tags$b("Upload View Map files")),
                   div(style = "position:absolute;right:1em;",
                       actionBttn(ns("reset_map"), style = "unite", label = "reset", icon = icon("redo"))
                   ), br(), br(), 
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload MAPpoly output")),
                       fileInput(ns("mappoly_in"), label = h6("File: mappoly_map.RData"), multiple = F)
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload polymapR output")),
                       fileInput(ns("polymapR.dataset"), label = h6("File: polymapR.dataset.RData"), multiple = F),
                       fileInput(ns("polymapR.map"), label = h6("File: polymapR.map.RData"), multiple = F)
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE, title = tags$h5(tags$b("Upload map informations standard format (.tsv or .tsv.gz)")),
                       box(
                         width = NULL, background = "red",
                         "This is the most difficult way of upload your data. Good luck!"
                       ),
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
                   ),
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title =   tags$h4(tags$b("Upload View QTL files")),
                   div(style = "position:absolute;right:1em;",
                       actionBttn(ns("reset_qtl"), style = "unite", label = "reset", icon = icon("redo"))
                   ), br(), br(), 
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload QTLpoly output")),
                       fileInput(ns("qtlpoly_data"), label = h6("File: QTLpoly_data.RData"), multiple = F),
                       fileInput(ns("qtlpoly_remim.mod"), label = h6("File: QTLpoly_remim.mod.RData"), multiple = F),
                       fileInput(ns("qtlpoly_est.effects"), label = h6("File: QTLpoly_est.effects.RData"), multiple = F),
                       fileInput(ns("qtlpoly_fitted.mod"), label = h6("File: QTLpoly_fitted.mod.RData"), multiple = F)
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload diaQTL output")),
                       fileInput(ns("diaQTL_scan1"), label = h6("File: diaQTL_scan1_list.RData"), multiple = F),
                       fileInput(ns("diaQTL_scan1.summaries"), label = h6("File: diaQTL_scan1.summaries_list.RData"), multiple = F),
                       fileInput(ns("diaQTL_fitQTL"), label = h6("File: diaQTL_fitQTL_list.RData"), multiple = F),
                       fileInput(ns("diaQTL_BayesCI"), label = h6("File: diaQTL_BayesCI_list.RData"), multiple = F)
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload polyqtlR output")),
                       fileInput(ns("polyqtlR_QTLscan_list"), label = h6("File: polyqtlR_QTLscan_list.RData"), multiple = F),
                       fileInput(ns("polyqtlR_IBD"), label = h6("File: polyqtlR_IBD.RData"), multiple = F),
                       fileInput(ns("polyqtlR_phenotypes"), label = h6("File: polyqtl_phenotypes.RData"), multiple = F),
                   ),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  title = tags$h5(tags$b("Upload QTL informations standard format (.tsv or .tsv.gz)")),
                       box(
                         width = NULL, background = "red",
                         "This is the most difficult way of upload your data. Good luck!"
                       ),
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
                   )
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = tags$h4(tags$b("View JBrowse")),
                   div(style = "position:absolute;right:1em;",
                       actionBttn(ns("reset_genome"), style = "unite", label = "reset", icon = icon("redo"))
                   ), br(), br(), 
                   tags$h5(tags$b("Upload genome information")),
                   p("Here you must upload the genome FASTA file compressed with bgzip, and the index files .fai and .gzi"),
                   box(
                     width = NULL, background = "red",
                     "Warning! The uploaded .fasta and .gff genome version should be the same one used to build the genetic map"
                   ),
                   fileInput(ns("fasta"), label = h6("File: genome_v2.fasta"), multiple = T),
                   tags$h5(tags$b("Upload .gff3 file with genes information")),
                   fileInput(ns("gff3"), label = h6("File: genome_v2.gff3"), multiple = T),
                   tags$h5(tags$b("Upload VCF file with genes information")),
                   box(
                     width = NULL, background = "red",
                     "Warning! The uploaded VCF file should be the same one used to build the genetic map"
                   ),
                   fileInput(ns("vcf"), label = h6("File: markers.vcf"), multiple = T),
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
  
  # Check environment
  observe({
    Objs <- Filter(function(x) inherits(get(x), 'viewpoly' ), ls(envir = .GlobalEnv) )
    dataset_choices <- as.list(c("tetra_map","hex_map", Objs))
    names(dataset_choices) <- c("Potato genetic map - Atlantic x B1829-5",
                                "Sweetpotato genetic map - Beauregard x Tanzania (BT)",
                                Objs)
    print(dataset_choices)
    updateRadioButtons(session, "example_map",
                       label="Check one of the availables datasets:",
                       choices = dataset_choices,
                       selected= dataset_choices[[1]])
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
  
  observeEvent(input$goQTL, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "qtl")
  })
  
  # Reset buttons
  observeEvent(input$reset_map,{
    input$dosages = input$phases = input$genetic_map = NULL
    input$mappoly_in = input$selected_mks = NULL
  })
  
  observeEvent(input$reset_qtl,{
    input$qtl_info = NULL
    input$blups = input$beta.hat = input$profile = input$effects = NULL
    input$probs = input$qtlpoly_data = input$qtlpoly_remim.mod = NULL
    input$qtlpoly_est.effects = input$qtlpoly_fitted.mod = NULL
    input$diaQTL_data = input$diaQTL_scan1 = NULL 
    input$diaQTL_scan1.summaries = input$diaQTL_fitQTL = NULL
    input$diaQTL_BayesCI = input$polyqtlR_QTLscan_list = NULL
    input$polyqtlR_IBD = input$polyqtlR_phenotypes = NULL
  })
  
  observeEvent(input$reset_genome,{
    input$fasta = NULL 
    input$gff3 = input$vcf = NULL    
  })
  
  return(list(
    
    loadExample = reactive({
      if(is.null(input$dosages) & is.null(input$phases) & is.null(input$genetic_map) &
         is.null(input$mappoly_in) &
         is.null(input$selected_mks) & 
         is.null(input$qtl_info) & 
         is.null(input$blups) & 
         is.null(input$beta.hat) & 
         is.null(input$profile) & 
         is.null(input$effects) & 
         is.null(input$probs) &
         is.null(input$qtlpoly_data) & 
         is.null(input$qtlpoly_remim.mod) &
         is.null(input$qtlpoly_est.effects) &
         is.null(input$qtlpoly_fitted.mod) &
         is.null(input$diaQTL_data) & 
         is.null(input$diaQTL_scan1) &
         is.null(input$diaQTL_scan1.summaries) &
         is.null(input$diaQTL_fitQTL) &
         is.null(input$diaQTL_BayesCI) &
         is.null(input$polyqtlR_QTLscan_list) & 
         is.null(input$polyqtlR_IBD) &
         is.null(input$polyqtlR_phenotypes) &
         is.null(input$fasta) &
         is.null(input$gff3) &
         is.null(input$vcf))
      prepare_examples(input$example_map, env.obj = get(input$example_map))
      else NULL
    }),
    
    loadMap_custom = reactive({
      if(!(is.null(input$dosages) & is.null(input$phases) & is.null(input$genetic_map))){
        req(input$dosages, input$phases, input$genetic_map)
        prepare_map_custom_files(input$dosages,
                                 input$phases,
                                 input$genetic_map)
      } else NULL
    }),
    
    loadMap_mappoly =  reactive({
      if(!is.null(input$mappoly_in))
        prepare_MAPpoly(input$mappoly_in)
      else NULL
    }),
    
    loadQTL_custom = reactive({
      if(!(is.null(input$selected_mks) & 
           is.null(input$qtl_info) & 
           is.null(input$blups) & 
           is.null(input$beta.hat) & 
           is.null(input$profile) & 
           is.null(input$effects) & 
           is.null(input$probs))) {
        req(input$selected_mks, input$qtl_info, input$blups,
            input$beta.hat, input$profile, input$effects,
            input$probs)
        prepare_qtl_custom_files(input$selected_mks,
                                 input$qtl_info,
                                 input$blups,
                                 input$beta.hat,
                                 input$profile,
                                 input$effects,
                                 input$probs)
      } else NULL
    }),
    
    loadQTL_qtlpoly = reactive({
      if(!(is.null(input$qtlpoly_data) & 
           is.null(input$qtlpoly_remim.mod) &
           is.null(input$qtlpoly_est.effects) &
           is.null(input$qtlpoly_fitted.mod))) {
        
        req(input$qtlpoly_data, 
            input$qtlpoly_remim.mod,
            input$qtlpoly_est.effects,
            input$qtlpoly_fitted.mod)
        
        prepare_QTLpoly(input$qtlpoly_data,
                        input$qtlpoly_remim.mod,
                        input$qtlpoly_est.effects,
                        input$qtlpoly_fitted.mod)
      } else NULL
    }),
    
    loadQTL_diaQTL = reactive({
      if(!(is.null(input$diaQTL_scan1) &
           is.null(input$diaQTL_scan1.summaries) &
           is.null(input$diaQTL_fitQTL) &
           is.null(input$diaQTL_BayesCI))) {
        
        req(input$diaQTL_scan1,
            input$diaQTL_scan1.summaries,
            input$diaQTL_fitQTL,
            input$diaQTL_BayesCI)
        
        prepare_diaQTL(input$diaQTL_scan1,
                       input$diaQTL_scan1.summaries,
                       input$diaQTL_fitQTL,
                       input$diaQTL_BayesCI)
      } else NULL
    }),
    
    loadQTL_polyqtlR = reactive({
      if(!(is.null(input$polyqtlR_QTLscan_list) & 
           is.null(input$polyqtlR_IBD) &
           is.null(input$polyqtlR_phenotypes))) {
        
        req(input$polyqtlR_QTLscan_list,
            input$polyqtlR_IBD,
            input$polyqtlR_phenotypes)
        
        prepare_polyqtlR(input$polyqtlR_QTLscan_list,
                         input$polyqtlR_IBD,
                         input$polyqtlR_phenotypes)
      } else NULL
    }),
    
    loadJBrowse_fasta = reactive({
      if(!is.null(input$fasta))
        input$fasta$datapath
      else NULL
    }),
    
    loadJBrowse_gff3 = reactive({
      if(!is.null(input$gff3))
        input$gff3$datapath
      else NULL
    }),
    
    loadJBrowse_vcf = reactive({
      if(!is.null(input$vcf))
        input$vcf$datapath
      else NULL
    })
  )
  )
}

## To be copied in the UI
# mod_upload_ui("upload_ui_1")

## To be copied in the server
# mod_upload_server("upload_ui_1")
