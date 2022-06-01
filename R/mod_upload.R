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
                 div(
                   actionButton(ns("goAbout"), "Go to About",icon("arrow-circle-left", verify_fa = FALSE), class = "btn btn-primary"), 
                   actionButton(ns("goQTL"), label = div("Go to QTL", icon("arrow-circle-right", verify_fa = FALSE)), class = "btn btn-primary")), br(),
                 div(style = "position:absolute;right:0em;",
                     actionButton(ns("reset_all"), "Reset all",icon("undo-alt", verify_fa = FALSE), class = "btn btn-danger"))
             ),
             tags$h2(tags$b("Input data")), br(),
             "Use this module to select an example dataset or to upload yours.", br(), br()
      ), br(),
      column(width = 12,
             fluidPage(
               box(id= ns("box_example"), width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status="primary", title = actionLink(inputId = ns("exampleID"), label = tags$b("Available example datasets")),
                   radioButtons(ns("example_map"), label = p("They contain the entire linkage map and QTL analysis but just a subset of individuals."), 
                                choices = c("Potato - Atlantic x B1829-5" = "tetra_map"),
                                selected = "tetra_map"), br(), br(), hr(),
                   tags$p("Access complete example datasets ", 
                          tags$a(href= "https://www.polyploids.org/input-tests","here"))
               )
             )
      ), br(),
      column(width = 12,
             fluidPage(
               box(id = ns("box_map"), width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = actionLink(inputId = ns("mapID"), label = tags$b("Upload linkage map files")),
                   div(style = "position:absolute;right:1em;",
                       actionBttn(ns("reset_map"), style = "jelly", color = "royal",  size = "sm", label = "reset", icon = icon("undo-alt", verify_fa = FALSE))
                   ), br(), br(), 
                   box(id = ns("box_mappoly"),width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE, status="primary",  title = actionLink(inputId = ns("mappolyID"), label = tags$b("Upload MAPpoly output")),
                       tags$p("Access further information about how to build a linkage maps with MAPpoly ", 
                              tags$a(href= "https://rpubs.com/mmollin/tetra_mappoly_vignette","here")), br(),
                       tags$p("Access a example code of how to obtain these inputs using MAPpoly functions ", 
                              tags$a(href= "https://cristianetaniguti.github.io/viewpoly_vignettes/VIEWpoly_tutorial.html#Upload_linkage_map_files","here")),
                       hr(),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_mappoly"), style = "jelly", color = "royal",  size = "sm", label = "submit MAPpoly", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
                       tags$p("Object of class `mappoly.map`."), 
                       fileInput(ns("mappoly_in"), label = h6("File: my_mappoly_list.RData"), multiple = F),
                   ),
                   box(id = ns("box_polymap"),width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary",title = actionLink(inputId = ns("polymapID"), label = tags$b("Upload polymapR output")),
                       tags$p("Access further information about how to build a linkage maps with polymapR ", 
                              tags$a(href= "https://cran.r-project.org/web/packages/polymapR/vignettes/Vignette_polymapR.html","here")), br(),
                       tags$p("Access a example code of how to obtain these inputs using polymapR functions ", 
                              tags$a(href= "https://cristianetaniguti.github.io/viewpoly_vignettes/VIEWpoly_tutorial.html#Upload_linkage_map_files","here")),
                       hr(),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_polymapR"), style = "jelly", color = "royal",  size = "sm", label = "submit polymapR", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
                       p("Indicates whether the genotype input is discrete or probabilistic."),
                       prettyRadioButtons(
                         inputId = ns("input.type"),
                         label = "Data type:", 
                         choices = c("discrete" = "discrete", "probabilistic" = "probabilistic"),
                         selected = "discrete",
                         inline = TRUE, 
                         status = "info",
                         fill = TRUE
                       ), br(),
                       p("Indicates the dataset specie ploidy."), 
                       prettyRadioButtons(
                         inputId = ns("ploidy"),
                         label = "Ploidy:", 
                         choices = c(4, 6),
                         selected = 4,
                         inline = TRUE, 
                         status = "info",
                         fill = TRUE
                       ), br(),
                       fileInput(ns("polymapR.dataset"), label = h6("File: polymapR.dataset.RData"), multiple = F),
                       fileInput(ns("polymapR.map"), label = h6("File: polymapR.map.RData"), multiple = F),
                   ),
                   box(id = ns("box_mapst"), width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE, status="primary",title = actionLink(inputId = ns("mapstID"), label = tags$b("Upload linkage map files with standard format (.csv, .tsv or .tsv.gz)")),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_map_custom"), style = "jelly", color = "royal",  size = "sm", label = "submit map custom", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
                       fileInput(ns("dosages"), label = h6("File: dosages.tsv"), multiple = F),
                       fileInput(ns("genetic_map"), label = h6("File: genetic_map.tsv"), multiple = F),
                       fileInput(ns("phases"), label = h6("File: phases.tsv"), multiple = F),
                       p("Upload here an TSV file with table with three columns: 1) marker ID; 2) genome position; 3) chromosome"),
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
               box(id = ns("box_qtl"),width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = actionLink(inputId = ns("qtlID"), label = tags$b("Upload QTL analysis files")),
                   div(style = "position:absolute;right:1em;",
                       actionBttn(ns("reset_qtl"), style = "jelly", color = "royal",  size = "sm", label = "reset", icon = icon("undo-alt", verify_fa = FALSE))
                   ), br(), br(), 
                   box(id= ns("box_qtlpoly"), width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = actionLink(inputId = ns("qtlpolyID"), label = tags$b("Upload QTLpoly output")),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_qtlpoly"), style = "jelly", color = "royal",  size = "sm", label = "submit QTLpoly", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
                       tags$p("Access further information about how to perform QTL analysis with QTLpoly ", 
                              tags$a(href= "https://guilherme-pereira.github.io/QTLpoly/1-tutorial","here")), br(),
                       tags$p("Access a example code of how to obtain these inputs using QTLpoly functions ", 
                              tags$a(href= "https://cristianetaniguti.github.io/viewpoly_vignettes/VIEWpoly_tutorial.html#Upload_QTL_analysis_files","here")),
                       hr(),
                       fileInput(ns("qtlpoly_data"), label = h6("File: QTLpoly_data.RData", br(), br(),"Object of class: qtlpoly.data"), multiple = F),
                       fileInput(ns("qtlpoly_remim.mod"), label = h6("File: QTLpoly_remim.mod.RData", br(), br(), "Object of class: qtlpoly.remim"), multiple = F),
                       fileInput(ns("qtlpoly_est.effects"), label = h6("File: QTLpoly_est.effects.RData", br(), br(),"Object of class: qtlpoly.effects"), multiple = F),
                       fileInput(ns("qtlpoly_fitted.mod"), label = h6("File: QTLpoly_fitted.mod.RData", br(), br(), "Object of class: qtlpoly.fitted"), multiple = F),
                   ),
                   box(id = ns("box_diaqtl"), width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = actionLink(inputId = ns("diaqtlID"), label = tags$b("Upload diaQTL output")),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_diaQTL"), style = "jelly", color = "royal",  size = "sm", label = "submit diaQTL", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
                       tags$p("Access further information about how to perform QTL analysis with diaQTL ", 
                              tags$a(href= "https://jendelman.github.io/diaQTL/diaQTL_Vignette.html","here")), br(),
                       tags$p("Access a example code of how to obtain these inputs using diaQTL functions ", 
                              tags$a(href= "https://cristianetaniguti.github.io/viewpoly_vignettes/VIEWpoly_tutorial.html#Upload_QTL_analysis_files","here")),
                       hr(),
                       fileInput(ns("diaQTL_scan1"), label = h6("File: diaQTL_scan1_list.RData"), multiple = F),
                       fileInput(ns("diaQTL_scan1.summaries"), label = h6("File: diaQTL_scan1.summaries_list.RData"), multiple = F),
                       fileInput(ns("diaQTL_BayesCI"), label = h6("File: diaQTL_BayesCI_list.RData"), multiple = F),
                       fileInput(ns("diaQTL_fitQTL"), label = h6("File: diaQTL_fitQTL_list.RData"), multiple = F),
                   ),
                   box(id = ns("box_polyqtl"),width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = actionLink(inputId = ns("polyqtlID"), label = tags$b("Upload polyqtlR output")),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_polyqtlR"), style = "jelly", color = "royal",  size = "sm", label = "submit polyqtlR", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
                       tags$p("Access further information about how to perform QTL analysis with polyqtlR ", 
                              tags$a(href= "https://cran.r-project.org/web/packages/polyqtlR/vignettes/polyqtlR_vignette.html","here")), br(),
                       tags$p("Access a example code of how to obtain these inputs using polyqtlR functions ", 
                              tags$a(href= "https://cristianetaniguti.github.io/viewpoly_vignettes/VIEWpoly_tutorial.html#Upload_QTL_analysis_files","here")),
                       hr(),
                       fileInput(ns("polyqtlR_effects"), label = h6("File: polyqtlR_effects.RData"), multiple = F), hr(),
                       fileInput(ns("polyqtlR_qtl_info"), label = h6("File: polyqtlR_qtl_info.RData"), multiple = F),
                       fileInput(ns("polyqtlR_QTLscan_list"), label = h6("File: polyqtlR_QTLscan_list.RData"), multiple = F),
                   ),
                   box(id = ns("box_qtlst"),width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE, status="primary", title = actionLink(inputId = ns("qtlstID"), label = tags$b("Upload QTL analysis results with standard format (.csv, .tsv or .tsv.gz)")),
                       div(style = "position:absolute;right:1em;",
                           actionBttn(ns("submit_qtl_custom"), style = "jelly", color = "royal",  size = "sm", label = "submit QTL custom", icon = icon("share-square", verify_fa = FALSE)), 
                       ), br(), br(),
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
               box(id = ns("box_genome"),width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="primary", title = actionLink(inputId = ns("genomeID"), label = tags$b("Upload Genome Browser files")),
                   tags$p("Access further information about the files expected in this section ", 
                          tags$a(href= "https://gmod.github.io/JBrowseR/articles/creating-urls.html","here")), br(),
                   
                   column(6,
                          tags$h5(tags$b("Upload genome information")),
                          box(
                            width = NULL, background = "red",
                            "Warning! The uploaded .fasta, .gff3, .vcf, .bam, .cram, .wig genome version must be the same one used to build the genetic map"
                          )
                   ), 
                   column(6,
                          div(style = "position:absolute;right:1em;",
                              actionBttn(ns("reset_genome"), style = "jelly", color = "royal",  size = "sm", label = "reset", icon = icon("undo-alt", verify_fa = FALSE)), br(),br(),
                              actionBttn(ns("submit_genome"), style = "jelly", color = "royal",  size = "sm", label = "submit", icon = icon("share-square", verify_fa = FALSE)), br(),br(),
                          ), 
                   ), br(), br(),
                   column(12,
                          br(), 
                          tags$h5(tags$b("Upload .fasta/.fasta.gz and .fasta.fai/.fasta.gz.fai,.fasta.gz.gzi file with assembly information. Using this option, a local HTTP server will be generated.")),
                          fileInput(ns("fasta"), label = h6("Files: genome_v2.fasta.gz, genome_v2.fasta.gz.fai, genome_v2.fasta.gz.gzi"), multiple = T), 
                          p("or"), 
                          tags$h5(tags$b("Add the URL of the hosted FASTA file location. The loading procedure is more efficient using this option.")),
                          textInput(ns("fasta_server"), label = h6("https://jbrowse.org/genomes/sars-cov2/fasta/sars-cov2.fa.gz"), value = NULL),
                          br(), hr(),
                          tags$h5(tags$b("Upload .gff3/.gff3.gz and .gff3.tbi/.gff3.gz.tbi file with annotation information")),
                          fileInput(ns("gff3"), label = h6("Files: genome_v2.gff3.gz, genome_v2.gff3.gz.tbi"), multiple = T), 
                          p("or"), 
                          tags$h5(tags$b("Add the URL of the hosted GFF3 file location. The loading procedure is more efficient using this option.")),
                          textInput(ns("gff3_server"), label = h6("https://jbrowse.org/genomes/sars-cov2/sars-cov2-annotations.sorted.gff.gz"), value = NULL),
                          br(), hr(),
                          tags$h5(tags$b("Upload VCF file with variants information")),
                          fileInput(ns("vcf"), label = h6("Files: markers.vcf, markers.vcf.tbi"), multiple = T), 
                          p("or"), 
                          tags$h5(tags$b("Add the URL of the hosted VCF file location. The loading procedure is more efficient using this option.")),
                          textInput(ns("vcf_server"), label = h6("https://some/path/file.vcf"), value = NULL),
                          br(), hr(),
                          tags$h5(tags$b("Upload .bam and .bam.bai or .cram and .cram.crai file with alignment information")),
                          fileInput(ns("align"), label = h6("Files: all_ind.bam, all_ind.bam.bai"), multiple = T), 
                          p("or"), 
                          tags$h5(tags$b("Add the URL of the hosted BAM or CRAM file location. The loading procedure is more efficient using this option.")),
                          textInput(ns("align_server"), label = h6("https://some/path/file.bam"), value = NULL),
                          br(), hr(),
                          tags$h5(tags$b("Upload .wig file with bigWig information")),
                          fileInput(ns("wig"), label = h6("File: data.wig"), multiple = F), 
                          p("or"), 
                          tags$h5(tags$b("Add the URL of the hosted WIG file location. The loading procedure is more efficient using this option.")),
                          textInput(ns("wig_server"), label = h6("https://some/path/file.wig"), value = NULL),
                   )
               )
             )
      ),
      column(width = 12,
             fluidPage(
               box(id = ns("box_viewpoly"),width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status="info", title = actionLink(inputId = ns("viewpolyID"), label = tags$b("Download VIEWpoly dataset")),
                   p("The uploaded data are converted to the viewpoly format. It keeps the map and the QTL information. Genome information is not stored."), br(),
                   textInput(ns("data.name"), label = p("Define the dataset name. Do not use spaces between words."), value = "dataset_name"), br(),
                   
                   downloadBttn(ns('export_viewpoly'), style = "gradient", color = "royal")
               )
             )
      ),
      column(width = 12,
             fluidPage(
               box(id = ns("box_viewpolyup"),width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status="info", title = actionLink(inputId = ns("viewpolyupID"), label = tags$b("Upload VIEWpoly dataset")),
                   column(8,
                          radioButtons(ns("viewpoly_env"), width = 500, label = "Check one of the availables datasets:", 
                                       choices = "There is no VIEWpoly object in your R environment. Load VIEWpoly object or convert other formats below.",
                                       selected =  "There is no VIEWpoly object in your R environment. Load VIEWpoly object or convert other formats below."), br(),
                   ),
                   column(4,
                          div(style = "position:absolute;right:1em;",
                              actionBttn(ns("reset_viewpoly"), style = "jelly", color = "royal",  size = "sm", label = "reset", icon = icon("undo-alt", verify_fa = FALSE)), br(), br(),
                              actionBttn(ns("submit_viewpoly"), style = "jelly", color = "royal",  size = "sm", label = "submit VIEWpoly file", icon = icon("share-square", verify_fa = FALSE))
                          )
                   ), 
                   column(12,
                          br(), br(), hr(),
                          p("Upload VIEWpoly RData file here:"), 
                          fileInput(ns("viewpoly_input"), label = h6("File: dataset_name.RData"), multiple = F)
                   )
               )
             )
      )
    )
  )
}

#' upload Server Functions
#'
#' @import vroom
#' @import shinyjs
#' @importFrom utils packageVersion
#' 
#' @noRd 
mod_upload_server <- function(input, output, session, parent_session){
  ns <- session$ns
  
  
  #Collapse boxes
  observeEvent(input$exampleID, {
    js$collapse(ns("box_example"))
  })
  
  observeEvent(input$mapID, {
    js$collapse(ns("box_map"))
  })
  
  observeEvent(input$mappolyID, {
    js$collapse(ns("box_mappoly"))
  })
  
  observeEvent(input$polymapID, {
    js$collapse(ns("box_polymap"))
  })
  
  observeEvent(input$mapstID, {
    js$collapse(ns("box_mapst"))
  })
  
  observeEvent(input$qtlID, {
    js$collapse(ns("box_qtl"))
  })
  
  observeEvent(input$qtlpolyID, {
    js$collapse(ns("box_qtlpoly"))
  })
  
  observeEvent(input$diaqtlID, {
    js$collapse(ns("box_diaqtl"))
  })
  
  observeEvent(input$polyqtlID, {
    js$collapse(ns("box_polyqtl"))
  })
  
  observeEvent(input$qtlstID, {
    js$collapse(ns("box_qtlst"))
  })
  
  observeEvent(input$genomeID, {
    js$collapse(ns("box_genome"))
  })
  
  observeEvent(input$viewpolyID, {
    js$collapse(ns("box_viewpoly"))
  })
  
  observeEvent(input$viewpolyupID, {
    js$collapse(ns("box_viewpolyup"))
  })

  # Check environment
  observe({
    Objs <- Filter(function(x) inherits(get(x), 'viewpoly' ), ls(envir = .GlobalEnv) )
    if(length(Objs) > 0){
      dataset_choices <- as.list(Objs)
      names(dataset_choices) <- Objs
      updateRadioButtons(session, "viewpoly_env",
                         label="Check one of the availables datasets:",
                         choices = dataset_choices,
                         selected= character(0))
    } else {
      updateRadioButtons(session, "viewpoly_env",
                         label="Check one of the availables datasets:",
                         choices = "There is no viewpoly object in your R environment. Load view viewpoly object or convert formats below",
                         selected= character(0))
    }
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
  
  observeEvent(input$goAbout, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "about")
  })
  
  # Reset buttons
  values <- reactiveValues(
    upload_state_map = 0,
    upload_state_mappoly = 0,
    upload_state_polymapR = 0,
    upload_state_map_custom = 0,
    upload_state_qtl = 0,
    upload_state_qtlpoly = 0,
    upload_state_diaQTL = 0,
    upload_state_polyqtlR = 0,
    upload_state_qtl_custom = 0,
    upload_state_genome = 0
  )
  
  observeEvent(input$reset_all, {
    values$upload_state_viewpoly <- 'reset'
    values$upload_state_map <- 'reset'
    values$upload_state_mappoly = 0
    values$upload_state_polymapR = 0
    values$upload_state_map_custom = 0
    values$upload_state_qtl <- 'reset'
    values$upload_state_qtlpoly = 0
    values$upload_state_diaQTL = 0
    values$upload_state_polyqtlR = 0
    values$upload_state_qtl_custom = 0
    values$upload_state_genome <- 'reset'
  })
  
  observeEvent(input$reset_viewpoly, {
    values$upload_state_viewpoly <- 'reset'
  })
  
  observeEvent(input$reset_map, {
    values$upload_state_map <- 'reset'
    values$upload_state_mappoly = 0
    values$upload_state_polymapR = 0
    values$upload_state_map_custom = 0
  })
  
  observeEvent(input$reset_qtl, {
    values$upload_state_qtl <- 'reset'
    values$upload_state_qtlpoly = 0
    values$upload_state_diaQTL = 0
    values$upload_state_polyqtlR = 0
    values$upload_state_qtl_custom = 0
  })
  
  observeEvent(input$reset_genome, {
    values$upload_state_genome <- 'reset'
  })
  
  observeEvent(input$submit_viewpoly, {
    values$upload_state_viewpoly <- 'uploaded'
  })
  
  observeEvent(input$submit_mappoly, {
    values$upload_state_mappoly <- 'uploaded'
    values$upload_state_map <- 0
  })
  
  observeEvent(input$submit_polymapR, {
    values$upload_state_polymapR <- 'uploaded'
    values$upload_state_map <- 0
  })
  
  observeEvent(input$submit_map_custom, {
    values$upload_state_map_custom <- 'uploaded'
    values$upload_state_map <- 0
  })
  
  observeEvent(input$submit_qtlpoly, {
    values$upload_state_qtlpoly <- 'uploaded'
    values$upload_state_qtl = 0
  })
  
  observeEvent(input$submit_diaQTL, {
    values$upload_state_diaQTL <- 'uploaded'
    values$upload_state_qtl = 0
  })
  
  observeEvent(input$submit_polyqtlR, {
    values$upload_state_polyqtlR <- 'uploaded'
    values$upload_state_qtl = 0
  })
  
  observeEvent(input$submit_qtl_custom, {
    values$upload_state_qtl_custom <- 'uploaded'
    values$upload_state_qtl = 0
  })
  
  observeEvent(input$submit_genome, {
    values$upload_state_genome <- 'uploaded'
  })
  
  input_map <- reactive({
    if (values$upload_state_map == 0 & 
        values$upload_state_mappoly == 0 & 
        values$upload_state_polymapR == 0 &
        values$upload_state_map_custom == 0) {
      return(NULL)
    } else if (values$upload_state_map == 'reset') {
      return(NULL)
    } else if(values$upload_state_mappoly == "uploaded"){
      validate(
        need(!is.null(input$mappoly_in), "Upload mappoly file before submit")
      )
      return(list(mappoly_in = input$mappoly_in))
    } else if(values$upload_state_polymapR == "uploaded"){
      validate(
        need(!is.null(input$polymapR.dataset), "Upload polymapR dataset file before submit"),
        need(!is.null(input$polymapR.map), "Upload polymapR map file before submit")
      )
      return(list(polymapR.dataset = input$polymapR.dataset,
                  polymapR.map = input$polymapR.map, 
                  input.type = input$input.type, 
                  ploidy = as.numeric(input$ploidy)))
    } else if(values$upload_state_map_custom == "uploaded"){
      validate(
        need(!is.null(input$dosages), "Upload custom dosages file before submit"),
        need(!is.null(input$phases), "Upload custom phases file before submit"),
        need(!is.null(input$genetic_map), "Upload custom genetic map file before submit")
      )
      return(list(dosages = input$dosages,
                  phases = input$phases,
                  genetic_map = input$genetic_map))
    }
  })
  
  input_qtl <- reactive({
    if (values$upload_state_qtl == 0 & 
        values$upload_state_qtlpoly == 0 & 
        values$upload_state_diaQTL == 0 &
        values$upload_state_polyqtlR == 0 &
        values$upload_state_qtl_custom == 0) {
      return(NULL)
    } else if (values$upload_state_qtl == 'reset') {
      return(NULL)
    } else if(values$upload_state_qtl_custom == "uploaded"){
      validate(
        need(!is.null(input$dosages), "Upload custom selected markers file before submit"),
        need(!is.null(input$phases), "Upload custom QTL info file before submit"),
        need(!is.null(input$blups), "Upload custom BLUPs file before submit"),
        need(!is.null(input$beta.hat), "Upload custom beta hat file before submit"),
        need(!is.null(input$profile), "Upload custom QTL profile file before submit"),
        need(!is.null(input$effects), "Upload custom effects file before submit"),
        need(!is.null(input$probs), "Upload custom genotype probabilities file before submit")
      )
      return(list(selected_mks = input$selected_mks,
                  qtl_info = input$qtl_info,
                  blups = input$blups,
                  beta.hat = input$beta.hat,
                  profile = input$profile,
                  effects = input$effects,
                  probs = input$probs))
    } else if(values$upload_state_qtlpoly == "uploaded"){
      validate(
        need(!is.null(input$qtlpoly_data), "Upload QTLpoly data file before submit"),
        need(!is.null(input$qtlpoly_remim.mod), "Upload QTLpoly remim.mod file before submit"),
        need(!is.null(input$qtlpoly_est.effects), "Upload QTLpoly estimated effects file before submit"),
        need(!is.null(input$qtlpoly_fitted.mod), "Upload QTLpoly fitted.mod file before submit")
      )
      return(list(
        qtlpoly_data = input$qtlpoly_data,
        qtlpoly_remim.mod = input$qtlpoly_remim.mod,
        qtlpoly_est.effects = input$qtlpoly_est.effects,
        qtlpoly_fitted.mod = input$qtlpoly_fitted.mod))
    } else if(values$upload_state_diaQTL == "uploaded"){
      validate(
        need(!is.null(input$qtlpoly_data), "Upload diaQTL data file before submit"),
        need(!is.null(input$qtlpoly_remim.mod), "Upload diaQTL scan1 file before submit"),
        need(!is.null(input$qtlpoly_est.effects), "Upload diaQTL scan1.summaries file before submit"),
        need(!is.null(input$qtlpoly_fitted.mod), "Upload diaQTL fitQTL file before submit"),
        need(!is.null(input$qtlpoly_fitted.mod), "Upload diaQTL BayesCI file before submit")
      )
      return(list(
        diaQTL_data = input$diaQTL_data,
        diaQTL_scan1 = input$diaQTL_scan1,
        diaQTL_scan1.summaries = input$diaQTL_scan1.summaries,
        diaQTL_fitQTL = input$diaQTL_fitQTL,
        diaQTL_BayesCI = input$diaQTL_BayesCI
      ))
    } else if(values$upload_state_polyqtlR == "uploaded"){
      validate(
        need(!is.null(input$qtlpoly_data), "Upload polyqtlR scan list file before submit"),
        need(!is.null(input$qtlpoly_remim.mod), "Upload polyqtlR QTL info file before submit"),
        need(!is.null(input$qtlpoly_est.effects), "Upload polyqtlR estimated effects file before submit")
      )
      return(list(
        polyqtlR_QTLscan_list = input$polyqtlR_QTLscan_list,
        polyqtlR_qtl_info = input$polyqtlR_qtl_info,
        polyqtlR_effects = input$polyqtlR_effects
      ))
    } 
  })
  
  input_genome  <- reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading fasta path..."))
      if (is.null(values$upload_state_genome)) {
        return(NULL)
      } else if (values$upload_state_genome == 'reset') {
        return(NULL)
      } else if(values$upload_state_genome == "uploaded"){
        validate(
          need(!is.null(input$fasta), "Upload reference genome (FASTA) file before submit.")
        )
        return(list(fasta = input$fasta,
                    fasta_server = input$fasta_server,
                    gff3 = input$gff3,
                    gff3_server = input$gff3_server,
                    vcf = input$vcf,
                    vcf_server = input$vcf_server,
                    align = input$align,
                    align_server = input$align_server,
                    wig = input$wig,
                    wig_server = input$wig_server))
      }
    })
  })
  
  # Wait system for the uploads
  loadExample = reactive({
    if(is.null(input_map()$dosages) & is.null(input_map()$phases) & is.null(input_map()$genetic_map) &
       is.null(input_map()$mappoly_in) &
       is.null(input_map()$polymapR.dataset) &
       is.null(input_map()$polymapR.map) &
       is.null(input_qtl()$selected_mks) & 
       is.null(input_qtl()$qtl_info) & 
       is.null(input_qtl()$blups) & 
       is.null(input_qtl()$beta.hat) & 
       is.null(input_qtl()$profile) & 
       is.null(input_qtl()$effects) & 
       is.null(input_qtl()$probs) &
       is.null(input_qtl()$qtlpoly_data) & 
       is.null(input_qtl()$qtlpoly_remim.mod) &
       is.null(input_qtl()$qtlpoly_est.effects) &
       is.null(input_qtl()$qtlpoly_fitted.mod) &
       is.null(input_qtl()$diaQTL_data) & 
       is.null(input_qtl()$diaQTL_scan1) &
       is.null(input_qtl()$diaQTL_scan1.summaries) &
       is.null(input_qtl()$diaQTL_fitQTL) &
       is.null(input_qtl()$diaQTL_BayesCI) &
       is.null(input_qtl()$polyqtlR_QTLscan_list) & 
       is.null(input_qtl()$polyqtlR_qtl_info) &
       is.null(input_qtl()$polyqtlR_effects) &
       is.null(input_genome()$fasta) &
       is.null(input_genome()$fasta_server) &
       is.null(input_genome()$gff3) &
       is.null(input_genome()$gff3_server) &
       is.null(input_genome()$vcf) &
       is.null(input_genome()$vcf_server) &
       is.null(input_genome()$align) &
       is.null(input_genome()$align_server) &
       is.null(input_genome()$wig) &
       is.null(input_genome()$wig_server) &
       is.null(input$viewpoly_input) &
       is.null(input$viewpoly_env))
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("Uploading example map data..."))
      prepare_examples(input$example_map)
    })
    else NULL
  })
  
  loadViewpoly = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading viewpoly file..."))
      if (is.null(values$upload_state_viewpoly)) {
        return(NULL)
      } else if (values$upload_state_viewpoly == 'reset') {
        return(NULL)
      } else if(values$upload_state_viewpoly == "uploaded"){
        validate(
          need(!is.null(input$viewpoly_input) | !is.null(input$viewpoly_env), 
               "Upload a viewpoly dataset or select one available in your R environment before submit.")
        )
        if(!is.null(input$viewpoly_input)){
          temp <- load(input$viewpoly_input$datapath)
          viewpoly.obj <- get(temp)
        } else if(!is.null(input$viewpoly_env)) {
          viewpoly.obj = get(input$viewpoly_env)
        }
        return(viewpoly.obj)
      } 
    })
  })
  
  loadMap_custom = reactive({
    if(!(is.null(input_map()$dosages) & is.null(input_map()$phases) & is.null(input_map()$genetic_map))){
      req(input_map()$dosages, input_map()$phases, input_map()$genetic_map)
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Uploading custom map data..."))
        prepare_map_custom_files(input_map()$dosages,
                                 input_map()$phases,
                                 input_map()$genetic_map)
      })
    } else NULL
  })
  
  loadMap_mappoly =  reactive({
    
    if(!is.null(input_map()$mappoly_in)){
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.3, detail = paste("Uploading MAPpoly data..."))
        prepare_MAPpoly(input_map()$mappoly_in)
      })
    } else NULL
  })
  
  loadMap_polymapR =  reactive({
    if(!(is.null(input_map()$polymapR.dataset) & 
         is.null(input_map()$polymapR.map))) {
      req(input_map()$polymapR.dataset, input_map()$polymapR.map)
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.1, detail = paste("Uploading polymapR data..."))
        prepare_polymapR(input_map()$polymapR.dataset, input_map()$polymapR.map, 
                         input$input.type, as.numeric(input$ploidy))
      })
    } else NULL
  })
  
  loadQTL_custom = reactive({
    if(!(is.null(input_qtl()$selected_mks) & 
         is.null(input_qtl()$qtl_info) & 
         is.null(input_qtl()$blups) & 
         is.null(input_qtl()$beta.hat) & 
         is.null(input_qtl()$profile) & 
         is.null(input_qtl()$effects) & 
         is.null(input_qtl()$probs))) {
      req(input_qtl()$selected_mks, input_qtl()$qtl_info, input_qtl()$blups,
          input_qtl()$beta.hat, input_qtl()$profile, input_qtl()$effects,
          input_qtl()$probs)
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Uploading custom QTL data..."))
        prepare_qtl_custom_files(input_qtl()$selected_mks,
                                 input_qtl()$qtl_info,
                                 input_qtl()$blups,
                                 input_qtl()$beta.hat,
                                 input_qtl()$profile,
                                 input_qtl()$effects,
                                 input_qtl()$probs)
      })
    } else NULL
  })
  
  loadQTL_qtlpoly = reactive({
    if(!(is.null(input_qtl()$qtlpoly_data) & 
         is.null(input_qtl()$qtlpoly_remim.mod) &
         is.null(input_qtl()$qtlpoly_est.effects) &
         is.null(input_qtl()$qtlpoly_fitted.mod))) {
      
      req(input_qtl()$qtlpoly_data, 
          input_qtl()$qtlpoly_remim.mod,
          input_qtl()$qtlpoly_est.effects,
          input_qtl()$qtlpoly_fitted.mod)
      
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.3, detail = paste("Uploading QTLpoly data..."))
        prepare_QTLpoly(input_qtl()$qtlpoly_data,
                        input_qtl()$qtlpoly_remim.mod,
                        input_qtl()$qtlpoly_est.effects,
                        input_qtl()$qtlpoly_fitted.mod)
      })
    } else NULL
  })
  
  loadQTL_diaQTL = reactive({
    if(!(is.null(input_qtl()$diaQTL_scan1) &
         is.null(input_qtl()$diaQTL_scan1.summaries) &
         is.null(input_qtl()$diaQTL_fitQTL) &
         is.null(input_qtl()$diaQTL_BayesCI))) {
      
      req(input_qtl()$diaQTL_scan1,
          input_qtl()$diaQTL_scan1.summaries,
          input_qtl()$diaQTL_fitQTL,
          input_qtl()$diaQTL_BayesCI)
      
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.3, detail = paste("Uploading diaQTL data..."))
        prepare_diaQTL(input_qtl()$diaQTL_scan1,
                       input_qtl()$diaQTL_scan1.summaries,
                       input_qtl()$diaQTL_fitQTL,
                       input_qtl()$diaQTL_BayesCI)
      })
    } else NULL
  })
  
  loadQTL_polyqtlR = reactive({
    if(!(is.null(input_qtl()$polyqtlR_QTLscan_list) & 
         is.null(input_qtl()$polyqtlR_qtl_info) &
         is.null(input_qtl()$polyqtlR_effects))) {
      
      req(input_qtl()$polyqtlR_QTLscan_list,
          input_qtl()$polyqtlR_qtl_info,
          input_qtl()$polyqtlR_effects)
      
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.3, detail = paste("Uploading polyqtlR data..."))
        prepare_polyqtlR(input_qtl()$polyqtlR_QTLscan_list,
                         input_qtl()$polyqtlR_qtl_info,
                         input_qtl()$polyqtlR_effects)
      })
    } else NULL
  })
  
  temp_dir <- reactive(tempdir())
  
  loadJBrowse_fasta = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading fasta path..."))
      if(!is.null(input_genome()$fasta) & !is.null(loadMap())){
        # keep fasta name
        for(i in 1:length(input_genome()$fasta$datapath)){
          file.rename(input_genome()$fasta$datapath[i], 
                      file.path(temp_dir(), input_genome()$fasta$name[i]))
        }
        file.path(temp_dir(), input_genome()$fasta$name[1]) 
      } else if(!is.null(input_genome()$fasta_server) & !is.null(loadMap())) {
        input_genome()$fasta_server
      } else if(!is.null(input_genome()$fasta) | !is.null(input_genome()$fasta_server)) {
        warning("Load map data first to use this feature.")
      } else NULL
    })
  })
  
  loadJBrowse_gff3 = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading gff3 path..."))
      if(!is.null(input_genome()$gff3)){
        for(i in 1:length(input_genome()$gff3$datapath)){
          file.rename(input_genome()$gff3$datapath[i], 
                      file.path(temp_dir(), input_genome()$gff3$name[i]))
        }
        file.path(temp_dir(), input_genome()$gff3$name[1]) 
      } else if(!is.null(input_genome()$gff3_server)) { 
        input_genome()$gff3_server
      } else NULL
    })
  })
  
  loadJBrowse_vcf = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading VCF path..."))
      if(!is.null(input_genome()$vcf)) {
        for(i in 1:length(input_genome()$vcf$datapath)){
          file.rename(input_genome()$vcf$datapath[i], 
                      file.path(temp_dir(), input_genome()$vcf$name[i]))
        }
        file.path(temp_dir(), input_genome()$vcf$name[1]) 
      } else if(!is.null(input_genome()$vcf_server)) {
        input_genome()$vcf_server
      } else NULL
    })
  })
  
  loadJBrowse_align = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading BAM or CRAM alignment data path..."))
      if(!is.null(input_genome()$align)) {
        for(i in 1:length(input_genome()$align$datapath)){
          file.rename(input_genome()$align$datapath[i], 
                      file.path(temp_dir(), input_genome()$align$name[i]))
        }
        file.path(temp_dir(), input_genome()$align$name[1]) 
      } else if(!is.null(input_genome()$align_server)) {
        input_genome()$align_server
      } else NULL
    })
  })
  
  loadJBrowse_wig = reactive({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.1, detail = paste("Uploading bigWig data path..."))
      if(!is.null(input_genome()$wig)) {
        for(i in 1:length(input_genome()$wig$datapath)){
          file.rename(input_genome()$wig$datapath[i], 
                      file.path(temp_dir(), input_genome()$wig$name[i]))
        }
        file.path(temp_dir(), input_genome()$wig$name[1]) 
      } else if(!is.null(input_genome()$wig_server)) {
        input_genome()$wig_server
      } else NULL
    })
  })
  
  loadMap = reactive({
    if(is.null(loadExample()) & 
       is.null(loadMap_custom()) & 
       is.null(loadMap_mappoly()) &
       is.null(loadMap_polymapR()) &
       is.null(loadViewpoly())){
      warning("Select one of the options in `upload` session")
      return(NULL)
    } else if(!is.null(loadViewpoly())){
      return(loadViewpoly()$map)
    } else if(!is.null(loadMap_custom())){
      return(loadMap_custom())
    } else if(!is.null(loadMap_mappoly())){
      return(loadMap_mappoly())
    } else if(!is.null(loadMap_polymapR())){
      return(loadMap_polymapR())
    } else if(!is.null(loadExample())){
      return(loadExample()$map)
    }
  })
  
  loadQTL = reactive({
    if(is.null(loadExample()) & 
       is.null(loadQTL_custom()) & 
       is.null(loadQTL_qtlpoly()) & 
       is.null(loadQTL_diaQTL()) &
       is.null(loadQTL_polyqtlR()) &
       is.null(loadViewpoly())){
      warning("Select one of the options in `upload` session")
      return(NULL)
    } else if(!is.null(loadViewpoly())){
      return(loadViewpoly()$qtl)
    } else if(!is.null(loadQTL_custom())){
      return(loadQTL_custom())
    } else if(!is.null(loadQTL_qtlpoly())){
      return(loadQTL_qtlpoly())
    } else if(!is.null(loadQTL_diaQTL())){
      return(loadQTL_diaQTL())
    } else if(!is.null(loadQTL_polyqtlR())){
      return(loadQTL_polyqtlR())
    } else if(!is.null(loadExample())){
      return(loadExample()$qtl)
    }
  })
  
  output$export_viewpoly <- downloadHandler(
    filename = function() {
      paste0("viewpoly.RData")
    },
    content = function(file) {
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.1, detail = paste("Saving viewpoly object..."))
        obj <- structure(list(map = loadMap(), 
                              qtl = loadQTL(), 
                              fasta = NULL, # It would save only the temporary path
                              gff3 = NULL, 
                              vcf = NULL,
                              align = NULL,
                              wig = NULL,
                              version = packageVersion("viewpoly")), 
                         class = "viewpoly")
        assign(input$data.name, obj)
        incProgress(0.5, detail = paste("Saving viewpoly object..."))
        save(list = input$data.name, file = file)
      })
    }
  )  
  
  return(list(loadMap = reactive(loadMap()), 
              loadQTL = reactive(loadQTL()), 
              loadJBrowse_fasta = reactive(loadJBrowse_fasta()), 
              loadJBrowse_gff3 = reactive(loadJBrowse_gff3()), 
              loadJBrowse_vcf = reactive(loadJBrowse_vcf()),
              loadJBrowse_align = reactive(loadJBrowse_align()),
              loadJBrowse_wig = reactive(loadJBrowse_wig()),
              loadExample = reactive(loadExample())
  ))
}

## To be copied in the UI
# mod_upload_ui("upload_ui_1")

## To be copied in the server
# mod_upload_server("upload_ui_1")
