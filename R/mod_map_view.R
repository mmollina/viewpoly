#' map_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_map_view_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      verticalLayout(
        fluidRow(
          tags$h4(tags$b("View Map")), hr(),
          column(2,
                 selectInput(inputId = ns("group"), label = p("Linkage group"), choices = 1:15, selected = 1),
                 checkboxInput(ns("op"), label = "Show SNP names", value = TRUE)
          ),
          column(1,
                 h4("Legend"),
                 includeHTML(system.file("ext/include.html", package="viewpoly"))
          ),
          # column(2,
          #        h4("Number of SNPs per dosage"),
          #        verbatimTextOutput(ns("text1")),
          #        HTML("Rows: Parent 1<br>Columns: Parent 2")
          # ),
          # column(2,
          #        h4("Summary"),
          #        verbatimTextOutput(ns("text2"))
          # ),
          # column(4,
          #        h4("Notes"),
          #        includeHTML(system.file("ext/include2.html", package="viewpoly"))
          # )
        ), hr(),
        wellPanel(
          sliderInput(ns("range"), "Map range (cM)", min = 0, max = 300,
                      value = c(0, 20), step = 1), style = "padding: 6px;"
        ),
        plotOutput(ns("plot1"), height = "500px"), hr(),
        JBrowseROutput(ns("browserOutput"))
      )
    )
  )
}


#' map_view Server Functions
#'
#' @import JBrowseR
#'
#' @noRd 
mod_map_view_server <-  function(input, output, session, loadMap, loadJBrowse){
  ns <- session$ns
  
  observe({
    group_choices <- as.list(1:length(loadMap()$dp))
    names(group_choices) <- 1:length(loadMap()$dp)
    
    updateSelectInput(session, "group",
                      label="Linkage group",
                      choices = group_choices,
                      selected= group_choices[[1]])
  })

  output$plot1 <- renderPlot({
    draw_map_shiny(left.lim = input$range[1], 
                   right.lim = input$range[2], 
                   ch = input$group,
                   dp = loadMap()$dp,
                   dq = loadMap()$dq, 
                   maps = loadMap()$maps, 
                   ph.p = loadMap()$ph.p, 
                   ph.q = loadMap()$ph.q,
                   snp.names = input$op)
    
    max_updated = reactive({
      map_summary(left.lim = input$range[1], right.lim = input$range[2], ch = input$group, maps = loadMap()$maps, dp = loadMap()$dp, dq = loadMap()$dq)[[5]]
    })
    observeEvent(max_updated, {
      updateSliderInput(inputId = "range", max = max_updated())
    })
  })
  
  
  # output$text1 <- renderPrint({
  #   map_summary(left.lim = input$range[1],
  #               right.lim = input$range[2],
  #               ch = input$group,
  #               maps = loadMap()$maps,
  #               dp = loadMap()$dp,
  #               dq = loadMap()$dq)[[1]]
  # })
  # 
  # output$text2 <- renderPrint({
  #   map_summary(left.lim = input$range[1],
  #               right.lim = input$range[2],
  #               ch = input$group,
  #               maps = loadMap()$maps,
  #               dp = loadMap()$dp,
  #               dq = loadMap()$dq)[2:4]
  # })
  
  # if(loadJBrowse()$example == "bt_map"){
  #   assembly <- assembly(system.file("ext/NSP306_trifida_chr_v3.fa.gz", package = "viewpoly"), bgzip = TRUE)
  #   annotations_track <- track_feature(system.file("ext/NSP306_trifida_v3.hc.gene_models.gff3", package = "viewpoly"),
  #                                      assembly)
  #   mk.pos <- readRDS(system.file("ext/mk_pos.rds", package = "viewpoly"))
  #   # Add tracks
  #   # variants_track <- track_variant()
  #   # alignments_track <- track_alignments()
  # } else {
  #   assembly <- assembly(loadJBrowse()$fasta, bgzip = TRUE)
  #   annotations_track <- track_feature(loadJBrowse()$gff3, assembly)
  #   # Add tracks
  # }

  assembly <- assembly(
    system.file("ext/NSP306_trifida_chr_v3.fa.gz", package = "viewpoly"),
    bgzip = TRUE
  )
  
  # create configuration for a JB2 GFF FeatureTrack
  annotations_track <- track_feature(
    system.file("ext/NSP306_trifida_v3.hc.gene_models.gff3", package = "viewpoly"),
    assembly
  )
  
  # create the tracks array to pass to browser
  tracks <- tracks(annotations_track)
  
  default_session <- default_session(
    assembly,
    c(annotations_track)
  )
  
  # link the UI with the browser widget
  output$browserOutput <- renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly,
      # pass our tracks here
      tracks = tracks,
      location = "Chr01:1..100",
      defaultSession = default_session
    )
  )
  
  # create the tracks array to pass to browser
  # tracks <- tracks(annotations_track) # ,variants_track, alignments_track)

  # # set up the default session for the browser
  # default_session <- default_session(
  #   assembly,
  #   c(annotations_track)
  # )
  # 
  # # Select default window
  # mk.cM <- data.frame(mk= names(loadMap()$maps[[input$select]]), cM = loadMap()$maps[[input$select]])
  # mk.pos <- filter(mk.pos, chr == input$select)
  # mks <- merge(mk.pos, mk.cM, by = c("mk"))
  # mks <- mks[order(mks$cM),]
  # mks.range <- which(mks$cM >= input$range[1] &  mks$cM <= input$range[2])
  # mks.range.1 <- mks$pos[mks.range[1]]
  # mks.range.2 <- mks$pos[mks.range[length(mks.range)]]

  # link the UI with the browser widget
  # output$browserOutput <- renderJBrowseR(
  #   if(!is.null(assembly)){
  #     JBrowseR(
  #       "View",
  #       assembly = assembly,
  #       tracks = tracks
  #     )
  #   } else {
  #     cat("Genome information not provided in upload session.")
  #   }
  # )
  
}

## To be copied in the UI
# mod_map_view_ui("map_view_ui_1")

## To be copied in the server
# mod_map_view_server("map_view_ui_1")
