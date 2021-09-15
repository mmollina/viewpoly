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
        actionButton(ns("create_server"), "Create local server",icon("refresh")), br(),
        actionButton(ns("server_off"), "Turn off local server",icon("refresh")), br(),
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
mod_map_view_server <- function(input, output, session, loadMap, loadJBrowse){
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
  
  button <- eventReactive(input$create_server, {
    
    if(!is.null(loadJBrowse()$fasta)){
      cat("genome")
      str(loadJBrowse()$fasta)
      cat("gff")
      str(loadJBrowse()$gff3)
      
      server_dir <- tempdir()
      
      path.fa <- paste0(server_dir, "/", loadJBrowse()$fasta$name[1])
      path.fai <- paste0(server_dir, "/", loadJBrowse()$fasta$name[2])
      path.gzi <- paste0(server_dir, "/", loadJBrowse()$fasta$name[3])
      path.gff <- paste0(server_dir, "/", loadJBrowse()$gff$name[1])
      path.tbi <- paste0(server_dir, "/", loadJBrowse()$gff$name[2])
      
      file.rename(loadJBrowse()$fasta$datapath[1], path.fa)
      file.rename(loadJBrowse()$fasta$datapath[2], path.fai)
      file.rename(loadJBrowse()$fasta$datapath[3], path.gzi)
      file.rename(loadJBrowse()$gff$datapath[1], path.gff)
      file.rename(loadJBrowse()$gff$datapath[2], path.tbi)
      cat("path")
      print(path.fa)
      
      mk.pos <- readRDS(loadJBrowse()$mks.pos$datapath)
      
    } else if(loadJBrowse()$example == "bt_map"){
      path.fa <- system.file("ext/Trifida.Chr01.fa.gz", package = "viewpoly")
      path.gff <- system.file("ext/Trifida.Chr01.sorted.gff3.gz", package = "viewpoly")
      mk.pos <- readRDS(system.file("ext/mk_pos.rds", package = "viewpoly"))
      # Add other tracks
      # variants_track <- track_variant()
      # alignments_track <- track_alignments()
    } 
    
    if(exists("data_server"))
      data_server$stop_server()
    
    data_server <- serve_data(dirname(path.fa), port = 5000)
    
    list(path.fa, path.gff, data_server, mk.pos)
  })
  
  #  Trying to fix server issue
  observeEvent(input$server_off, {
    httpuv::stopAllServers()
  })
  
  # link the UI with the browser widget
  output$browserOutput <- renderJBrowseR({
    
    assembly <- assembly(
      paste0("http://127.0.0.1:5000/", basename(button()[[1]])), 
      bgzip = TRUE
    )
    
    # create configuration for a JB2 GFF FeatureTrack
    annotations_track <- track_feature(
      paste0("http://127.0.0.1:5000/", basename(button()[[2]])), 
      assembly
    )
    
    # create the tracks array to pass to browser
    tracks <- tracks(annotations_track)
    
    # # Select default window
    group <- as.numeric(input$group)
    mk.cM <- data.frame(mk= names(loadMap()$maps[[group]]), cM = loadMap()$maps[[group]])
    mk.pos <- filter(button()[[4]], chr == group)
    mks <- merge(mk.pos, mk.cM, by = c("mk"))
    mks <- mks[order(mks$cM),]
    mks.range <- which(mks$cM >= input$range[1] &  mks$cM <= input$range[2])
    mks.range.1 <- mks$pos[mks.range[1]]
    mks.range.2 <- mks$pos[mks.range[length(mks.range)]]
    
    default_session <- default_session(
      assembly,
      c(annotations_track)
    )
    
    JBrowseR(
      "View",
      assembly = assembly,
      # pass our tracks here
      tracks = tracks,
      location = paste0("Chr01:", mks.range.1,"..",mks.range.2),
      defaultSession = default_session
    )
  })
}

## To be copied in the UI
# mod_map_view_ui("map_view_ui_1")

## To be copied in the server
# mod_map_view_server("map_view_ui_1")
