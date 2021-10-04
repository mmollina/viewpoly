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
          inlineCSS(".form-group {margin-bottom: 0;}
                                .irs-with-grid {bottom: 0px;}
                                .irs-grid {height: 13px;}
                                .irs-grid-text {height: 0px;}
                                "
          ),
          column(width = 12,
                 div(style = "position:absolute;right:1em;", 
                     actionButton(ns("server_off"), "Exit",icon("times-circle"), class = "btn btn-danger"),
                 )
          ),
          tags$h2(tags$b("View Map")), br(), hr(),
          
          # column(3,
          #        h4("Legend"),
          #        includeHTML(system.file("ext/include.html", package="viewpoly"))
          # ),
          column(4,
                 checkboxGroupInput(ns("phenotypes"),
                                    label = h4("Phenotypes"),
                                    choices = "This will be updated",
                                    selected = "This will be updated"), br(),
                 
          ),
          column(5,
                 selectInput(inputId = ns("group"), label = p("Linkage group"), choices = 1:15, selected = 1),
                 checkboxInput(ns("op"), label = "Show SNP names", value = TRUE), br(), br(), br(),
                 actionButton(ns("create_server"), "Create local server",icon("refresh"))
          ),
        ), hr(),
        wellPanel(
          sliderInput(ns("range"), "Map range (cM)", min = 0, max = 300,
                      value = c(0, 20), step = 1), 
          uiOutput(ns("interval"))
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
mod_map_view_server <- function(input, output, session, loadMap, loadJBrowse, loadQTL){
  ns <- session$ns
  
  #  Trying to fix server issue
  observeEvent(input$server_off, {
    httpuv::stopAllServers()
  })
  
  
  observe({
    # Dynamic linkage group number
    group_choices <- as.list(1:length(loadMap()$dp))
    names(group_choices) <- 1:length(loadMap()$dp)
    
    updateSelectInput(session, "group",
                      label="Linkage group",
                      choices = group_choices,
                      selected= group_choices[[1]])
    
    # Dynamic QTLs
    pheno_choices <- as.list(unique(loadQTL()$pheno))
    names(pheno_choices) <- unique(loadQTL()$pheno)
    
    updateCheckboxGroupInput(session, "phenotypes",
                             label = "Phenotypes",
                             choices = pheno_choices,
                             selected=unlist(pheno_choices)[1])
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
  
  
  qtl.int <- reactive({
    data <- loadQTL() %>% filter(pheno %in% input$phenotypes & LG == input$group)

    if(dim(data)[1] == 0) stop("No QTLs available in this group")
    
    data <- data[order(data$Pos_lower, data$Pos_upper),]
    command <- paste0(round(data$Pos_lower,0), ":", round(data$Pos_upper, 0))
    seqs <- list()
    for(i in 1:length(command))
      seqs[[i]] <- eval(parse(text = command[i]))
    
    max_updated <- map_summary(left.lim = input$range[1], 
                               right.lim = input$range[2], 
                               ch = input$group, maps = loadMap()$maps, 
                               dp = loadMap()$dp, dq = loadMap()$dq)[[5]]
    
    qtls_pos <- Reduce(union, seqs)
    chr_all <- 0:max_updated
    
    idx.comp <-  chr_all %in% qtls_pos
    int <- chr_all[sequence(rle(idx.comp)$length) == 1]
    
    int <- (int*100)/max_updated
    # add start and end
    ints_all <- unique(c(0,int, 100))
    # add qtls 
    qtls <- unique(sort(data$Pos))
    
    ints_all <- diff(ints_all)
    ints_all[length(ints_all)] <- ints_all[length(ints_all)] -1.5
    
    divs <- vector()
    for(i in 1:length(ints_all)){
      if(idx.comp[1]){ # If 0 is included in some qtl
        if(i %% 2 != 0){
          divs_temp <- paste0("display:inline-block; width: ", ints_all[i] ,"% ; background-color: blue;")
        } else {
          divs_temp <- paste0("display:inline-block; width: ", ints_all[i] ,"% ; background-color: gray;")
        }
      } else {
        if(i %% 2 != 0){
          divs_temp <- paste0("display:inline-block; width: ", ints_all[i] ,"% ; background-color: gray;")
        } else {
          divs_temp <- paste0("display:inline-block; width: ", ints_all[i] ,"% ; background-color: blue;")
        }
      }
      divs <- c(divs, divs_temp)
    }

    if(!is.null(input$phenotypes)){
      divs_lst <- list()
      for(i in 1:length(divs)){
        divs_lst[[i]] <- div(id= paste0("belowslider",i), style= divs[i], p())
      }
      p(divs_lst, "QTLs")
    }
  })
  
  output$interval <- renderUI({ 
    qtl.int()
  })
  
  
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
