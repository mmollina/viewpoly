#' genes_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
#' 
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_genes_view_ui <- function(id){
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
                     actionButton(ns("server_off"), "Exit",icon("times-circle"), class = "btn btn-danger"),  br(), br(),
                     actionButton(ns("goMap"), "Next",icon("arrow-circle-right"), class = "btn btn-success")
                 )
          ),
          tags$h2(tags$b("View Genome Browser")), br(), hr(),
          
          column(6,
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", title = h4("Select phenotypes"),
                            pickerInput(ns("phenotypes"),
                                        label = h4("Phenotypes:"),
                                        choices = "This will be updated",
                                        selected = "This will be updated",
                                        options = list(
                                          `actions-box` = TRUE, 
                                          size = 10,
                                          `selected-text-format` = "count > 3"
                                        ), 
                                        multiple = TRUE),
                        ), br(),
                 ),
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", title = h4("Select linkage group"),
                            selectInput(inputId = ns("group"), label = p("Linkage group:"), choices = 1:15, selected = 1),
                        ), br(),
                 )
          ),
        ), hr(),
        wellPanel(
          sliderInput(ns("range"), "Map range (cM)", min = 0, max = 300,
                      value = c(0, 20), step = 1), 
          uiOutput(ns("interval"))
        ),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("QTL profile"),
            column(2,
                   downloadBttn(ns('bn_download'), style = "gradient", color = "royal")
            ),
            column(10,
                   radioButtons(ns("fformat"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
            ), br(),
            column(12,
                   hr(),
                   plotlyOutput(ns("plot_qtl"))
            )
        ), br(),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("Linkage Map position (cM) x Physical position (Mp)"),
            column(2,
                   downloadBttn(ns('bn_download_phi'), style = "gradient", color = "royal")
            ),
            column(10,
                   radioButtons(ns("fformat_phi"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
            ), br(),
            column(12,
                   hr(),
                   plotlyOutput(ns("plot_pos"))
            )
        ), br(),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("JBrowseR"),
            actionButton(ns("create_server"), "Open JBrowseR",icon("sync")), br(),
            JBrowseROutput(ns("browserOutput"))
        ),br(),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("Genes table"),
            DT::dataTableOutput(ns("genes_ano"))
        )
      )
    )
  )
}

#' genes_view Server Functions
#'
#' @importFrom JBrowseR serve_data renderJBrowseR assembly track_feature tracks default_session JBrowseR JBrowseROutput 
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom plotly event_data layout
#' @importFrom shinyjs inlineCSS
#' @importFrom httpuv stopAllServers
#'
#' @noRd 
mod_genes_view_server <- function(input, output, session, 
                                  loadMap, loadQTL,
                                  loadJBrowse_fasta, loadJBrowse_gff3, loadJBrowse_vcf, 
                                  loadExample,
                                  parent_session){
  ns <- session$ns
  
  #  Trying to fix server issue
  observeEvent(input$server_off, {
    stopAllServers()
  })
  
  observe({
    # Dynamic linkage group number
    group_choices <- as.list(1:length(loadMap()$d.p1))
    names(group_choices) <- 1:length(loadMap()$d.p1)
    
    updateSelectInput(session, "group",
                      label="Linkage group",
                      choices = group_choices,
                      selected= group_choices[[1]])
    
    # Dynamic QTL
    if(!is.null(loadQTL())){
      pheno_choices <- as.list(unique(loadQTL()$profile$pheno))
      names(pheno_choices) <- unique(loadQTL()$profile$pheno)
      
      updatePickerInput(session, "phenotypes",
                        label = "Select phenotypes",
                        choices = pheno_choices,
                        selected=unlist(pheno_choices)[1])
    }
  })
  
  observeEvent(input$goMap, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "map")
  })
  
  # Plot QTL bar
  qtl.int <- reactive({
    if(!is.null(loadQTL())){
      data <- loadQTL()$qtl_info %>% filter(pheno %in% input$phenotypes & LG == input$group)
      
      if(dim(data)[1] == 0) stop("No QTL available in this group")
      
      data <- data[order(data$Pos_lower, data$Pos_upper),]
      command <- paste0(round(data$Pos_lower,0), ":", round(data$Pos_upper, 0))
      seqs <- list()
      for(i in 1:length(command))
        seqs[[i]] <- eval(parse(text = command[i]))
      
      maps <- lapply(loadMap()$maps, function(x) {
        y <- x$l.dist
        names(y) <- x$mk.names
        y
      })
      
      max_updated <- map_summary(left.lim = input$range[1], 
                                 right.lim = input$range[2], 
                                 ch = input$group, maps = maps, 
                                 d.p1 = loadMap()$d.p1, d.p2 = loadMap()$d.p2)[[5]]
      
      qtls_pos <- Reduce(union, seqs)
      chr_all <- 0:max_updated
      
      idx.comp <-  chr_all %in% qtls_pos
      int <- chr_all[sequence(rle(idx.comp)$length) == 1]
      
      int <- (int*100)/max_updated
      # add start and end
      ints_all <- unique(c(0,int, 100))
      # add qtls 
      qtls <- (unique(sort(data$Pos))*100)/max_updated
      qtls <- sort(c(qtls -0.3, qtls +0.3))
      
      labs <- c(rep("int", length(ints_all)), rep(c("red","#34495E "), length(qtls/2)))
      labs <- labs[order(c(ints_all, qtls))]
      labs[which(labs == "red")-1] <- "#34495E "
      labs[which(labs == "int")] <- "#D5D8DC"
      labs <- labs[-length(labs)]
      
      ints_all <- diff(sort(c(ints_all, qtls)))
      
      # Each interval add small blank space to the scale - need to remove
      reduce <- cumsum(ints_all)[length(cumsum(ints_all))] - 99.7
      ints_all[which(labs != "red")] <- ints_all[which(labs != "red")] - reduce
      
      # Add gradient colors
      if(length(labs[which(labs == "red")]) < 3){
        qtl.colors <- brewer.pal(7, name = "OrRd")[-c(1:5)][1:length(labs[which(labs == "red")])]
      } else {
        qtl.colors <- brewer.pal(length(labs[which(labs == "red")]), name = "OrRd")
      }
      
      labs[which(labs == "red")][order(as.numeric(data$Pval), decreasing = T)] <- qtl.colors
      
      divs <- paste0("display:inline-block; width: ", ints_all ,"% ; background-color: ", labs, ";")
      if(!is.null(input$phenotypes)){
        divs_lst <- list()
        for(i in 1:length(divs)){
          divs_lst[[i]] <- div(id= paste0("belowslider",i), style= divs[i], p())
        }
      }
      p(divs_lst, "QTL")
    } else {
      stop("Upload the QTL information in upload session to access this feature.")
    }
  })
  
  output$interval <- renderUI({ 
    qtl.int()
  })
  
  # Plot QTL profile
  output$plot_qtl <- renderPlotly({
    if(!is.null(loadQTL())){
      idx <- which(unique(loadQTL()$profile$pheno) %in% input$phenotypes)
      pl <- plot_profile(profile = loadQTL()$profile, qtl_info = loadQTL()$qtl_info, selected_mks = loadQTL()$selected_mks,
                         pheno.col = idx,
                         lgs.id = as.numeric(input$group),
                         range.min = input$range[1],
                         range.max = input$range[2], 
                         by_range=T, 
                         software = loadQTL()$software)
      ggplotly(source = "qtl_profile", pl, tooltip=c("Trait","Position (cM)")) %>% layout(legend = list(orientation = 'h', y = -0.3))
    } else 
      stop("Upload the QTL information in upload session to access this feature.")
  })
  
  # Reactive to change page with click
  s <- reactive({
    event_data("plotly_click", source = "qtl_profile")
  })
  
  observeEvent(s(), {
    updateNavbarPage(session = parent_session, inputId = "viewpoly", selected = "qtl")
  }) 
  
  # Open JBrowser server 
  button <- eventReactive(input$create_server, {
    
    if(!is.null(loadJBrowse_fasta())){
      server_dir <- tempdir()
      
      path.fa <- paste0(server_dir, "/", loadJBrowse_fasta()$name[1])
      path.fai <- paste0(server_dir, "/", loadJBrowse_fasta()$name[2])
      path.gzi <- paste0(server_dir, "/", loadJBrowse_fasta()$name[3])
      
      file.rename(loadJBrowse_fasta()$datapath[1], path.fa)
      file.rename(loadJBrowse_fasta()$datapath[2], path.fai)
      file.rename(loadJBrowse_fasta()$datapath[3], path.gzi)
    } 
    
    if(!is.null(loadJBrowse_gff3())){
      path.gff <- paste0(server_dir, "/", loadJBrowse_gff3()$name[1])
      path.tbi <- paste0(server_dir, "/", loadJBrowse_gff3()$name[2])
      file.rename(loadJBrowse_gff3()$datapath[1], path.gff)
      file.rename(loadJBrowse_gff3()$datapath[2], path.tbi)
    }
    
    if(is.null(loadJBrowse_fasta()) & !is.null(loadExample())){
      path.fa <- loadExample()$fasta
      path.gff <- loadExample()$gff3
      # Add other tracks
      # variants_track <- track_variant()
      # alignments_track <- track_alignments()
    } else {
      stop("Upload the genome information in upload session to access this feature.")
    }
    
    if(exists("data_server"))
      data_server$stop_server()
    
    data_server <- serve_data(dirname(path.fa), port = 5000)
    
    list(path.fa = path.fa, path.gff = path.gff, data_server = data_server)
  })
  
  # Link the UI with the browser widget
  output$browserOutput <- renderJBrowseR({
    
    assembly <- assembly(
      paste0("http://127.0.0.1:5000/", basename(button()$path.fa)), 
      bgzip = TRUE
    )
    
    ## create configuration for a JB2 GFF FeatureTrack
    annotations_track <- track_feature(
      paste0("http://127.0.0.1:5000/", basename(button()$path.gff)), 
      assembly
    )
    
    ## create the tracks array to pass to browser
    tracks <- tracks(annotations_track)
    
    ## select default window
    group <- as.numeric(input$group)
    mk.pos <- loadMap()$maps[[group]]
    mks <- mk.pos[order(mk.pos$l.dist),]
    mks.range <- which(mks$l.dist >= input$range[1] &  mks$l.dist <= input$range[2])
    mks.range.1 <- mks$g.dist[mks.range[1]]
    mks.range.2 <- mks$g.dist[mks.range[length(mks.range)]]
    
    if(mks.range.1 > mks.range.2) stop("Inverted region. Check graphic `Genomic position (bp) x Linkage Map position (cM)`")
    
    default_session <- default_session(
      assembly,
      c(annotations_track)
    )
    
    theme <- JBrowseR::theme("#6c81c0", "#22284c")
    JBrowseR(
      "View",
      assembly = assembly,
      tracks = tracks,
      location = paste0(unique(mks$g.chr),":", mks.range.1,"..",mks.range.2), 
      defaultSession = default_session,
      theme = theme
    )
  })
  
  output$genes_ano  <- DT::renderDataTable(server = FALSE, {
    if(!is.null(button()$path.gff)) {
      group <- as.numeric(input$group)
      mks<- loadMap()$maps[[group]]
      mks <- mks[order(mks$l.dist),]
      mks.range <- which(mks$l.dist >= input$range[1] &  mks$l.dist <= input$range[2])
      mks.range.1 <- mks$g.dist[mks.range[1]]
      mks.range.2 <- mks$g.dist[mks.range[length(mks.range)]]
      
      df <- vroom(button()$path.gff, delim = "\t", skip = 3, col_names = F)
      colnames(df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
      df <- df %>% filter(start > mks.range.1 & end < mks.range.2)
      DT::datatable(df, extensions = 'Buttons',
                    options = list(
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    } else 
      stop("Upload annotation file (.gff3) in the upload session to access this feature.")
  })
  
  output$plot_pos <- renderPlotly({
    map.lg <- loadMap()$maps[[as.numeric(input$group)]]
    
    map.lg$high <- map.lg$g.dist
    map.lg$high[round(map.lg$l.dist,5) < input$range[1] | round(map.lg$l.dist,5) > input$range[2]] <- "black"
    map.lg$high[round(map.lg$l.dist,5) >= input$range[1] & round(map.lg$l.dist,5) <= input$range[2]] <- "red"
    
    map.lg$high <- as.factor(map.lg$high)
    p <- ggplot(map.lg, aes(x=l.dist, y = g.dist/1000, colour = high, text = paste("Marker:", mk.names, "\n", 
                                                                                   "Genetic:", round(l.dist,2), "cM \n",
                                                                                   "Genomic:", g.dist/1000, "Mb"))) +
      geom_point() + scale_color_manual(values=c('black','red')) + 
      theme(legend.position = "none") + 
      labs(x = "Linkage map (cM)", y = "Reference genome (Mb)") +
      theme_bw()
    
    max_updated = reactive({
      dist <- loadMap()$maps[[as.numeric(input$group)]]$l.dist
      max.range <- max(dist)
      max.range
    })
    
    observeEvent(max_updated, {
      updateSliderInput(inputId = "range", max = round(max_updated(),2))
    })
    
    ggplotly(p, tooltip="text") %>% layout(showlegend = FALSE)
  })
  
  ## Downloads
  
  # QTL profile
  fn_downloadname <- reactive({
    seed <- sample(1:1000,1)
    if(input$fformat=="png") filename <- paste0("profile","_",seed,".png")
    if(input$fformat=="tiff") filename <- paste0("profile","_",seed,".tif")
    if(input$fformat=="jpeg") filename <- paste0("profile","_",seed,".jpg")
    if(input$fformat=="pdf") filename <- paste0("profile","_",seed,".pdf")
    return(filename)
  })
  
  # download profile 
  fn_download <- function()
  {
    idx <- which(unique(loadQTL()$profile$pheno) %in% input$phenotypes)
    pl <- plot_profile(profile = loadQTL()$profile, qtl_info = loadQTL()$qtl_info, selected_mks = loadQTL()$selected_mks,
                       pheno.col = idx,
                       lgs.id = as.numeric(input$group),
                       range.min = input$range[1],
                       range.max = input$range[2], 
                       by_range=T, 
                       software = loadQTL()$software)
    ggsave(pl, filename = fn_downloadname(), 
           width = 12.7, height = 8, units = "in")    
  }
  
  # download handler
  output$bn_download <- downloadHandler(
    filename = fn_downloadname,
    content = function(file) {
      fn_download()
      file.copy(fn_downloadname(), file, overwrite=T)
    }
  )
  
  # Download cMxMb
  fn_downloadname_phi <- reactive({
    seed <- sample(1:1000,1)
    if(input$fformat_phi=="png") filename <- paste0("linkageXphisical","_",seed,".png")
    if(input$fformat_phi=="tiff") filename <- paste0("linkageXphisical","_",seed,".tif")
    if(input$fformat_phi=="jpeg") filename <- paste0("linkageXphisical","_",seed,".jpg")
    if(input$fformat_phi=="pdf") filename <- paste0("linkageXphisical","_",seed,".pdf")
    return(filename)
  })
  
  # download  
  fn_download_phi <- function()
  {
    map.lg <- loadMap()$maps[[as.numeric(input$group)]]
    
    map.lg$high <- map.lg$g.dist
    map.lg$high[round(map.lg$l.dist,5) < input$range[1] | round(map.lg$l.dist,5) > input$range[2]] <- "black"
    map.lg$high[round(map.lg$l.dist,5) >= input$range[1] & round(map.lg$l.dist,5) <= input$range[2]] <- "red"
    
    map.lg$high <- as.factor(map.lg$high)
    p <- ggplot(map.lg, aes(x=l.dist, y = g.dist/1000, colour = high, text = paste("Marker:", mk.names, "\n", 
                                                                                   "Genetic:", round(l.dist,2), "cM \n",
                                                                                   "Genomic:", g.dist/1000, "Mb"))) +
      geom_point() + scale_color_manual(values=c('black','red')) + 
      labs(x = "Linkage map (cM)", y = "Reference genome (Mb)") +
      theme_bw() + theme(legend.position = "none") 
    
    ggsave(p, filename = fn_downloadname_phi(), 
           width = 12.7, height = 8, units = "in")    
  }
  
  # download handler
  output$bn_download_phi <- downloadHandler(
    filename = fn_downloadname_phi,
    content = function(file) {
      fn_download_phi()
      file.copy(fn_downloadname_phi(), file, overwrite=T)
    }
  )
  
}

## To be copied in the UI
# mod_genes_view_ui("genes_view_ui_1")

## To be copied in the server
# mod_genes_view_server("genes_view_ui_1")
