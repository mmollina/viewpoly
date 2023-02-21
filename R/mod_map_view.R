#' map_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS useShinyjs
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS tagList  
#' 
#' @noRd 
#'
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
                     div(
                       actionButton(ns("goGenes"), "Go to Genome",icon("arrow-circle-left", verify_fa = FALSE), class = "btn btn-primary"),
                       actionButton(ns("goHidecan"), label = div("Go to HIDECAN", icon("arrow-circle-right", verify_fa = FALSE)), class = "btn btn-primary"))
                 )
          ),
          tags$h2(tags$b("VIEWmap")), br(), hr(),
          
          column(6,
                 column(12,
                        box(
                          background = "light-blue",
                          "Required inputs (*)", br(),
                        )
                 ),
                 column(6,
                        box(width = 12, solidHeader = TRUE,  status="info", title = "Select phenotypes *",
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
                        box(width = 12, solidHeader = TRUE, status="info", title = "Select linkage group *",
                            selectInput(inputId = ns("group"), label = p("Linkage group:"), choices = 1:15, selected = 1),
                            checkboxInput(ns("op"), label = "Show SNP names", value = TRUE)
                        ), br(),
                 )
          ),
        ), hr(),
        wellPanel(
          sliderInput(ns("range"), "Map range (cM)", min = 0, max = 300,
                      value = c(0, 20), step = 1), 
          uiOutput(ns("interval"))
        ),
        box(id= ns("box_profile"), width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = actionLink(inputId = ns("profileID"), label = "QTL profile"),
            column(12,
                   box(
                     width = 5, background = "light-blue",
                     "* QTL analysis files or viewpoly object or example dataset (check `Input data` tab)", 
                   )
            ), 
            column(3,
                   tags$head(tags$style(".butt{background-color:#add8e6; border-color: #add8e6; color: #337ab7;}")),
                   useShinyjs(),
                   downloadButton(ns('bn_download'), "Download", class = "butt")
            ),
            column(3,
                   radioButtons(ns("fformat"), "File type", choices=c("png","tiff","jpeg","pdf", "RData"), selected = "png", inline = T)
            ),                     
            column(2,
                   numericInput(ns("width_profile"), "Width (mm)", value = 180),
            ),
            column(2,
                   numericInput(ns("height_profile"), "Height (mm)", value = 120),
            ),
            column(2,
                   numericInput(ns("dpi_profile"), "DPI", value = 300)
            ), br(),
            column(12,
                   hr(),
                   plotlyOutput(ns("plot_qtl")), 
            )
        ), br(),
        box(id = ns("box_map"), width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = actionLink(inputId = ns("mapID"), label = "Map"),
            column(12,
                   box(
                     width = 5, background = "light-blue",
                     "* Linkage map files or viewpoly object or example dataset (check `Input data` tab)", 
                   )
            ), 
            column(3,
                   downloadButton(ns('bn_download_map'), "Download", class = "butt")
            ),
            column(3,
                   radioButtons(ns("fformat_map"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
            ),                     
            column(2,
                   numericInput(ns("width_map"), "Width (mm)", value = 180),
            ),
            column(2,
                   numericInput(ns("height_map"), "Height (mm)", value = 120),
            ),
            column(2,
                   numericInput(ns("dpi_map"), "DPI", value = 300)
            ), br(),
            column(12,
                   hr(),
                   plotOutput(ns("plot_map"), height = "500px"), br(),
                   includeHTML(system.file(package = "viewpoly", "ext/include.html")), br(), br(),
                   box(id = ns("box_phaplo"),width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = actionLink(inputId = ns("phaploID"), label = "Parents haplotypes table"),
                       DT::dataTableOutput(ns("parents_haplo"))
                   )
            )
        ),
        box(id = ns("box_mapsumm"), width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = actionLink(inputId = ns("mapsummID"), label = "Map summary"),
            column(12,
                   box(
                     width = 5, background = "light-blue",
                     "* Linkage map files or viewpoly object or example dataset (check `Input data` tab)", 
                   )
            ), 
            column(12,
                   DT::dataTableOutput(ns("summary")), br(), hr()
            ),
            column(3,
                   downloadButton(ns('bn_download_summary'), "Download", class = "butt")
            ),
            column(3,
                   radioButtons(ns("fformat_summary"), "File type", choices=c("png","tiff","jpeg","pdf", "RData"), selected = "png", inline = T), br(),
            ),                     
            column(2,
                   numericInput(ns("width_summary"), "Width (mm)", value = 180),
            ),
            column(2,
                   numericInput(ns("height_summary"), "Height (mm)", value = 120),
            ),
            column(2,
                   numericInput(ns("dpi_summary"), "DPI", value = 300)
            ), br(),
            column(12,
                   hr(),
                   plotOutput(ns("map_summary"))
            )
        )
      )
    )
  )
}

#' map_view Server Functions
#'
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom dplyr `%>%`
#' @importFrom shinyjs js
#' @noRd 
mod_map_view_server <- function(input, output, session, 
                                loadMap, loadQTL,
                                parent_session){
  ns <- session$ns
  
  pheno <- LG <- NULL
  
  #Collapse boxes
  observeEvent(input$profileID, {
    js$collapse(ns("box_profile"))
  })
  
  observeEvent(input$mapID, {
    js$collapse(ns("box_map"))
  })
  
  observeEvent(input$mapsummID, {
    js$collapse(ns("box_mapsumm"))
  })
  
  observeEvent(input$phaploID, {
    js$collapse(ns("box_phaplo"))
  })
  
  observeEvent(input$goGenes, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "genes")
  })
  
  observeEvent(input$goHidecan, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "hidecan")
  })
  
  observe({
    # Dynamic linkage group number
    if(!is.null(loadMap())){
      group_choices <- as.list(1:length(loadMap()$d.p1))
      names(group_choices) <- 1:length(loadMap()$d.p1)
    } else if(!is.null(loadQTL())){
      group_choices <- as.list(1:length(unique(loadQTL()$selected_mks$LG)))
      names(group_choices) <- 1:length(unique(loadQTL()$selected_mks$LG))
    } else {
      group_choices <- as.list("Upload map or QTL data in `upload` session.")
      names(group_choices) <-  "Upload map or QTL data in `upload` session."
    }
    
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
    } else {
      updatePickerInput(session, "phenotypes",
                        label = "Phenotype:",
                        choices = "Upload QTL information to update",
                        selected= "Upload QTL information to update")
    }
  })
  
  # Plot QTL bar
  qtl.int <- reactive({
    if(!is.null(loadQTL())){
      data <- loadQTL()$qtl_info %>% filter(pheno %in% input$phenotypes & LG == input$group)
      
      if(dim(data)[1] == 0) return(p(" "))
      
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
      # add qtl
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
      OrRd <- c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000")
      if(length(labs[which(labs == "red")]) < 3){
        qtl.colors <- OrRd[1:7][-c(1:5)][1:length(labs[which(labs == "red")])]
      } else {
        qtl.colors <- OrRd[1:length(labs[which(labs == "red")])]
      }
      
      labs[which(labs == "red")][order(as.numeric(data$Pval), decreasing = T)] <- qtl.colors
      
      divs <- paste0("display:inline-block; width: ", ints_all ,"% ; background-color: ", labs, ";")
      if(!is.null(input$phenotypes)){
        divs_lst <- list()
        for(i in 1:length(divs)){
          divs_lst[[i]] <- div(id= paste0("belowslider",i), style= divs[i], p())
        }
        p(divs_lst, "QTL")
      }
    } else 
      return(p(" "))
  })
  
  output$interval <- renderUI({ 
    qtl.int()
  })
  
  # Plot QTL profile
  output$plot_qtl <- renderPlotly({
    validate(
      need(!is.null(loadQTL()), "Upload the QTL information in upload session to access this feature.")
    )
    idx <- which(unique(loadQTL()$profile$pheno) %in% input$phenotypes)
    pl <- plot_profile(profile = loadQTL()$profile, 
                       qtl_info = loadQTL()$qtl_info, 
                       selected_mks = loadQTL()$selected_mks,
                       pheno.col = idx,
                       lgs.id = as.numeric(input$group),
                       range.min = input$range[1],
                       range.max = input$range[2], by_range=T)
    
    ggplotly(source = "qtl_profile", pl, tooltip=c("Trait", "Position (cM)")) %>% 
      layout(legend = list(orientation = 'h', y = -0.3), 
             modebar = list(
               remove = c("toImage", 
                          "hovercompare", 
                          "hoverCompareCartesian")),
             clickmode ="none",
             dragmode = FALSE)
  })
  
  # Plot map
  output$plot_map <- renderPlot({
    validate(
      need(!is.null(loadMap()$ph.p1), "Upload map information in the upload session to access this feature.")
    )
    maps <- lapply(loadMap()$maps, function(x) {
      y <- x$l.dist
      names(y) <- x$mk.names
      y
    })
    draw_map_shiny(left.lim = input$range[1], 
                   right.lim = input$range[2], 
                   ch = input$group,
                   d.p1 = loadMap()$d.p1,
                   d.p2 = loadMap()$d.p2, 
                   maps = maps, 
                   ph.p1 = loadMap()$ph.p1, 
                   ph.p2 = loadMap()$ph.p2,
                   snp.names = input$op)
    
    max_updated = reactive({
      map_summary(left.lim = input$range[1], right.lim = input$range[2], 
                  ch = input$group, maps = maps, 
                  d.p1 = loadMap()$d.p1, 
                  d.p2 = loadMap()$d.p2)[[5]]
    })
    
    observeEvent(max_updated, {
      updateSliderInput(inputId = "range", max = round(max_updated(),2))
    })
  })
  
  output$parents_haplo  <- DT::renderDataTable(server = FALSE, {
    validate(
      need(!is.null(loadMap()$ph.p1), "Upload map information in the upload session to access this feature.")
    )
    group <- as.numeric(input$group)
    mks<- loadMap()$maps[[group]]
    mks <- mks[order(mks$l.dist),]
    mks.range <- which(mks$l.dist >= input$range[1] &  mks$l.dist <= input$range[2])
    p1 <- loadMap()$ph.p1[[group]][mks.range,]
    #colnames(p1) <- paste0("p1.",1:dim(p1)[2])
    p2 <- loadMap()$ph.p2[[group]][mks.range,]
    #colnames(p2) <- paste0("p2.",1:dim(p2)[2])
    p.haplo <- cbind(p1,p2)
    
    DT::datatable(p.haplo, extensions = 'Buttons', 
                  options = list(
                    dom = 'Bfrtlp',
                    buttons = c('copy', 'csv', 'excel', 'pdf')
                  ),
                  class = "display")
  })
  
  # Map summary
  output$summary  <- DT::renderDataTable(server = FALSE, {
    validate(
      need(!is.null(loadMap()$ph.p1), "Upload map information in the upload session to access this feature.")
    )
    summary <- summary_maps(loadMap())
    DT::datatable(summary, extensions = 'Buttons', 
                  options = list(
                    scrollX = TRUE,
                    dom = 'Bfrtlp',
                    buttons = c('copy', 'csv', 'excel', 'pdf') 
                  ),
                  class = "display")
  })
  
  output$map_summary <- renderPlot({
    validate(
      need(!is.null(loadMap()$ph.p1), "Upload map information in the upload session to access this feature.")
    )
    plot_map_list(loadMap())
  })
  
  ### Downloads
  
  # QTL profile
  fn_downloadname <- reactive({
    seed <- sample(1:1000,1)
    if(input$fformat=="png") filename <- paste0("profile","_",seed,".png")
    if(input$fformat=="tiff") filename <- paste0("profile","_",seed,".tiff")
    if(input$fformat=="jpeg") filename <- paste0("profile","_",seed,".jpg")
    if(input$fformat=="pdf") filename <- paste0("profile","_",seed,".pdf")
    if(input$fformat=="RData") filename <- paste0("profile","_",seed,".RData")
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
    
    if(input$fformat!="RData"){
      ggsave(pl, filename = fn_downloadname(), 
             width = input$width_profile, height = input$height_profile, 
             units = "mm", dpi = input$dpi_profile)    
    } else save(pl, file = fn_downloadname())
  }
  
  observe({
    if (!is.null(loadQTL()) & input$width_profile > 1 & input$height_profile > 1 & input$dpi_profile > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download")
    } else {
      shinyjs::disable("bn_download")
    }
  })
  
  # download handler
  output$bn_download <- downloadHandler(
    filename = fn_downloadname,
    content = function(file) {
      fn_download()
      file.copy(fn_downloadname(), file, overwrite=T)
      file.remove(fn_downloadname())
    }
  )
  
  # Map file
  fn_downloadname_map <- reactive({
    png <- tiff <- jpeg <- pdf <- dev.off <- NULL
    seed <- sample(1:1000,1)
    if(input$fformat_map=="png") filename <- paste0("map","_",seed,".png")
    if(input$fformat_map=="tiff") filename <- paste0("map","_",seed,".tiff")
    if(input$fformat_map=="jpeg") filename <- paste0("map","_",seed,".jpg")
    if(input$fformat_map=="pdf") filename <- paste0("map","_",seed,".pdf")
    return(filename)
  })
  
  # download profile 
  fn_download_map <- function()
  {
    maps <- lapply(loadMap()$maps, function(x) {
      y <- x$l.dist
      names(y) <- x$mk.names
      y
    })
    
    if(input$fformat_map == "png"){
      png(fn_downloadname_map(), width = input$width_map, height = input$height_map, units = "mm", res = input$dpi_map)
    } else  if(input$fformat_map == "tiff"){
      tiff(fn_downloadname_map(), width = input$width_map, height = input$height_map, units = "mm", res = input$dpi_map)
    } else  if(input$fformat_map == "jpeg"){
      jpeg(fn_downloadname_map(), width = input$width_map, height = input$height_map, units = "mm", res = input$dpi_map)
    } else  if(input$fformat_map == "pdf"){
      pdf(fn_downloadname_map(), width = input$width_map, height = input$height_map, units = "mm", res = input$dpi_map)
    } 
    
    draw_map_shiny(left.lim = input$range[1], 
                   right.lim = input$range[2], 
                   ch = input$group,
                   d.p1 = loadMap()$d.p1,
                   d.p2 = loadMap()$d.p2, 
                   maps = maps, 
                   ph.p1 = loadMap()$ph.p1, 
                   ph.p2 = loadMap()$ph.p2,
                   snp.names = input$op)   
    
    dev.off()
  }
  
  observe({
    if (!is.null(loadMap()) & input$width_map > 1 & input$height_map > 1 & input$dpi_map > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download_map")
    } else {
      shinyjs::disable("bn_download_map")
    }
  })
  
  # download handler
  output$bn_download_map <- downloadHandler(
    filename = fn_downloadname_map,
    content = function(file) {
      fn_download_map()
      file.copy(fn_downloadname_map(), file, overwrite=T)
      file.remove(fn_downloadname_map())
    }
  )
  
  # Summary map
  fn_downloadname_summary <- reactive({
    png <- tiff <- jpeg <- pdf <- dev.off <- NULL
    seed <- sample(1:1000,1)
    if(input$fformat_summary=="png") filename <- paste0("map","_",seed,".png")
    if(input$fformat_summary=="tiff") filename <- paste0("map","_",seed,".tiff")
    if(input$fformat_summary=="jpeg") filename <- paste0("map","_",seed,".jpg")
    if(input$fformat_summary=="pdf") filename <- paste0("map","_",seed,".pdf")
    if(input$fformat_summary=="RData") filename <- paste0("map","_",seed,".RData")
    return(filename)
  })
  
  # download profile 
  fn_download_summary <- function()
  {
    png <- tiff <- jpeg <- pdf <- dev.off <- NULL
    
    maps <- lapply(loadMap()$maps, function(x) {
      y <- x$l.dist
      names(y) <- x$mk.names
      y
    })
    
    if(input$fformat_summary != "RData"){
      if(input$fformat_summary == "png"){
        png(fn_downloadname_summary(),  width = input$width_summary, height = input$height_summary, units = "mm", res = input$dpi_summary)
      } else  if(input$fformat_summary == "tiff"){
        tiff(fn_downloadname_summary(),  width = input$width_summary, height = input$height_summary, units = "mm", res = input$dpi_summary)
      } else  if(input$fformat_summary == "jpeg"){
        jpeg(fn_downloadname_summary(),  width = input$width_summary, height = input$height_summary, units = "mm", res = input$dpi_summary)
      } else  if(input$fformat_summary == "pdf"){
        pdf(fn_downloadname_summary(),  width = input$width_summary, height = input$height_summary, units = "mm", res = input$dpi_summary)
      }
      
      plot_map_list(loadMap())   
      
      dev.off()
    } else {
      p <-  plot_map_list(loadMap())   
      save(p, file = fn_downloadname_summary())
    }
  }
  
  observe({
    if (!is.null(loadMap()) & input$width_summary > 1 & input$height_summary > 1 & input$dpi_summary > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download_summary")
    } else {
      shinyjs::disable("bn_download_summary")
    }
  })
  
  # download handler
  output$bn_download_summary <- downloadHandler(
    filename = fn_downloadname_summary,
    content = function(file) {
      fn_download_summary()
      file.copy(fn_downloadname_summary(), file, overwrite=T)
      file.remove(fn_downloadname_summary())
    }
  )
}

## To be copied in the UI
# mod_map_view_ui("map_view_ui_1")

## To be copied in the server
# mod_map_view_server("map_view_ui_1")
