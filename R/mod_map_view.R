#' map_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
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
                         actionButton(ns("goGenes"), "Go to Genome",icon("arrow-circle-left", verify_fa = FALSE), class = "btn btn-primary"))
                 )
          ),
          tags$h2(tags$b("VIEWmap")), br(), hr(),
          
          column(6,
                 column(12,
                        box(
                          width = 3, background = "light-blue",
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
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "QTL profile",
            column(12,
                   box(
                     width = 5, background = "light-blue",
                     "* QTL analysis files or viewpoly object or example dataset (check `Input data` tab)", 
                   )
            ), 
            column(1,
                   downloadBttn(ns('bn_download'), style = "gradient", color = "royal")
            ),
            column(3,
                   radioButtons(ns("fformat"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
            ),                     
            column(2,
                   numericInput(ns("width_profile"), "Width (mm)", value = 180, width = '40%'),
            ),
            column(2,
                   numericInput(ns("height_profile"), "Height (mm)", value = 120, width = '40%'),
            ),
            column(2,
                   numericInput(ns("dpi_profile"), "DPI", value = 300, width = '30%')
            ), br(),
            column(12,
                   hr(),
                   plotlyOutput(ns("plot_qtl")), 
            )
        ), br(),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = "Map",
            column(12,
                   box(
                     width = 5, background = "light-blue",
                     "* Linkage map files or viewpoly object or example dataset (check `Input data` tab)", 
                   )
            ), 
            column(1,
                   downloadBttn(ns('bn_download_map'), style = "gradient", color = "royal")
            ),
            column(3,
                   radioButtons(ns("fformat_map"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
            ),                     
            column(2,
                   numericInput(ns("width_map"), "Width (mm)", value = 180, width = '40%'),
            ),
            column(2,
                   numericInput(ns("height_map"), "Height (mm)", value = 120, width = '40%'),
            ),
            column(2,
                   numericInput(ns("dpi_map"), "DPI", value = 300, width = '30%')
            ), br(),
            column(12,
                   hr(),
                   plotOutput(ns("plot_map"), height = "500px"), br(),
                   includeHTML("www/include.html"), br(), br(),
                   box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "Parents haplotypes table",
                       DT::dataTableOutput(ns("parents_haplo"))
                   )
            )
        ),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "Map summary",
            column(12,
                   box(
                     width = 5, background = "light-blue",
                     "* Linkage map files or viewpoly object or example dataset (check `Input data` tab)", 
                   )
            ), 
            column(12,
                   DT::dataTableOutput(ns("summary")), br(), hr()
            ),
            column(1,
                   downloadBttn(ns('bn_download_summary'), style = "gradient", color = "royal")
            ),
            column(3,
                   radioButtons(ns("fformat_summary"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T), br(),
            ),                     
            column(2,
                   numericInput(ns("width_summary"), "Width (mm)", value = 180, width = '40%'),
            ),
            column(2,
                   numericInput(ns("height_summary"), "Height (mm)", value = 120, width = '40%'),
            ),
            column(2,
                   numericInput(ns("dpi_summary"), "DPI", value = 300, width = '30%')
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
#'
#' @noRd 
mod_map_view_server <- function(input, output, session, 
                                loadMap, loadQTL,
                                parent_session){
  ns <- session$ns
  
  pheno <- LG <- NULL
  
  observeEvent(input$exit, {
    stopApp()
  })
  
  observeEvent(input$goGenes, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "genes")
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
    if(!is.null(loadQTL())){
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
    } else 
      stop(safeError("Upload the QTL information in upload session to access this feature."))
  })
  
  # Plot map
  output$plot_map <- renderPlot({
    if(!is.null(loadMap()$ph.p1)) {
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
    } else 
      stop(safeError("Upload map information in the upload session to access this feature."))
  })
  
  output$parents_haplo  <- DT::renderDataTable(server = FALSE, {
    if(!is.null(loadMap()$ph.p1)) {
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
    } else 
      stop(safeError("Upload map information in the upload session to access this feature."))
  })
  
  # Map summary
  output$summary  <- DT::renderDataTable(server = FALSE, {
    if(!is.null(loadMap()$ph.p1)) {
      summary <- summary_maps(loadMap())
      
      DT::datatable(summary, extensions = 'Buttons', 
                    options = list(
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf') 
                    ),
                    class = "display")
    } else 
      stop(safeError("Upload map information in the upload session to access this feature."))
  })
  
  output$map_summary <- renderPlot({
    if(!is.null(loadMap()$ph.p1)) {
      plot_map_list(loadMap())
    } else 
      stop(safeError("Upload map information in the upload session to access this feature."))
  })
  
  ### Downloads
  
  # QTL profile
  fn_downloadname <- reactive({
    seed <- sample(1:1000,1)
    if(input$fformat=="png") filename <- paste0("profile","_",seed,".png")
    if(input$fformat=="tiff") filename <- paste0("profile","_",seed,".tiff")
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
           width = input$width_profile, height = input$height_profile, 
           units = "mm", dpi = input$dpi_profile)    
  }
  
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
  }
  
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
