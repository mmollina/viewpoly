#' qtl_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @import shinydashboard
#' @import shinyWidgets
#' 
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_qtl_view_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      verticalLayout(
        fluidRow(
          column(width = 12,
                 div(style = "position:absolute;right:1em;", 
                     actionButton(ns("exit"), "Exit",icon("times-circle"), class = "btn btn-danger"), br(), br(),
                     actionButton(ns("goGenes"), "Next",icon("arrow-circle-right"), class = "btn btn-success")
                 )
          ),
          tags$h2(tags$b("View QTL")), br(), hr(),
          column(6,
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", title = h4("Select linkage group/s"),
                            pickerInput(ns("group"),
                                        label = h6("Linkage group/s:"),
                                        choices = "This will be updated",
                                        selected = "This will be updated",
                                        options = list(
                                          `actions-box` = TRUE, 
                                          size = 10,
                                          `selected-text-format` = "count > 3"
                                        ), 
                                        multiple = TRUE)
                        )
                 ),
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", title = h4("Select phenotype/s"),
                            pickerInput(ns("phenotypes"),
                                        label = h6("Phenotype/s:"),
                                        choices = "This will be updated",
                                        selected = "This will be updated",
                                        options = list(
                                          `actions-box` = TRUE, 
                                          size = 10,
                                          `selected-text-format` = "count > 3"
                                        ), 
                                        multiple = TRUE)
                        )
                 )
          ),
          column(12,
                 box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("QTL profile"),
                     column(2,
                            downloadBttn(ns('bn_download'), style = "gradient", color = "royal")
                     ),
                     column(10,
                            radioButtons(ns("fformat"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
                     ), br(), 
                     column(12,
                            hr(),
                            plotOutput(ns("plot_qtl"), 
                                       click=ns("plot_click"), brush = ns("plot_brush"))
                     ),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Effects"),
                         div(style = "position:absolute;right:3em;",
                             radioButtons(ns("effects_design"), "Design", 
                                          choices = c("Additive (bar)" = "bar", "Additive (circle)" = "circle", "Alleles combination" = "digenic"), 
                                          selected = "bar", inline= T)
                         ), br(), br(), 
                         column(2,
                                downloadBttn(ns('bn_download_effects'), style = "gradient", color = "royal")
                         ),
                         column(10,
                                radioButtons(ns("fformat_effects"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
                         ), br(),
                         column(12,
                                hr(),
                                uiOutput(ns("plot.ui"))
                         )
                     ), br(),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Progeny haplotypes"),
                         column(12,
                                pickerInput(ns("haplo"),
                                            label = h6("Select haplotypes"),
                                            choices = "Select QTL in the profile graphic to update",
                                            selected = "Select QTL in the profile graphic to update",
                                            options = pickerOptions(
                                              size = 15,
                                              `selected-text-format` = "count > 3",
                                              `live-search`=TRUE,
                                              actionsBox = TRUE,
                                              dropupAuto = FALSE,
                                              dropdownAlignRight = TRUE
                                            ), 
                                            multiple = TRUE), br(),
                                actionBttn(ns("haplo_submit"), style = "jelly", color = "royal",  size = "sm", label = "submit selected haplotypes", icon = icon("share-square")), 
                                br(), hr()),
                         column(2,
                                downloadBttn(ns('bn_download_haplo'), style = "gradient", color = "royal")
                         ),
                         column(10,
                                radioButtons(ns("fformat_haplo"), "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", inline = T)
                         ), br(),
                         column(12,
                                hr(),
                                uiOutput(ns("plot_haplo.ui"))
                         )
                     ),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Breeding values"),
                         DT::dataTableOutput(ns("breeding_values"))
                     ),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("QTL summary"),
                         DT::dataTableOutput(ns("info"))
                     )
                 )
          )
        )
      )
    )
  )
}

#' qtl_view Server Functions
#'
#' @importFrom ggpubr ggarrange
#' @import shinydashboard
#' 
#' @noRd 
mod_qtl_view_server <- function(input, output, session, 
                                loadMap, loadQTL,
                                parent_session){
  ns <- session$ns
  
  observeEvent(input$exit, {
    stopApp()
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
    updatePickerInput(session, "group",
                      label="Linkage group/s:",
                      choices = group_choices,
                      selected= group_choices[[1]])
    
    
    # Dynamic QTL
    if(!is.null(loadQTL())){
      pheno_choices <- as.list(unique(loadQTL()$profile$pheno))
      names(pheno_choices) <- unique(loadQTL()$profile$pheno)
      
      updatePickerInput(session, "phenotypes",
                        label = "Phenotype/s:",
                        choices = pheno_choices,
                        selected=unlist(pheno_choices)[1])
    } else {
      updatePickerInput(session, "phenotypes",
                        label = "Phenotype/s:",
                        choices = "Upload QTL information to update",
                        selected= "Upload QTL information to update")
    }
  })
  
  observeEvent(input$goGenes, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "genes")
  })
  
  qtl.data <- reactive({
    if(!is.null(loadQTL())){
      idx <- which(unique(loadQTL()$profile$pheno) %in% input$phenotypes)
      pl <- plot_profile(profile = loadQTL()$profile, 
                         qtl_info = loadQTL()$qtl_info, 
                         selected_mks = loadQTL()$selected_mks,
                         pheno.col = idx,
                         lgs.id = as.numeric(input$group), 
                         by_range=F, plot = F)
    } else
      stop("Upload the QTL information in upload session to access this feature.")
  })
  
  output$plot_qtl <- renderPlot({
    only_plot_profile(pl.in = qtl.data())
  })
  
  effects.data <- reactive({
    if(!is.null(loadQTL())){
      if(!is.null(input$plot_brush)){
        df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
      } else {
        stop("Select a point or region on QTL profile graphic.") 
      }
      data <- data_effects(qtl_info = loadQTL()$qtl_info, 
                           effects = loadQTL()$effects,
                           pheno.col = as.character(df$Trait), 
                           lgs = df$LG, 
                           position = df$`Position (cM)`,
                           groups = as.numeric(input$group),
                           software = loadQTL()$software,
                           design = input$effects_design)
    } else 
      stop("Upload the QTL information in upload session to access this feature.")
  })
  
  output$effects <- renderPlot({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("building graphic..."))
      plot_effects(effects.data(), software = loadQTL()$software, design = input$effects_design)
    })
  })
  
  plotHeight <- reactive({
    if(!is.null(loadQTL())){
      if(!is.null(input$plot_brush)){
        dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
      } else {
        stop("Select a point or region on QTL profile graphic.")
      }
      counts <- nrow(dframe)
      counts <- ceiling(counts/4)
      if(counts == 0) counts <- 1
      if(loadQTL()$software == "polyqtlR") {
        size <- counts*650 
      } else if(input$effects_design == "bar" | input$effects_design == "digenic"){ 
        size <- counts*350
      } else if(input$effects_design == "circle"){
        counts <- length(unique(dframe$LG))
        counts <- ceiling(counts/2)
        if(counts == 0) counts <- 1
        size <- counts*500
      }
      size
    } else 
      stop("Upload the QTL information in upload session to access this feature.")
  })
  
  output$plot.ui <- renderUI({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("building graphic..."))
      plotOutput(ns("effects"), height =  plotHeight())
    })
  })
  
  observeEvent(output$plot.ui,{
    if(!is.null(loadQTL())){
      if(loadQTL()$software == "polyqtlR" | loadQTL()$software == "diaQTL") {
        dframe <- NULL
        updatePickerInput(session, "haplo",
                          label = "Select haplotypes",
                          choices = paste0("Feature not implemented for software: ", loadQTL()$software),
                          selected= paste0("Feature not implemented for software: ", loadQTL()$software))
      } else if(!is.null(input$plot_brush)){
        dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
      } else {
        dframe <- NULL
        updatePickerInput(session, "haplo",
                          label = "Select haplotypes",
                          choices = "Select QTL in the profile graphic to update",
                          selected= "Select QTL in the profile graphic to update")
      }
    } else {
      dframe <- NULL
      updatePickerInput(session, "haplo",
                        label = "Select haplotypes",
                        choices = "Upload the QTL information in upload session to access this feature.",
                        selected= "Upload the QTL information in upload session to access this feature.")
    }
    
    if(!is.null(dframe)){
      if(input$effects_design == "digenic" | input$effects_design == "circle") {
        updatePickerInput(session, "haplo",
                          label = "Select haplotypes",
                          choices = "Select `bar` design to access this feature.",
                          selected= "Select `bar` design to access this feature.")
      } else {
        haplo_choices <- paste0("Trait:", dframe$Trait, "_LG:", dframe$LG, "_Pos:", dframe$`Position (cM)`)
        alleles <- effects.data()[[1]]$data$Alleles[!grepl("_",effects.data()[[1]]$data$Alleles)]
        alleles <- rep(alleles, length(haplo_choices))
        haplo_choices <- rep(haplo_choices, each = length(alleles)/length(haplo_choices))
        haplo_choices <- paste0(haplo_choices, "_homolog:", alleles)
        haplo_choices <- as.list(haplo_choices)
        names(haplo_choices) <- unlist(haplo_choices)
        updatePickerInput(session, "haplo",
                          label = "Select haplotypes",
                          choices = haplo_choices,
                          selected= haplo_choices[1:3])
      }
    }
  })
  
  haplo_data <- eventReactive(input$haplo_submit, {
    if(all(input$haplo == paste0("Feature not implemented for software: ", loadQTL()$software))) stop(paste0("Feature not implemented for software: ", loadQTL()$software))
    if(all(input$haplo == "Select QTL in the profile graphic to update")) stop("Select QTL in the profile graphic to update")
    if(all(input$haplo == "Select `bar` design to access this feature.")) stop("Select `bar` design to access this feature.")
    p <- select_haplo(input$haplo, loadQTL()$probs, loadQTL()$selected_mks, effects.data())
    counts <- ceiling(length(p)/3)
    if(counts == 0) counts <- 1
    size <- counts*450
    list(p, size)
  })
  
  output$haplotypes <- renderPlot({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.3, detail = paste("building graphic..."))
      nrow.lst <- ceiling(length(haplo_data()[[1]])/3)
      if(nrow.lst == 0) nrow.lst <- 1
      p.all <- ggarrange(plotlist = haplo_data()[[1]], ncol = 3, nrow = nrow.lst, common.legend = TRUE)
    })
    p.all
  })
  
  output$plot_haplo.ui <- renderUI({
    plotOutput(ns("haplotypes"), height = haplo_data()[[2]])
  })
  
  output$info <- DT::renderDataTable(server = FALSE, {
    if(!is.null(loadQTL())){
      if(!is.null(input$plot_brush)){
        dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
      } else {
        stop("Select a point or region on graphic.")
      }
      dframe <- dframe[,-c(dim(dframe)[2]-1,dim(dframe)[2])]
      if(loadQTL()$software == "QTLpoly"){
        colnames(dframe)[c(2,4,5,6,7)] <- c("Linkage group", "Lower interval (cM)", "Upper interval (cM)", "p-value", "h2")
      } else if(loadQTL()$software == "diaQTL") {
        colnames(dframe)[c(2,4,5,6)] <- c("Linkage group", "Lower interval (cM)", "Upper interval (cM)", "LL")
      } else if(loadQTL()$software == "polyqtlR"){
        dframe <- dframe[,-c(4,5)]
        colnames(dframe)[c(2,4)] <- c("Linkage group", "Threshold")
      }
      DT::datatable(dframe, extensions = 'Buttons',
                    options = list(
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    } else 
      stop("Upload the QTL information in upload session to access this feature.")
  })
  
  # Breeding values
  output$breeding_values <- DT::renderDataTable(server = FALSE, {
    if(!is.null(loadQTL())){
      if(loadQTL()$software == "QTLpoly"){
        if(!is.null(input$plot_brush)){
          dframe <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
        } else {
          stop("Select a point or region on graphic.")
        }
        
        pos <- split(dframe$`Position (cM)`, dframe$Trait)
        dt <- breeding_values(loadQTL()$qtl_info, loadQTL()$probs, 
                              loadQTL()$selected_mks, loadQTL()$blups, 
                              loadQTL()$beta.hat, pos)
        rownames(dt) <- NULL
        DT::datatable(dt, extensions = 'Buttons',
                      options = list(
                        dom = 'Bfrtlp',
                        buttons = c('copy', 'csv', 'excel', 'pdf')
                      ),
                      class = "display")
      } else stop(paste("Feature not implemented for software:",loadQTL()$software))
    } else 
      stop("Upload the QTL information in upload session to access this feature.")
  })
  
  # Download profile
  # create filename
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
    p <- only_plot_profile(pl.in = qtl.data())
    ggsave(p, filename = fn_downloadname(), 
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
  
  # Download effects
  # create filename
  fn_downloadname_effects <- reactive({
    
    seed <- sample(1:1000,1)
    if(input$fformat_effects=="png") filename <- paste0("effects","_",seed,".png")
    if(input$fformat_effects=="tiff") filename <- paste0("effects","_",seed,".tif")
    if(input$fformat_effects=="jpeg") filename <- paste0("effects","_",seed,".jpg")
    if(input$fformat_effects=="pdf") filename <- paste0("effects","_",seed,".pdf")
    return(filename)
  })
  
  # download 
  fn_download_effects <- function()
  {
    if(!is.null(input$plot_brush)){
      df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    } else {
      stop("Select a point or region on QTL profile graphic.") 
    }
    data <- data_effects(qtl_info = loadQTL()$qtl_info, 
                         effects = loadQTL()$effects,
                         pheno.col = as.character(df$Trait), 
                         lgs = df$LG, 
                         position = df$`Position (cM)`,
                         groups = as.numeric(input$group),
                         software = loadQTL()$software,
                         design = input$effects_design)
    
    plots <- plot_effects(data, software = loadQTL()$software, design = input$effects_design)
    
    ggsave(plots, filename = fn_downloadname_effects(), height = plotHeight()/3, width = plotHeight(),units = "mm", bg = "white")    
  }
  
  # download handler
  output$bn_download_effects <- downloadHandler(
    filename = fn_downloadname_effects,
    content = function(file) {
      fn_download_effects()
      file.copy(fn_downloadname_effects(), file, overwrite=T)
    }
  )
  
  # Download haplotypes
  # create filename
  fn_downloadname_haplo <- reactive({
    
    seed <- sample(1:1000,1)
    if(input$fformat_haplo=="png") filename <- paste0("haplotypes","_",seed,".png")
    if(input$fformat_haplo=="tiff") filename <- paste0("haplotypes","_",seed,".tif")
    if(input$fformat_haplo=="jpeg") filename <- paste0("haplotypes","_",seed,".jpg")
    if(input$fformat_haplo=="pdf") filename <- paste0("haplotypes","_",seed,".pdf")
    return(filename)
  })
  
  # download 
  fn_download_haplo <- function()
  {
    p <- select_haplo(input$haplo, loadQTL()$probs, loadQTL()$selected_mks, effects.data())
    plots <- ggarrange(plotlist = p, ncol = 3, common.legend = TRUE)
    
    ggsave(plots, filename = fn_downloadname_haplo(), height = plotHeight()/2, width = plotHeight(),units = "mm", bg = "white")    
  }
  
  # download handler
  output$bn_download_haplo <- downloadHandler(
    filename = fn_downloadname_haplo,
    content = function(file) {
      fn_download_haplo()
      file.copy(fn_downloadname_haplo(), file, overwrite=T)
    }
  )
}

## To be copied in the UI
# mod_qtl_view_ui("qtl_view_ui_1")

## To be copied in the server
# mod_qtl_view_server("qtl_view_ui_1")
