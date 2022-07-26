#' qtl_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @import shinydashboard
#' @import shinyWidgets
#' @importFrom shinyjs inlineCSS useShinyjs
#' 
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_qtl_view_ui <- function(id){
  ns <- NS(id)
  tagList(
    useShinyjs(),
    extendShinyjs(text = jscode, functions = "collapse"),
    fluidPage(
      verticalLayout(
        fluidRow(
          column(width = 12,
                 div(style = "position:absolute;right:1em;", 
                     div(
                       actionButton(ns("goUploads"), "Go to Input data", icon("arrow-circle-left", verify_fa = FALSE), class = "btn btn-primary"),
                       actionButton(ns("goGenes"), label = div("Go to Genome", icon("arrow-circle-right", verify_fa = FALSE)), class = "btn btn-primary"))
                 )
          ),
          tags$h2(tags$b("VIEWqtl")), br(), hr(),
          column(6,
                 column(12,
                        box(
                          background = "light-blue",
                          "Required inputs (*)", br(),
                        )
                 ),
                 column(6,
                        box(width = 12, solidHeader = TRUE, status="info", title = "Select linkage group/s *",
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
                        box(width = 12, solidHeader = TRUE, status="info", title = "Select phenotype/s *",
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
                 box(id = ns("box_profile"), width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = actionLink(inputId = ns("profileID"), label = "QTL profile"),
                     column(12,
                            box(
                              background = "light-blue",
                              "* QTL analysis files or viewpoly object or example dataset (check `Input data` tab)"
                            )
                     ), 
                     column(12,
                            column(3,
                                   useShinyjs(),
                                   tags$head(tags$style(".butt{background-color:#add8e6; border-color: #add8e6; color: #337ab7;}")),
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
                            )), br(), 
                     column(12,
                            hr(),
                            plotOutput(ns("plot_qtl"), 
                                       click=ns("plot_click"), brush = ns("plot_brush"))
                     ),
                     box(id = ns("box_effects"), width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = actionLink(inputId = ns("effectsID"), label = "Effects"),
                         column(12,
                                box(
                                  background = "light-blue",
                                  "* QTL analysis files or viewpoly object or example dataset (check `Input data` tab)", br(),
                                  "* Selection of QTL/s (triangle/s at the bottom of QTL profile graphic)"
                                )
                         ), 
                         div(style = "position:absolute;right:3em;",
                             radioButtons(ns("effects_design"), "Design", 
                                          choices = c("Additive (bar)" = "bar", "Additive (circle)" = "circle", "Alleles combination" = "digenic"), 
                                          selected = "bar")
                         ), br(), br(), 
                         column(3,
                                downloadButton(ns('bn_download_effects'), "Download", class = "butt")
                         ),
                         
                         column(3,
                                numericInput(ns("width_effects"), "Width (mm)", value = 180),
                         ),
                         column(3,
                                numericInput(ns("height_effects"), "Height (mm)", value = 120),
                         ),
                         column(3,
                                numericInput(ns("dpi_effects"), "DPI", value = 300)
                         ), br(), 
                         column(12, 
                                column(6,
                                       radioButtons(ns("fformat_effects"), "File type", choices=c("png","tiff","jpeg","pdf", "RData"), selected = "png", inline = T)
                                ),
                                column(6,
                                       textInput(ns("parents_name"), "Parents name", value = "Example: P1, P2")
                                ),
                         ),
                         column(12,
                                hr(),
                                uiOutput(ns("plot.ui"))
                         )
                     ), br(),
                     box(id = ns("box_haplo"),width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = actionLink(inputId = ns("haploID"), label = "Progeny haplotypes"),
                         column(12,
                                box(
                                  background = "light-blue",
                                  "* QTLpoly analysis files or viewpoly object or example dataset (check `Input data` tab)", br(),
                                  "* Selection of QTL/s (triangle/s at the bottom of QTL profile graphic)"
                                )
                         ), 
                         column(12,
                                actionBttn(ns("haplo_update"), style = "jelly", color = "royal",  size = "sm", label = "update available haplotypes", icon = icon("refresh", verify_fa = FALSE)), 
                                br(), br(),
                                pickerInput(ns("haplo"),
                                            label = h6("Select haplotypes*"),
                                            choices = "Click on `update available haplotype` to update",
                                            selected = "Click on `update available haplotype` to update",
                                            options = pickerOptions(
                                              size = 15,
                                              `selected-text-format` = "count > 3",
                                              `live-search`=TRUE,
                                              actionsBox = TRUE,
                                              dropupAuto = FALSE,
                                              dropdownAlignRight = TRUE
                                            ), 
                                            multiple = TRUE), br(),
                                pickerInput(ns("haplo_exclude"),
                                            label = h6("Exclude haplotypes (optional)"),
                                            choices = "Click on `update available haplotype` to update",
                                            selected = "Click on `update available haplotype` to update",
                                            options = pickerOptions(
                                              size = 15,
                                              `selected-text-format` = "count > 3",
                                              `live-search`=TRUE,
                                              actionsBox = TRUE,
                                              dropupAuto = FALSE,
                                              dropdownAlignRight = TRUE
                                            ), 
                                            multiple = TRUE), br(),
                                actionBttn(ns("haplo_submit"), style = "jelly", color = "royal",  size = "sm", label = "submit selected haplotypes*", icon = icon("share-square", verify_fa = FALSE)), 
                                br(), hr()),
                         column(3,
                                downloadButton(ns('bn_download_haplo'), "Download", class = "butt")
                         ),
                         column(3,
                                radioButtons(ns("fformat_haplo"), "File type", choices=c("png","tiff","jpeg","pdf", "RData"), selected = "png", inline = T)
                         ),                     
                         column(2,
                                numericInput(ns("width_haplo"), "Width (mm)", value = 180),
                         ),
                         column(2,
                                numericInput(ns("height_haplo"), "Height (mm)", value = 120),
                         ),
                         column(2,
                                numericInput(ns("dpi_haplo"), "DPI", value = 300)
                         ), br(), 
                         column(12,
                                hr(),
                                htmlOutput(ns("ind_names")), hr(),
                                uiOutput(ns("plot_haplo.ui"))
                         )
                     ),
                     box(id = ns("box_bree"), width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = actionLink(inputId = ns("breeID"), label = "Breeding values"),
                         column(12,
                                box(
                                  background = "light-blue",
                                  "* QTLpoly analysis files or viewpoly object or example dataset (check `Input data` tab)", br(),
                                  "* Selection of QTL/s (triangle/s at the bottom of QTL profile graphic)"
                                )
                         ), 
                         column(12,
                                DT::dataTableOutput(ns("breeding_values"))
                         )
                     ), br(), br(),
                     box(id = ns("box_summary"),width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = actionLink(inputId = ns("summaryID"), label = "QTL summary"),
                         column(12,
                                box(
                                  background = "light-blue",
                                  "* QTL analysis files or viewpoly object or example dataset (check `Input data` tab)", br(),
                                  "* Selection of QTL/s (triangle/s at the bottom of QTL profile graphic)"
                                )
                         ), 
                         column(12,
                                DT::dataTableOutput(ns("info"))
                         )
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
#' @importFrom shinyjs js
#' 
#' @noRd 
mod_qtl_view_server <- function(input, output, session, 
                                loadMap, loadQTL,
                                parent_session){
  ns <- session$ns
  
  #Collapse boxes
  observeEvent(input$profileID, {
    js$collapse(ns("box_profile"))
  })
  
  observeEvent(input$effectsID, {
    js$collapse(ns("box_effects"))
  })
  
  observeEvent(input$haploID, {
    js$collapse(ns("box_haplo"))
  })
  
  observeEvent(input$breeID, {
    js$collapse(ns("box_bree"))
  })
  
  observeEvent(input$summaryID, {
    js$collapse(ns("box_summary"))
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
    
    if(length(group_choices) < 5) the_choice <- group_choices[[1]] else the_choice <- group_choices[[5]]
    
    updatePickerInput(session, "group",
                      label="Linkage group/s:",
                      choices = group_choices,
                      selected= the_choice)
    
    
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
  
  observeEvent(input$goUploads, {
    updateTabsetPanel(session = parent_session, inputId = "viewpoly",
                      selected = "upload")
  })
  
  qtl.data <- reactive({
    validate(
      need(length(input$phenotypes) != 0, "Select at least one phenotype"),
      need(length(input$group) != 0, "Select at least one linkage group"),
      need(!is.null(loadQTL()), "Upload the QTL information in upload session to access this feature.")
    )
    idx <- which(unique(loadQTL()$profile$pheno) %in% input$phenotypes)
    
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.3, detail = paste("building graphic..."))
      pl <- plot_profile(profile = loadQTL()$profile, 
                         qtl_info = loadQTL()$qtl_info, 
                         selected_mks = loadQTL()$selected_mks,
                         pheno.col = idx,
                         lgs.id = as.numeric(input$group), 
                         by_range=F, plot = F)
    })
  })
  
  output$plot_qtl <- renderPlot({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.3, detail = paste("building graphic..."))
      only_plot_profile(pl.in = qtl.data())
    })
  })
  
  effects.data <- reactive({
    validate(
      need(!is.null(loadQTL()), "Upload the QTL information in upload session to access this feature."),
      need(!is.null(input$plot_brush), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
    df <- try(brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat"))
    validate(
      need(dim(df)[1] > 0, "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
    
    if(!grepl("Example" ,input$parents_name)) {
      parents <- unlist(strsplit(input$parents_name, ","))
    } else parents <- NULL
    
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("Getting data..."))
      data <- data_effects(qtl_info = loadQTL()$qtl_info, 
                           effects = loadQTL()$effects,
                           pheno.col = as.character(df$Trait), 
                           lgs = df$LG, 
                           position = df$`Position (cM)`,
                           groups = as.numeric(input$group),
                           software = loadQTL()$software,
                           design = input$effects_design,
                           parents = parents)
    })
  })
  
  output$effects <- renderPlot({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("building graphic..."))
      plot_effects(effects.data(), software = loadQTL()$software, design = input$effects_design)
    })
  })
  
  plotHeight <- reactive({
    
    validate(
      need(!is.null(loadQTL()), "Upload the QTL information in upload session to access this feature."),
      need(!is.null(input$plot_brush), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
    dframe <- try(brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat"))
    validate(
      need(!inherits(dframe, "try-error"), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
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
  })
  
  output$plot.ui <- renderUI({
    withProgress(message = 'Working:', value = 0, {
      incProgress(0.5, detail = paste("building graphic..."))
      plotOutput(ns("effects"), height =  plotHeight())
    })
  })
  
  observeEvent(input$haplo_update,{
    if(!is.null(loadQTL())){
      if(loadQTL()$software == "polyqtlR" | loadQTL()$software == "diaQTL") {
        dframe <- NULL
        updatePickerInput(session, "haplo",
                          label = "Select haplotypes",
                          choices = paste0("Feature not implemented for software: ", loadQTL()$software),
                          selected= paste0("Feature not implemented for software: ", loadQTL()$software))
        
        updatePickerInput(session, "haplo_exclude",
                          label = "Exclude haplotypes (optional)",
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
        
        updatePickerInput(session, "haplo_exclude",
                          label = "Exclude haplotypes (optional)",
                          choices = "Select QTL in the profile graphic to update",
                          selected= "Select QTL in the profile graphic to update")
      }
    } else {
      dframe <- NULL
      updatePickerInput(session, "haplo",
                        label = "Select haplotypes",
                        choices = "Upload the QTL information in upload session to access this feature.",
                        selected= "Upload the QTL information in upload session to access this feature.")
      
      updatePickerInput(session, "haplo_exclude",
                        label = "Exclude haplotypes (optional)",
                        choices = "Upload the QTL information in upload session to access this feature.",
                        selected= "Upload the QTL information in upload session to access this feature.")
    }
    if(!is.null(dframe)){
      if(input$effects_design == "digenic" | input$effects_design == "circle") {
        updatePickerInput(session, "haplo",
                          label = "Select haplotypes",
                          choices = "Select `bar` design to access this feature.",
                          selected= "Select `bar` design to access this feature.")
        
        updatePickerInput(session, "haplo_exclude",
                          label = "Exclude haplotypes (optional)",
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
        
        updatePickerInput(session, "haplo_exclude",
                          label = "Exclude haplotypes (optional)",
                          choices = haplo_choices,
                          selected= NULL)
      }
    }
  })
  
  haplo_data <- eventReactive(input$haplo_submit, {
    validate(
      need(all(input$haplo != paste0("Feature not implemented for software: ", loadQTL()$software)), paste0("Feature not implemented for software: ", loadQTL()$software)),
      need(all(input$haplo != "Click on `update available haplotype` to update"), "Click on `update available haplotype` to update"),
      need(all(input$haplo != "Select QTL in the profile graphic to update"), "Select QTL in the profile graphic to update"),
      need(all(input$haplo != "Select `bar` design to access this feature."), "Select `bar` design to access this feature.")
    )

    list.p <- select_haplo(input$haplo, loadQTL()$probs, loadQTL()$selected_mks, effects.data(), exclude.haplo = input$haplo_exclude)
    p <- list.p[[1]]
    inds <- list.p[[2]]
    counts <- ceiling(length(p)/3)
    if(counts == 0) counts <- 1
    size <- counts*450
    list(p, size, inds)
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
  
  output$ind_names <- renderUI({
    x <- paste0("<strong>Progeny individuals with selected haplotypes  </strong>: ", paste(haplo_data()[[3]], collapse = ", "))
    HTML(x)
  })
  
  output$info <- DT::renderDataTable(server = FALSE, {
    validate(
      need(!is.null(loadQTL()), "Upload the QTL information in upload session to access this feature."),
      need(!is.null(input$plot_brush), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
    dframe <- try(brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat"))
    validate(
      need(!inherits(dframe, "try-error"), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
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
  })
  
  # Breeding values
  output$breeding_values <- DT::renderDataTable(server = FALSE, {
    validate(
      need(!is.null(loadQTL()), "Upload the QTL information in upload session to access this feature."),
      need(loadQTL()$software == "QTLpoly", paste("Feature not implemented for software:",loadQTL()$software)),
      need(!is.null(input$plot_brush), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
    dframe <- try(brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat"))
    validate(
      need(!inherits(dframe, "try-error"), "Select at least one triangle on the bottom of the QTL profile graphic. The triangles refer to QTL peaks detected. You can click and brush your cursor to select more than one.")
    )
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
  })
  
  # Download profile
  # create filename
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
    p <- only_plot_profile(pl.in = qtl.data())
    
    if(input$fformat!="RData"){
      ggsave(p, filename = fn_downloadname(), 
             width = input$width_profile, height = input$height_profile, units = "mm", dpi = input$dpi_profile)    
    } else save(p, file = fn_downloadname())
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
  
  # Download effects
  # create filename
  fn_downloadname_effects <- reactive({
    
    seed <- sample(1:1000,1)
    if(input$fformat_effects=="png") filename <- paste0("effects","_",seed,".png")
    if(input$fformat_effects=="tiff") filename <- paste0("effects","_",seed,".tiff")
    if(input$fformat_effects=="jpeg") filename <- paste0("effects","_",seed,".jpg")
    if(input$fformat_effects=="pdf") filename <- paste0("effects","_",seed,".pdf")
    if(input$fformat_effects=="RData") filename <- paste0("effects","_",seed,".RData")
    return(filename)
  })
  
  # download 
  fn_download_effects <- function()
  {
    validate(
      need(!is.null(input$plot_brush), "Select a point or region on QTL profile graphic.")
    )
    
    df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    
    if(!grepl("Example" ,input$parents_name)) {
      cat("here2")
      parents <- unlist(strsplit(input$parents_name, ","))
    } else parents <- NULL
    
    data <- data_effects(qtl_info = loadQTL()$qtl_info, 
                         effects = loadQTL()$effects,
                         pheno.col = as.character(df$Trait), 
                         lgs = df$LG, 
                         parents = parents,
                         position = df$`Position (cM)`,
                         groups = as.numeric(input$group),
                         software = loadQTL()$software,
                         design = input$effects_design)
    
    plots <- plot_effects(data, software = loadQTL()$software, design = input$effects_design)
    
    if(input$fformat_effects!="RData"){
      ggsave(plots, filename = fn_downloadname_effects(), height = input$height_effects, 
             width = input$width_effects, units = "mm", bg = "white", dpi = input$dpi_effects)    
    } else save(data, file = fn_downloadname_effects())
  }
  
  shinyjs::disable("bn_download_effects")
  
  # To make observeEvent watch more than one input
  toListen <- reactive({
    list(input$plot_brush, input$plot_brush, input$width_effects, input$height_effects, input$dpi_effects)
  })
  
  observeEvent(toListen(),{
    df <- brushedPoints(qtl.data()[[2]], input$plot_brush, xvar = "x", yvar = "y.dat")
    
    if (dim(df)[1] > 0 & !is.null(loadQTL()) & !is.null(input$plot_brush) & input$width_effects > 1 & input$height_effects > 1 & input$dpi_effects > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download_effects")
    } else {
      shinyjs::disable("bn_download_effects")
    }
  })
  
  # download handler
  output$bn_download_effects <- downloadHandler(
    filename = fn_downloadname_effects,
    content = function(file) {
      fn_download_effects()
      file.copy(fn_downloadname_effects(), file, overwrite=T)
      file.remove(fn_downloadname_effects())
    }
  )
  
  # Download haplotypes
  shinyjs::disable("bn_download_haplo")
  # create filename
  fn_downloadname_haplo <- reactive({
    
    seed <- sample(1:1000,1)
    if(input$fformat_haplo=="png") filename <- paste0("haplotypes","_",seed,".png")
    if(input$fformat_haplo=="tiff") filename <- paste0("haplotypes","_",seed,".tiff")
    if(input$fformat_haplo=="jpeg") filename <- paste0("haplotypes","_",seed,".jpg")
    if(input$fformat_haplo=="pdf") filename <- paste0("haplotypes","_",seed,".pdf")
    if(input$fformat_haplo=="RData") filename <- paste0("haplotypes","_",seed,".RData")
    return(filename)
  })
  
  # download 
  fn_download_haplo <- function()
  {
    p <- select_haplo(input$haplo, loadQTL()$probs, loadQTL()$selected_mks, effects.data())
    plots <- ggarrange(plotlist = p, ncol = 3, common.legend = TRUE)
    
    if(input$fformat_haplo!="RData"){
      ggsave(plots, filename = fn_downloadname_haplo(), height = input$height_haplo, 
             width = input$width_haplo, units = "mm", bg = "white", dpi = input$dpi_haplo)    
    } else save(p, file = fn_downloadname_haplo())
  }
  
  observe({
    if (input$haplo_submit & length(grep("Trait",input$haplo)) > 0 & !is.null(input$plot_brush) & input$height_haplo > 1 & input$width_haplo > 1 & input$dpi_haplo > 1) {
      Sys.sleep(1)
      # enable the download button
      shinyjs::enable("bn_download_haplo")
    } else {
      shinyjs::disable("bn_download_haplo")
    }
  })
  
  # download handler
  output$bn_download_haplo <- downloadHandler(
    filename = fn_downloadname_haplo,
    content = function(file) {
      fn_download_haplo()
      file.copy(fn_downloadname_haplo(), file, overwrite=T)
      file.remove(fn_downloadname_haplo())
    }
  )
}

## To be copied in the UI
# mod_qtl_view_ui("qtl_view_ui_1")

## To be copied in the server
# mod_qtl_view_server("qtl_view_ui_1")
