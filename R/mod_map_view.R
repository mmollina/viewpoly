#' map_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
#' @importFrom RColorBrewer brewer.pal
#' @import plotly
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
          
          column(4,
                 box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("Phenotypes"),
                 checkboxGroupInput(ns("phenotypes"),
                                    label = h4("Phenotypes"),
                                    choices = "This will be updated",
                                    selected = "This will be updated")
                 ), br(),
          ),
          column(5,
                 selectInput(inputId = ns("group"), label = p("Linkage group"), choices = 1:15, selected = 1),
                 checkboxInput(ns("op"), label = "Show SNP names", value = TRUE),
          ),
        ), hr(),
        wellPanel(
          sliderInput(ns("range"), "Map range (cM)", min = 0, max = 300,
                      value = c(0, 20), step = 1), 
          uiOutput(ns("interval"))
        ),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("LOD curve"),
            plotlyOutput(ns("plot_qtl")), 
        ), br(),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("Marker dose"),
            plotOutput(ns("plot_map"), height = "500px")
        )
      )
    )
  )
}

#' map_view Server Functions
#'
#' @importFrom shinyjs inlineCSS
#'
#' @noRd 
mod_map_view_server <- function(input, output, session, loadMap, loadJBrowse, loadQTL, parent_session){
  ns <- session$ns
  
  observe({
    # Dynamic linkage group number
    group_choices <- as.list(1:length(loadMap()$dp))
    names(group_choices) <- 1:length(loadMap()$dp)
    
    updateSelectInput(session, "group",
                      label="Linkage group",
                      choices = group_choices,
                      selected= group_choices[[1]])
    
    # Dynamic QTLs
    pheno_choices <- as.list(unique(loadQTL()$qtl_info$pheno))
    names(pheno_choices) <- unique(loadQTL()$qtl_info$pheno)
    
    updateCheckboxGroupInput(session, "phenotypes",
                             label = "Phenotypes",
                             choices = pheno_choices,
                             selected=unlist(pheno_choices)[1])
  })
  
  # Plot map
  output$plot_map <- renderPlot({
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
  
  # Plot QTL bar
  qtl.int <- reactive({
    data <- loadQTL()$qtl_info %>% filter(pheno %in% input$phenotypes & LG == input$group)
    
    if(dim(data)[1] == 0) stop("No QTL available in this group")
    
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
      p(divs_lst, "QTL")
    }
  })
  
  output$interval <- renderUI({ 
    qtl.int()
  })
  
  
  # Plot QTL profile
  output$plot_qtl <- renderPlotly({
    idx <- which(unique(loadQTL()$qtl_info$pheno) %in% input$phenotypes)
    pl <- plot_profile(loadQTL()$profile, loadQTL()$qtl_info, loadQTL()$selected_mks,
                       pheno.col = idx,
                       lgs.id = as.numeric(input$group),
                       range.min = input$range[1],
                       range.max = input$range[2], by_range=T)
    ggplotly(source = "qtl_profile", pl) %>% layout(legend = list(orientation = 'h', y = -0.3))
  })
  
  # Reactive to change page with click
  s <- reactive({
    event_data("plotly_click", source = "qtl_profile")
  })
  
  observeEvent(s(), {
    updateNavbarPage(session = parent_session, inputId = "viewpoly", selected = "qtl")
  })
}

## To be copied in the UI
# mod_map_view_ui("map_view_ui_1")

## To be copied in the server
# mod_map_view_server("map_view_ui_1")
