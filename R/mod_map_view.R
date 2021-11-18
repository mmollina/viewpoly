#' map_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
#' @importFrom plotly plotlyOutput
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
          
          column(6,
                 column(6,
                        box(width = 12, solidHeader = FALSE, collapsible = FALSE,  collapsed = FALSE,
                            pickerInput(ns("phenotypes"),
                                        label = h4("Phenotypes"),
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
                        box(width = 12, solidHeader = FALSE, collapsible = FALSE,  collapsed = FALSE,
                            selectInput(inputId = ns("group"), label = p("Linkage group"), choices = 1:15, selected = 1),
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
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = h4("QTL profile"),
            plotlyOutput(ns("plot_qtl")), 
        ), br(),
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("Haplotypes and dosages"),
            plotOutput(ns("plot_map"), height = "500px"),
            box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("Parents haplotypes"),
                DT::dataTableOutput(ns("parents_haplo"))
            ),
            box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = FALSE, status="primary", title = h4("Progeny haplotypes"),
                DT::dataTableOutput(ns("progeny_haplo"))
            )
        )
      )
    )
  )
}

#' map_view Server Functions
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly ggplotly renderPlotly
#'
#' @noRd 
mod_map_view_server <- function(input, output, session, 
                                loadExample,
                                loadMap_custom, loadMap_mappoly,
                                loadQTL_custom, loadQTL_qtlpoly, loadQTL_diaQTL, loadQTL_polyqtlR,
                                parent_session){
  ns <- session$ns
  
  loadMap = reactive({
    if(is.null(loadExample()) & is.null(loadMap_custom()) & is.null(loadMap_mappoly())){
      warning("Select one of the options in `upload` session")
      return(NULL)
    } else if(!is.null(loadMap_custom())){
      return(loadMap_custom())
    } else if(!is.null(loadMap_mappoly())){
      return(loadMap_mappoly())
    } else if(!is.null(loadExample())){
      return(loadExample()$map)
    }
  })
  
  loadQTL = reactive({
    if(is.null(loadExample()) & 
       is.null(loadQTL_custom()) & 
       is.null(loadQTL_diaQTL()) &
       is.null(loadQTL_polyqtlR())) {
      warning("Select one of the options in `upload` session")
      return(NULL)
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
  
  # Plot map
  output$plot_map <- renderPlot({
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
    } else 
      stop("Upload the QTL information in upload session to access this feature.")
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
      
      ggplotly(source = "qtl_profile", pl, tooltip=c("Trait", "Position (cM)")) %>% layout(legend = list(orientation = 'h', y = -0.3))
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
      stop("Upload map information in the upload session to access this feature.")
  })
  
}

## To be copied in the UI
# mod_map_view_ui("map_view_ui_1")

## To be copied in the server
# mod_map_view_server("map_view_ui_1")
