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
        titlePanel("Sweetpotato genetic map - Beauregard x Tanzania (BT)"),
        selectInput(inputId = ns('select'), label = '', choices = paste("Linkage Group", 1:15), selected = "Linkage Group 1"),
        plotOutput(ns("plot1"), height = "500px"),
        wellPanel(
          sliderInput(ns("range"), "Centimorgans", 0, 300,
                      value = c(0, 20), step = 1) , style = "padding: 6px;"
        ),
        checkboxInput(ns("op"), label = "Show SNP names?", value = TRUE),
        splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), h4("Legend"),h4("Number of SNPs per dosage"),"",h4("Summary"),"", h4("   Notes")),
        splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), 
                    includeHTML(system.file("ext/include.html", package = "viewpoly")), verbatimTextOutput("text1"), "", verbatimTextOutput("text2"), "",
                    includeHTML(system.file("ext/include2.html", package="viewpoly"))),
        splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), "", HTML("rows: Beauregard"),"","", "",""),
        splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), "", HTML("columns: Tanzania"),"","", "","")
      )
    )
  )
}
    
#' map_view Server Functions
#'
#' @noRd 
mod_map_view_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    loadMap<-reactive({
      load(file = system.file("ext/maps.rda", package = "viewpoly"))
      return(maps)
    })
    loadPh<-reactive({
      load(file = system.file("ext/phs.rda", package = "viewpoly"))
      return(list(ph.br = ph.br, ph.tz = ph.tz))
    })
    loadDose<-reactive({
      load(file = system.file("ext/ds.rda", package = "viewpoly"))
      return(list(d.br = d.br, d.tz = d.tz))
    })
    output$plot1 <- renderPlot({
      draw_map_shiny(left.lim = input$range[1], 
                     right.lim = input$range[2], 
                     ch = input$select,
                     d.br = loadDose()$d.br,
                     d.tz = loadDose()$d.tz, 
                     maps = loadMap(), 
                     ph.br = loadPh()$ph.br, 
                     ph.tz = loadPh()$ph.tz,
                     snp.names = input$op)
      output$text1 <- renderPrint({
        map_summary(left.lim = input$range[1],
                    right.lim = input$range[2],
                    ch = input$select,
                    maps = loadMap(),
                    d.br = loadDose()$d.br,
                    d.tz = loadDose()$d.tz)[1]
      })
      output$text2 <- renderPrint({
        map_summary(left.lim = input$range[1],
                    right.lim = input$range[2],
                    ch = input$select,
                    maps = loadMap(),
                    d.br = loadDose()$d.br,
                    d.tz = loadDose()$d.tz)[2:4]
      })
    })
  })
}
    
## To be copied in the UI
# mod_map_view_ui("map_view_ui_1")
    
## To be copied in the server
# mod_map_view_server("map_view_ui_1")
