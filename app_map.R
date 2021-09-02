##################################################################################################################
#1. Delete application at https://www.shinyapps.io
#2. Delete rsconnect
#3. Deploy application with devtools::install_github("mmollina/mappoly") 
#4. Comment devtools::install_github("mmollina/mappoly")
# OBS: Version of mappoly used in the deploiment (github) should be the same as the one installed in your computer
##################################################################################################################
#require(devtools)
#devtools::install_github("mmollina/mappoly")
require(mappoly)
library("shiny")
library("RColorBrewer")
require("mappoly")
source("utils.R")

ui <- fluidPage(
  verticalLayout(
    titlePanel("Sweetpotato genetic map - Beauregard x Tanzania (BT)"),
    selectInput(inputId = 'select', label = '', choices = paste("Linkage Group", 1:15), selected = "Linkage Group 1"),
    plotOutput("plot1", height = "500px"),
    wellPanel(
      sliderInput("range", "Centimorgans", 0, 300,
                  value = c(0, 20), step = 1) , style = "padding: 6px;"
    ),
    checkboxInput("op", label = "Show SNP names?", value = TRUE),
    splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), h4("Legend"),h4("Number of SNPs per dosage"),"",h4("Summary"),"", h4("   Notes")),
    splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), includeHTML("include.html"), verbatimTextOutput("text1"), "", verbatimTextOutput("text2"), "",includeHTML("include2.html")),
    splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), "", HTML("rows: Beauregard"),"","", "",""),
    splitLayout(cellWidths = c("10%", "15%","5%", "15%", "5%","55%"), "", HTML("columns: Tanzania"),"","", "","")
  )
)

server <- function(input, output) {
  loadMap<-reactive({
    load(file = "maps.rda")
    return(maps)
  })
  loadPh<-reactive({
    load(file = "phs.rda")
    return(list(ph.br = ph.br, ph.tz = ph.tz))
  })
  loadDose<-reactive({
    load(file = "ds.rda")
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
}

shinyApp(ui = ui, server = server)
