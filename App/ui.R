#########################################
### Copyright Cari van Schalkwyk 2017 ###
### Released under GPL-3              ###
#########################################

library(shiny)

# Define UI for application
shinyUI(fluidPage(
  
  titlePanel("Cost-efficiency of pooled PCR testing for infant HIV diagnosis"),
  
  sidebarLayout(
    sidebarPanel(
      h2("Input values"),
      p("Choose values appropriate for your lab"),
      numericInput("samples","Average number of samples tested per day:", value=50, step=1, min=0),
      sliderInput("prevalence","Estimated positivity (%):",min=0, max=10, value=4, step=0.5),
      selectInput("type", "Sample Type", 
                  choices = c("Dried blood spot", "Whole blood")),
      numericInput("reagents","Optional: Cost of reagent (in US$):", value=0, step=0.25, min=0),
      numericInput("min_pools","Minimum samples required to run batch (pooling)", value=20, step=1, min=0),
      numericInput("min_indiv","Minimum samples required to run batch (individual testing)", value=10, step=1, min=0),
      br(),
      actionButton("button","Calculate"),
      br(),
      br(),
      br(),
      
      img(src = "SACEMA_logo.png", height = 100, width = 200), img(src = "US_logo.png", height = 100, width = 200),
      br(),
      img(src = "nhls.png", height = 100, width = 100)
    ),
    mainPanel(
      h2("Output values"),
      #verbatimTextOutput("text1"),
      h3(textOutput("text1")),
      br(),
      
      plotOutput("plot1", height = 800, width = 800),
      tags$head(tags$style("#text1{color: blue;
                                 font-size: 20px;
                                
                                 }"))
      
    )
  )
))