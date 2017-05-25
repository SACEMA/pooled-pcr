#########################################
### Copyright Cari van Schalkwyk 2017 ###
### Released under GPL-3              ###
#########################################

library(shiny)
library(plot3D)

source("create_image.R")

shinyServer(function(input, output) {
  do_it <- eventReactive(input$button, {create_plot(input$samples, input$prevalence, input$reagents,input$type, input$min_pools, input$min_indiv) })
  
  output$text1 <- renderText({
    cost_do <- do_it()
    cost_do$cost_text
  })
  
  output$plot1 <- renderPlot({
    plot_do <- do_it()
    
    hist3D (x = plot_do$prevvec, y = plot_do$possible_ps, z = plot_do$short,
            bty = "g", phi = 20,  theta = -60,
            ylab = "pool size", xlab = list("positivity (%)", cex=2), zlab = "cost as (%) of testing individual samples", border = "black",
            ticktype = "detailed", d = 5, cex.axis =2, cex.lab=2)
    arrows3D(input$prevalence, plot_do$ps, plot_do$min_cost+50, input$prevalence, plot_do$ps, plot_do$min_cost, colvar = input$prevalence^2,
             col = "black", lwd = 10, length=0.75, add=TRUE)
    
  })
})