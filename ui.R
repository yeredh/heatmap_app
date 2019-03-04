#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(RColorBrewer)
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  #titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      #inputPanel(
      selectInput(
        "mat_display", label = "Display",
        choices = c("lower_tri","full"),
        selected = "lower_tri"
      ),
      selectInput(
        "color_pal", label = "Color Palette",
        choices = c("pheatmap",rownames(brewer.pal.info)[brewer.pal.info$category=="seq"]),
        selected = "pheatmap"
      ),
      selectInput(
        "dat", label = "Data Set",
        choices = c("probes","genes"),
        selected = "genes"
      )
      #)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
