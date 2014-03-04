library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("PopED (Population Experimental Design)"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("db", "Model and design:",
                list("1-comp PK, Direct effect Emax PD" = "db.1", 
                     "Warfarin" = "db.2")),
    
    ##  create plot of model 
    ## plot_model_prediction(poped.db,IPRED=T,DV=T)
    
    checkboxInput("IPRED", "Show IPRED", FALSE),
    
    checkboxInput("DV", "Show DV", FALSE),
    checkboxInput("separate.groups", "Separate Groups", FALSE)
    
    
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    #h3(textOutput("caption")),
    
    plotOutput("modelPlot")
  )
))
