library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("PopED (Population Experimental Design)"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("model", "Model:",
                list("PK: 1-cpt, 1st order abs." = "one.comp.oral",
                     "1-comp PK, Direct effect Emax PD" = "db.1", 
                     "Warfarin" = "db.2")),
    
    selectInput("ruv_model", "Residual Unexplained Variability Model:",
                list("Additive + Proportional" = "feps.add.prop")),
    conditionalPanel(
      condition = "input$model == one.comp.oral",
      textInput("xt", "Sample times:", "1,2,8,16,24" )
    )    
    ##  create plot of model 
    ## plot_model_prediction(poped.db,IPRED=T,DV=T)
    
    #     checkboxInput("IPRED", "Show IPRED", FALSE),
    #     
    #     checkboxInput("DV", "Show DV", FALSE),
    #     checkboxInput("separate.groups", "Separate Groups", FALSE)
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    #h3(textOutput("caption")),
    #     checkboxInput("smooth", "Smooth"),
    #     conditionalPanel(
    #       condition = "input.smooth == true",
    #       selectInput("smoothMethod", "Method",
    #                   list("lm", "glm", "gam", "loess", "rlm"))
    #     ),
     tabsetPanel(
       tabPanel("Plot", 
                plotOutput("modelPlot"),
                checkboxInput("IPRED", "Show IPRED", FALSE),
                checkboxInput("DV", "Show DV", FALSE),
                checkboxInput("separate.groups", "Separate Groups", FALSE)), 
       tabPanel("Summary", verbatimTextOutput("summary")),
       tabPanel("Table", tableOutput("table")))
    #plotOutput("modelPlot")
  )
))
