library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("PopED (Population Experimental Design)"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    h2("Model Definition"),
    selectInput("model", "Model:",
                list("PK: 1-cpt, 1st order abs." = "one.comp.oral")),
    
    selectInput("ruv_model", "Residual Unexplained Variability Model:",
                list("Additive + Proportional" = "feps.add.prop")),
    h2("Design Definition"),
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
       tabPanel("Plot model/design", 
                plotOutput("modelPlot"),
                #submitButton("Update View"),
                checkboxInput("IPRED", "Show IPRED", FALSE),
                checkboxInput("DV", "Show DV", FALSE),
                checkboxInput("separate.groups", "Separate Groups", FALSE)), 
       tabPanel("Evaluate design", verbatimTextOutput("summary")),
       tabPanel("Optimize design", tableOutput("table")))
    #plotOutput("modelPlot")
  )
))
