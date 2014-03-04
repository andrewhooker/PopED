library(shiny)
#library(datasets)

# We tweak the "am" field to have nicer factor labels. Since this doesn't
# rely on any user inputs we can do this once at startup and then use the
# value throughout the lifetime of the application
#mpgData <- mtcars
#mpgData$am <- factor(mpgData$am, labels = c("Automatic", "Manual"))

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  updateModel <- reactive({
    if(input$db=="db.1") {
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.model.POPED.R")
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.design.POPED.R") 
      poped.db <- poped.db.2
      facet_scales="free"
      #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    }
    if(input$db=="db.2"){
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.design.all_in_one.POPED.R") # 4-group, add+prop, as pfim
      poped.db <- create.poped.database(warfarin.design.1.red.input())
      poped.db$ga <- rbind(50,60,70,80)
      facet_scales="fixed"
      #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    }
  })
  
  # Return the formula text for printing as a caption
  #output$caption <- renderText({
  #  "Model predictions"
  #})
  
  # Generate a plot of the requested variable against mpg and only 
  # include outliers if requested
  output$modelPlot <- renderPlot({
    if(input$db=="db.1") {
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.model.POPED.R")
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.design.POPED.R") 
      poped.db <- poped.db.2
      facet_scales="free"
      #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    }
    if(input$db=="db.2"){
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.design.all_in_one.POPED.R") # 4-group, add+prop, as pfim
      poped.db <- create.poped.database(warfarin.design.1.red.input())
      poped.db$ga <- rbind(50,60,70,80)
      facet_scales="fixed"
      #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    }
    #updateModel()
    print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales=facet_scales))
    #print(plot_model_prediction(poped.db.2,IPRED=TRUE,DV=TRUE))
  })
})