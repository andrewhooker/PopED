library(shiny)
#library(datasets)
source("MODEL.2.PK.one.comp.oral.R")

# We tweak the "am" field to have nicer factor labels. Since this doesn't
# rely on any user inputs we can do this once at startup and then use the
# value throughout the lifetime of the application
#mpgData <- mtcars
#mpgData$am <- factor(mpgData$am, labels = c("Automatic", "Manual"))

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  updateDesign <- reactive({
    xt_txt <- input$xt
    xt <-  eval(parse(text=paste("c(",xt_txt,")")))
    return(list(xt=xt))
  })
  
  updateModel_2 <- reactive({
    if(input$model=="one.comp.oral"){ 
      ff <- "PK.one.cpt.oral.ff.model"
      sfg  <- "PK.one.cpt.oral.params.trans2"
      bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
      d_vals <- c(V=0.09,KA=0.09,CL=0.25^2)
      sigma_vals <- diag(c(0.04,5e-6))
      notfixed_bpop <- c(1,1,1,0)
      notfixed_sigma <- c(0,0)
      notfixed_d <- c(1,1,1)
    }
    return(list(ff=ff,sfg=sfg,feps=input$ruv_model,bpop=bpop_vals,d=d_vals,sigma=sigma_vals,
                notfixed_bpop=notfixed_bpop,notfixed_d=notfixed_d,notfixed_sigma=notfixed_sigma))
  })
  updateModel <- reactive({
    
    if(input$model=="one.comp") {
      source("Model.2.PK.one.comp.oral.R")
      input$model <- poped.db.1
      #input$facet_scales="fixed"
      #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    }
    if(input$model=="db.1") {
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.model.POPED.R")
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.design.POPED.R") 
      poped.db <- poped.db.2
      facet_scales="free"
      #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    }
    if(input$model=="db.3"){
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
    #     if(input$model=="one.comp") {
    #       source("Model.2.PK.one.comp.oral.R")
    #       poped.db <- poped.db.1
    #       facet_scales="fixed"
    #       #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    #     }
    if(input$model=="db.1") {
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.model.POPED.R")
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.design.POPED.R") 
      poped.db <- poped.db.2
      facet_scales="free"
      #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    }
    if(input$model=="db.2"){
      source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.design.all_in_one.POPED.R") # 4-group, add+prop, as pfim
      poped.db <- create.poped.database(warfarin.design.1.red.input())
      poped.db$ga <- rbind(50,60,70,80)
      facet_scales="fixed"
      #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    }
    model <- updateModel_2()
    design <- updateDesign()
    poped.db <- create.poped.database(list(),
                                      ff_file=model$ff,
                                      fError_file=model$feps,
                                      fg_file=model$sfg,
                                      groupsize=20,
                                      sigma=model$sigma,
                                      bpop=model$bpop,  
                                      d=model$d, 
                                      xt=design$xt,
                                      maxxt=336, 
                                      minxt=0, 
                                      a=c(20,24))  
    
    #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales=facet_scales))
    print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))#print(plot_model_prediction(poped.db.2,IPRED=TRUE,DV=TRUE))
  })
})