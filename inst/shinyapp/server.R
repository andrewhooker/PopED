library(PopED)
#source("models/MODEL.2.PK.one.comp.oral.R")

ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list(y=y,poped.db=poped.db))
  })
}

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional 
  returnArgs <- ff(model_switch,xt,parameters,poped.db) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y*(1+epsi[,1])
  
  return(list(y=y,poped.db=poped.db)) 
}


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
    ff <- input$struct_model
    feps <- input$ruv_model
    sfg <- "sfg"
    if(input$struct_model=="ff.PK.1.comp.oral.sd.CL"){ 
      feps <- "feps.prop"
      notfixed_sigma <- c(1)
      notfixed_d <- c(1,1,1)
      
      bpop_vals=c(CL=0.15, V=8, KA=1.0, Favail=1) 
      notfixed_bpop=c(1,1,1,0)
      d_vals=c(CL=0.07, V=0.02, KA=0.6) 
      sigma_vals=0.01
      groupsize=32
      #xt=c( 0.5,1,2,6,24,36,72,120),
      minxt=0
      maxxt=120
      a=70
    }
    return(list(ff=ff,sfg=sfg,feps=feps,bpop=bpop_vals,d=d_vals,sigma=sigma_vals,
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
    #     if(input$model=="db.1") {
    #       source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.model.POPED.R")
    #       source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.design.POPED.R") 
    #       poped.db <- poped.db.2
    #       facet_scales="free"
    #       #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    #     }
    #     if(input$model=="db.2"){
    #       source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.design.all_in_one.POPED.R") # 4-group, add+prop, as pfim
    #       poped.db <- create.poped.database(warfarin.design.1.red.input())
    #       poped.db$ga <- rbind(50,60,70,80)
    #       facet_scales="fixed"
    #       #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #     }
    model <- updateModel_2()
    design <- updateDesign()
    poped.db <- create.poped.database(ff_file=model$ff,
                                      fError_file=model$feps,
                                      fg_file=model$sfg,
                                      groupsize=32,
                                      sigma=model$sigma,
                                      bpop=model$bpop,  
                                      d=model$d, 
                                      xt=design$xt,
                                      maxxt=120, 
                                      minxt=0, 
                                      a=70) 
    
    poped.db.1 <- create.poped.database(ff_file="ff",
                                      fg_file="sfg",
                                      fError_file="feps.prop",
                                      bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(CL=0.07, V=0.02, KA=0.6), 
                                      sigma=0.01,
                                      groupsize=32,
                                      xt=c( 0.5,1,2,6,24,36,72,120),
                                      minxt=0,
                                      maxxt=120,
                                      a=70)
    
    
    #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales=facet_scales))
    print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))#print(plot_model_prediction(poped.db.2,IPRED=TRUE,DV=TRUE))
  })
})