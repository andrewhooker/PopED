library(PopED)

library(rhandsontable)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output, session) {
  
  values = reactiveValues()
  
  data = reactive({
    if (!is.null(input$hot)) {
      DF = hot_to_r(input$hot)
      names(DF) <- c("name",
                     "pop_val" ,
                     "pop_fixed",
                     "bsv_model",
                     "variance" ,
                     "var_fixed",
                     "covariate")
    } else {
      if (is.null(values[["DF"]])){
        name <- codetools::findGlobals(eval(parse(text=input$struct_PK_model)),merge=F)$variables  
        covariate <- name %in% c("Dose","DOSE","dose","tau","TAU","Tau")
        bsv_model <- rep("exp",length(name))
        
        
        df <- data.frame(name=name,
                         pop_val = rep(1,length(name)),
                         pop_fixed=FALSE,
                         bsv_model=factor(bsv_model,levels = c("exp","add","prop","none"),ordered=TRUE),
                         variance = rep(0.09,length(name)),
                         var_fixed=FALSE,
                         covariate=covariate, 
                         stringsAsFactors = FALSE)
        
        #no_eta <- 1:length(parameter_names_ff)*FALSE
        #names(no_eta) <- parameter_names_ff
        #no_eta[parameter_names_ff %in% no_etas]  <- TRUE
        
        
        #         DF = data.frame(val = 1:10, bool = TRUE, nm = LETTERS[1:10],
        #                         dt = seq(from = Sys.Date(), by = "days", length.out = 10),
        #                         stringsAsFactors = F)
        DF = df
      } else{
        DF = values[["DF"]]
      }
    }
    
    DF[DF[,"covariate"],"variance"] <- 0
    #DF[DF[,"covariate"],"var_fixed"] <- TRUE
    DF[DF[,"covariate"],"bsv_model"] <- "none"
    
    DF[DF[,"bsv_model"]=="none","variance"] <- 0
    #DF[DF[,"bsv_model"]=="none","var_fixed"] <- TRUE

    values[["DF"]] = DF
    DF
  })
  
  output$hot <- renderRHandsontable({
    DF = data()
    if (!is.null(DF))
      rhandsontable(DF, useTypes = TRUE, 
                    #stretchH = "all", stretchV="all", 
                    overflow="visible",
                    colHeaders = c("Parameter names","Pop. value","Fix pop. value",
                                   "BSV model", "BSV Value", "Fix BSV value", "Treat as \ncovariate"),
                    highlightCol = TRUE, highlightRow = TRUE)
  })
  
  
  
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  updateDesign <- reactive({
    xt <- list()
    #a <- list()
    groupsize <- list()
    num_groups <- input$num_groups
    for(i in 1:num_groups){
      xt_txt <- input[[paste0("xt_",i)]]
      #a_txt <- input[[paste0("a_",i)]]
      groupsize_txt <- input[[paste0("groupsize_",i)]]
      xt <-  c(xt,list(eval(parse(text=paste("c(",xt_txt,")")))))
      #a <-  c(a,list(eval(parse(text=paste("c(",a_txt,")")))))
      groupsize <-  c(groupsize,list(eval(parse(text=paste("c(",groupsize_txt,")")))))
    }
    return(list(xt=xt,a=a,groupsize=groupsize))
  })
  
  output$group_designs <- renderUI({
    out <- list()
    num_groups <- input$num_groups
    for(i in 1:num_groups){
      out <- c(out,list(h2(paste0("Group ", i))))
      out <- c(out,list(textInput(paste0("groupsize_",i), 
                                  paste0("Number of individuals in group ",i,":"), "" )))
      out <- c(out,list(textInput(paste0("xt_",i), paste0("Sample times:"))))
      #       if(num_groups > 1){
      #         out <- c(out,list(actionButton(paste0("remove_group_",i),paste0("Remove Group ",i)))) 
      #         #out <- c(out,list(renderPrint({ input[[paste0("remove_group_",i)]] })))
      #       }
    }
    #out <- c(out,list(renderPrint({ input$new_group })))
    #out <- c(out,list(actionButton("new_group","Add a new group")))  
    return(as.list(out))
  })
  
  output$parameter_vales <- renderUI({
    out <- list()
    parameter_names <- codetools::findGlobals(eval(parse(text=input$struct_PK_model)),merge=F)$variables  
    df <- data.frame(par_names=parameter_names)
    df$covariate <- df$par_names %in% c("Dose","DOSE","dose","tau","TAU","Tau")
    
    out <- c(out,list(fluidRow(
      
      column(3, wellPanel(
        
        h3("Prameter")              
      )),
      
      column(3, wellPanel(
        h3("Fixed effect value")
        # This outputs the dynamic UI component
        #uiOutput("ui")
      )),
      
      column(3, wellPanel(
        h3("Random effect value")
      )
      ))))
    
    
    for(i in 1:length(df$par_names)){
      if(!df$covariate[i]){
        out <- c(out,list(fluidRow(
          
          column(3, wellPanel(
            
            h3(df$par_names[i])              
          )),
          
          column(3, wellPanel(
            textInput(df$par_names[i],NULL)
            # This outputs the dynamic UI component
            #uiOutput("ui")
          )),
          
          column(3, wellPanel(
            textInput(df$par_names[i],NULL)
          )
          ))))
        #out <- c(out,list(h3(df$par_names[i])))
        #out <- c(out,list(textInput(df$par_names[i],"Fixed effect value")))
        #out <- c(out,list(textInput(df$par_names[i],"Random effect value")))
        
        #         out <- c(out,list(selectInput("struct_PK_model", "Structural PK Model:",
        #                                       list(
        #                                         "1-cpt, 1st order abs., single dose, CL param." = "ff.PK.1.comp.oral.sd.CL",
        #                                         "1-cpt, 1st order abs., single dose, KE param." = "ff.PK.1.comp.oral.sd.KE",
        #                                         "1-cpt, 1st order abs., multi. dose, CL param." = "ff.PK.1.comp.oral.md.CL",
        #                                         "1-cpt, 1st order abs., multi. dose, KE param" = "ff.PK.1.comp.oral.md.KE"
        #                                       ))))
      }
    }
    
    
    return(as.list(out))
  })
  
  number_of_groups <- reactive({
    max_groups <- input$new_group + 1
    num_groups <- max_groups
    cat("\n\n NEW search:\n")
    for(i in 1:max_groups){
      if(length(input[[paste0("remove_group_",i)]])!=0) {
        cat("remove_group",i, input[[paste0("remove_group_",i)]], "\n")
        num_groups <- num_groups - input[[paste0("remove_group_",i)]]
      }
    }
    return(num_groups)
  })
  
  updateModel <- reactive({
    struct_pk_model <- input$struct_pk_model
    struct_pd_model <- input$struct_pd_model
    ruv_pk_model <- input$ruv_pk_model
    ruv_pd_model <- input$ruv_pd_model
    bsv_pk_model <- input$bsv_pk_model
    bsv_pd_model <- input$bsv_pd_model
    
    sfg <- build_sfg(model=input$struct_pk_model,etas=input$bsv_pk_model)
    #environment(eval(parse(text=input$struct_model)))
    #parent.env(environment())
    
    #browser()
    
    nbpop <- find.largest.index(func.str=sfg,lab="bpop") 
    #bpop_vals=c(CL=0.15, V=8, KA=1.0, Favail=1)
    #bpop_vals=c(CL=1, V=1, KA=1, Favail=1)
    #bpop_vals <- rep(1,nbpop)
    nb <- find.largest.index(func.str=sfg,lab="b")    
    
    
    
    if(input$struct_model=="ff.PK.1.comp.oral.sd.CL"){ 
      bpop_vals=c(CL=0.15, V=8, KA=1.0, Favail=1) 
      notfixed_bpop=c(1,1,1,0)
      d_vals=c(CL=0.07, V=0.02, KA=0.6) 
      sigma_vals=c(0.1,0.1)
      groupsize=32
      #xt=c( 0.5,1,2,6,24,36,72,120),
      minxt=0
      maxxt=120
      a=70
    }
    return(list(bpop=bpop_vals,d=d_vals,sigma=sigma_vals,
                notfixed_bpop=notfixed_bpop,
                sfg=sfg))
  })
  
  create_sfg <- reactive({
    build_sfg(model=input$struct_PK_model,etas=input$bsv_pk_model)
  })
  # Return the formula text for printing as a caption
  #output$caption <- renderText({
  #  "Model predictions"
  #})
  
  #get_parameters <- 
  #codetools::findGlobals(ff.PK.1.comp.oral.sd.CL,merge=F)
  
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
    #       poped.db$design$a <- rbind(50,60,70,80)
    #       facet_scales="fixed"
    #       #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #     }
    
    #model <- updateModel()
    #design <- updateDesign()
    
    
    #     ff <- function(model_switch,xt,parameters,poped.db){
    #       ##-- Model: One comp first order absorption
    #       with(as.list(parameters),{
    #         y=xt
    #         y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    #         return(list(y=y,poped.db=poped.db))
    #       })
    #     }
    #     
    #     sfg <- function(x,a,bpop,b,bocc){
    #       ## -- parameter definition function 
    #       parameters=c(CL=bpop[1]*exp(b[1]),
    #                    V=bpop[2]*exp(b[2]),
    #                    KA=bpop[3]*exp(b[3]),
    #                    Favail=bpop[4],
    #                    DOSE=a[1])
    #       return(parameters) 
    #     }
    
    #     feps <- function(model_switch,xt,parameters,epsi,poped.db){
    #       ## -- Residual Error function
    #       ## -- Proportional 
    #       returnArgs <- ff(model_switch,xt,parameters,poped.db) 
    #       y <- returnArgs[[1]]
    #       poped.db <- returnArgs[[2]]
    #       y = y*(1+epsi[,1])
    #       
    #       return(list(y=y,poped.db=poped.db)) 
    #     }
    
    # input bpop, not_fixed (for all), d_vec, sigma
    # get design variables
    # get design space
    
    
    df <- data()
    df_2 <- df[df$covariate==F,]
    bpop <- df_2[["pop_val"]]
    names(bpop) <- df_2[["name"]]
    bpop_notfixed <- !df_2[["pop_fixed"]]
    names(bpop_notfixed) <- df_2[["name"]]
    par_names <- df_2[["name"]]
    
    sfg <- build_sfg(model=NULL,
                     par_names = df[["name"]],
                     covariates = df[["covariate"]],
                     bsv_model=df[["bsv_model"]],
                     covariates=c("DOSE"),
                     etas="exp", # can be exp, prop, add
                     no_etas=c("F","Favail"),
                     env = parent.frame())
    
    browser()
    
    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1)
    bpop_notfixed <- c(CL=1, V=1, KA=1, Favail=0) 
    d_vec <- c(CL=0.07, V=0.02, KA=0.6)
    sigma <- c(prop=0.01,add=0.1)
    parameter_names <- codetools::findGlobals(eval(parse(text=input$struct_PK_model)),merge=F)$variables  
    new_bpop <- c()
    new_bpop_notfixed <- c()
    new_d_vec <- c()
    for(var in parameter_names){
      if(var %in% names(bpop)) new_bpop <- c(new_bpop,bpop[var])
      if(var %in% names(bpop_notfixed))  new_bpop_notfixed <- c(new_bpop_notfixed,bpop_notfixed[var])
      if(var %in% names(d_vec)) new_d_vec <- c(new_d_vec,d_vec[var])
    }
    
    new_sigma <- sigma
    
    design <- updateDesign()
    
    ## -- Define initial design  and design space
    poped.db <- create.poped.database(ff_file=input$struct_PK_model,
                                      #ff_file="ff",
                                      #fg_fun=model$sfg,
                                      #fg_fun=sfg,
                                      fg_fun=create_sfg(),
                                      #fError_file="feps",
                                      fError_file=input$ruv_pk_model,
                                      #bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                      bpop=new_bpop,  
                                      notfixed_bpop=new_bpop_notfixed,
                                      #d=c(CL=0.07, V=0.02, KA=0.6), 
                                      d=new_d_vec, 
                                      sigma=new_sigma,
                                      groupsize=design$groupsize[[1]],
                                      #xt=c(0.5,1,2,6,24,36,72,120),
                                      #xt=design$xt,
                                      #xt=eval(parse(text=paste("c(",input$xt,")"))),
                                      xt=design$xt[[1]],
                                      minxt=0,
                                      maxxt=120,
                                      a=70)
    #plot_model_prediction(poped.db)
    #print(plot_model_prediction(poped.db))
    
    plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups)
    #print(plot_model_prediction(poped.db.1,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #print(plot_model_prediction(poped.db.2,IPRED=TRUE,DV=TRUE))
  })
})