library(PopED)

library(rhandsontable)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output, session) {
  
  
  
  
  model_name <- reactive({
    mod_name <- NULL
    if(input$struct_PK_model!="NULL") mod_name <- paste(input$struct_PK_model,"sd",input$param_PK_model,sep=".")
    if(input$struct_PD_model!="NULL"){
      if(input$struct_PK_model=="NULL") mod_name <- paste0(input$struct_PD_model)
      if(input$struct_PK_model!="NULL") mod_name <- paste(mod_name,input$link_fcn,input$struct_PD_model,sep=".")
    }
    return(mod_name)
  })
  
  values = reactiveValues()
  setHot = function(x) values[["hot"]] = x
  
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
        mod_name <- model_name()
        name <- codetools::findGlobals(eval(parse(text=mod_name)),merge=F)$variables  
        covariate <- name %in% c("Dose","DOSE","dose","tau","TAU","Tau")
        bsv_model <- rep("exp",length(name))
        
        
        df <- data.frame(name=name,
                         pop_val = runif(length(name)),
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
    DF[DF[,"variance"]==0,"var_fixed"] <- TRUE
    
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
                                   "BSV model", "BSV Value", "Fix BSV value", "Treat as \ndesign variable")
                    #highlightCol = TRUE, highlightRow = TRUE
      )
  })
  
  
  param_names <- reactive({
    mod_name <- model_name()
    name <- codetools::findGlobals(eval(parse(text=mod_name)),merge=F)$variables  
    names_par <- name[!name %in% c("Dose","DOSE","dose","tau","TAU","Tau")]
  })
  
  
  #values = reactiveValues()
  setHot2 = function(x) values[["hot2"]] = x
  
  output$hot2 = renderRHandsontable({
    if (!is.null(input$hot2)) {
      DF = hot_to_r(input$hot2)
    } else {
      
      par_name <- param_names()
      bsv_model <- rep("Exponential",length(par_name))
      
      
      df <- data.frame(name=par_name,
                       bsv_model=factor(bsv_model,levels = c("Exponential","Additive","Proportional","None"),ordered=TRUE),
                       stringsAsFactors = FALSE)
      
      df$bsv_model[df$name %in% c("Favail","F")] <- "None"
      DF = df 
    }
    
    setHot2(DF)
    rhandsontable(DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible")
  })
  
  setHot3 = function(x) values[["hot3"]] = x
  
  output$hot3 = renderRHandsontable({
    if (!is.null(input$hot3)) {
      DF = hot_to_r(input$hot3)
    } else {
      par_names <- param_names()
      
      df <- data.frame(name=par_names,
                       value = runif(length(par_names)),
                       fixed=FALSE,
                       stringsAsFactors = FALSE)
      
      df$fixed[df$name %in% c("Favail","F")] <- TRUE
      df$value[df$name %in% c("Favail","F")] <- 1
      DF = df 
    }
    
    setHot3(DF)
    rhandsontable(DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible")
  })
  
  setHot4 = function(x) values[["hot4"]] = x
  
  output$hot4 = renderRHandsontable({
    
    par_names <- values[["hot2"]] %>% dplyr::filter(bsv_model!="None") %>% dplyr::select(name)
    
    if (!is.null(input$hot4)) {
      MAT = hot_to_r(input$hot4)
      update <- FALSE
      if(length(dimnames(MAT)[[1]])!=length(par_names[[1]])){
        update <- TRUE
      } else {
        if(any(dimnames(MAT)[[1]]!=par_names[[1]])) update <- TRUE        
      }
      if(update){
        MAT1 <-  zeros(nrow(par_names))
        diag(MAT1) <- 0.09
        dimnames(MAT1) <- c(par_names,par_names)
        old_names <- dimnames(MAT)[[1]]
        still_here_old_names <- old_names[old_names %in% par_names[[1]]]
        MAT1[still_here_old_names,still_here_old_names] <- MAT[still_here_old_names,still_here_old_names]
        MAT <- MAT1
      }
    } else {
      
      #par_names <- values[["hot2"]] %>% dplyr::filter(bsv_model!="None") %>% dplyr::select(name)
      MAT <-  zeros(nrow(par_names))
      diag(MAT) <- 0.09
      dimnames(MAT) <- c(par_names,par_names)
      
    }
    
    setHot4(MAT)
    rhandsontable(MAT) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible") %>%
      hot_cols(renderer = "
               function (instance, td, row, col, prop, value, cellProperties) {
               Handsontable.renderers.TextRenderer.apply(this, arguments);
               if (row == col) {
               td.style.background = 'lightgrey';
               } else if (col > row) {
               td.style.background = 'grey';
               td.style.color = 'grey';
               } else if (value != 0) {
               td.style.background = 'lightgreen';
               } else if (value > 0.75) {
               td.style.background = 'lightgreen';
               }
               }")
  })
  
  setHot5 = function(x) values[["hot5"]] = x
  
  output$hot5 = renderRHandsontable({
    
    bsv_parameters <- values[["hot4"]]
    
    if (!is.null(input$hot5)) {
      MAT = hot_to_r(input$hot5)
      update <- FALSE
      if(length(MAT)!=length(bsv_parameters)){
        update <- TRUE
      } else {
        if(any(dimnames(MAT)[[1]]!=dimnames(bsv_parameters)[[1]])) update <- TRUE        
      }
      if(update){
        MAT1 <-  bsv_parameters*FALSE
        MAT1[bsv_parameters==0] <- TRUE
        MAT1 <- as.data.frame(MAT1)
        MAT1 <- sapply(MAT1,as.logical)
        MAT1[upper.tri(MAT1)] <- NA
        MAT1 <- data.frame(MAT1,stringsAsFactors = FALSE)
        rownames(MAT1) <- names(MAT1)
        
        new_names <- dimnames(MAT1)[[1]] 
        old_names <- dimnames(MAT)[[1]]
        still_here_old_names <- old_names[old_names %in% new_names]
        MAT1[still_here_old_names,still_here_old_names] <- MAT[still_here_old_names,still_here_old_names]
        MAT <- MAT1
      }
    } else {
      
      #par_names <- values[["hot2"]] %>% dplyr::filter(bsv_model!="None") %>% dplyr::select(name)
      MAT <-  bsv_parameters*FALSE
      MAT[bsv_parameters==0] <- TRUE
      MAT <- as.data.frame(MAT)
      MAT <- sapply(MAT,as.logical)
      MAT[upper.tri(MAT)] <- NA
      MAT <- data.frame(MAT,stringsAsFactors = FALSE)
      rownames(MAT) <- names(MAT)
      
      
    }
    
    setHot5(MAT)
    rhandsontable(MAT) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible")
  })
  
  setHot6 = function(x) values[["hot6"]] = x
  
  output$hot6 = renderRHandsontable({
    par_names <- c()
    if(input$struct_PK_model!="NULL"){
      if (input$ruv_pk_model=="feps.add.prop") par_names <- c(par_names,"PK_prop","PK_add")
      if (input$ruv_pk_model=="feps.prop") par_names <- c(par_names,"PK_prop")
      if (input$ruv_pk_model=="feps.add") par_names <- c(par_names,"PK_add")
    }
    if(input$struct_PD_model!="NULL"){
      if (input$ruv_pd_model=="feps.add.prop") par_names <- c(par_names,"PD_prop","PD_add")
      if (input$ruv_pd_model=="feps.prop") par_names <- c(par_names,"PD_prop")
      if (input$ruv_pd_model=="feps.add") par_names <- c(par_names,"PD_add")
    }
    
    if (!is.null(input$hot6)) {
      MAT = hot_to_r(input$hot6)
      update <- FALSE
      if(length(dimnames(MAT)[[1]])!=length(par_names)){
        update <- TRUE
      } else {
        if(any(dimnames(MAT)[[1]]!=par_names)) update <- TRUE        
      }
      if(update){
        MAT1 <-  zeros(length(par_names))
        diag(MAT1) <- 0.01
        dimnames(MAT1) <- c(list(par_names),list(par_names))
        old_names <- dimnames(MAT)[[1]]
        still_here_old_names <- old_names[old_names %in% par_names]
        MAT1[still_here_old_names,still_here_old_names] <- MAT[still_here_old_names,still_here_old_names]
        MAT <- MAT1
      }
    } else {
      
      
      
      #par_names <- values[["hot2"]] %>% dplyr::filter(bsv_model!="None") %>% dplyr::select(name)
      MAT <-  zeros(length(par_names))
      diag(MAT) <- 0.01
      dimnames(MAT) <- c(list(par_names),list(par_names))
      
    }
    
    setHot6(MAT)
    rhandsontable(MAT) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible",stretchH = "right") %>%
      hot_cols(renderer = "
               function (instance, td, row, col, prop, value, cellProperties) {
               Handsontable.renderers.TextRenderer.apply(this, arguments);
               if (row == col) {
               td.style.background = 'lightgrey';
               } else if (col > row) {
               td.style.background = 'grey';
               td.style.color = 'grey';
               } else if (value != 0) {
               td.style.background = 'lightgreen';
               } else if (value > 0.75) {
               td.style.background = 'lightgreen';
               }
               }")
  })
  
  setHot7 = function(x) values[["hot7"]] = x
  
  output$hot7 = renderRHandsontable({
    
    bsv_parameters <- values[["hot6"]]
    
    if (!is.null(input$hot7)) {
      MAT = hot_to_r(input$hot7)
      update <- FALSE
      if(length(MAT)!=length(bsv_parameters)){
        update <- TRUE
      } else {
        if(any(dimnames(MAT)[[1]]!=dimnames(bsv_parameters)[[1]])) update <- TRUE        
      }
      if(update){
        MAT1 <-  bsv_parameters*FALSE
        MAT1[bsv_parameters==0] <- TRUE
        MAT1 <- as.data.frame(MAT1)
        MAT1 <- sapply(MAT1,as.logical)
        MAT1[upper.tri(MAT1)] <- NA
        MAT1 <- data.frame(MAT1,stringsAsFactors = FALSE)
        rownames(MAT1) <- names(MAT1)
        
        new_names <- dimnames(MAT1)[[1]] 
        old_names <- dimnames(MAT)[[1]]
        still_here_old_names <- old_names[old_names %in% new_names]
        MAT1[still_here_old_names,still_here_old_names] <- MAT[still_here_old_names,still_here_old_names]
        MAT <- MAT1
      }
    } else {
      
      #par_names <- values[["hot2"]] %>% dplyr::filter(bsv_model!="None") %>% dplyr::select(name)
      MAT <-  bsv_parameters*FALSE
      MAT[bsv_parameters==0] <- TRUE
      MAT <- as.data.frame(MAT)
      MAT <- sapply(MAT,as.logical)
      MAT[upper.tri(MAT)] <- NA
      MAT <- data.frame(MAT,stringsAsFactors = FALSE)
      names(MAT) <- colnames(bsv_parameters)
      rownames(MAT) <- names(MAT)
      
      
    }
    
    setHot7(MAT)
    rhandsontable(MAT) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible") 
  })
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  updateDesign <- reactive({
    xt <- list()
    groupsize <- c()
    
    DF <- data()
    cov_names <- DF[DF["covariate"]==T,"name"]
    a <- list()
    if(length(cov_names)==0) a <- NULL
    
    num_groups <- input$num_groups
    
    for(i in 1:num_groups){
      xt_txt <- input[[paste0("xt_",i)]]
      xt <-  c(xt,list(eval(parse(text=paste("c(",xt_txt,")")))))
      
      groupsize_txt <- input[[paste0("groupsize_",i)]]
      groupsize <-  c(groupsize,eval(parse(text=groupsize_txt)))
      
      # find covariates
      if(length(cov_names)!=0){
        cov_vals <- c()
        for(j in cov_names){
          a_txt <- input[[paste0(j,"_",i)]]          
          cov_vals <-  c(cov_vals,eval(parse(text=paste0(j,"=",a_txt))))
        }
        a <- c(a,cov_vals)
      }
    }
    return(list(xt=xt,a=a,groupsize=groupsize))
  })
  
  get_dose_type <- reactive({
    dose_type <- input$dose_type
    return(dose_type)
  })
  
  
  setHot8 = function(x) values[["hot8"]] = x
  
  output$hot8 = renderRHandsontable({
    if (!is.null(input$hot8)) {
      DF = hot_to_r(input$hot8)
    } else {
      
      
      df <- data.frame(group = 1L, 
                         amount = 20,
                       time = 0,
                       duration = 0,
                       n = 1L,
                       tau = 0,
                       stringsAsFactors = FALSE)
      
      DF = df 
    }
    
    setHot8(DF)
    rhandsontable(DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, overflow="visible")
  })
  
  output$group_designs <- renderUI({
    out <- list()
    num_groups <- input$num_groups
    DF <- data()
    #if(any(DF$covariate))
      for(i in 1:num_groups){
        out <- c(out,list(h3(paste0("Group ", i))))
        out <- c(out,list(textInput(paste0("groupsize_",i), 
                                    paste0("Number of individuals in group ",i,":"), "" )))
        if(input$struct_PK_model!="NULL"){
          out <- c(out,list(
            #wellPanel(
            h3(paste0("Regimen")),
            textInput(paste0("amt_",i),
                      paste0("Dose amount(s):")),
            textInput(paste0("d_time_",i),
                      paste0("Dose time(s):"),
                      value="0"),
            selectInput(paste0("dose_type_",i), "Dose type",
                        list(
                          "Bolus" = "bolus",
                          "Infusion" = "infusion"
                        ))
            #)
            
            # conditionalPanel(
            #   condition = "input.dose_type == 'infusion'",
            #   sliderInput("breakCount", "Break Count", min=1, max=1000, value=10)
            # )
          ))
          #if(!is.null(input$dose_type)){
          out <- c(out,list(
            conditionalPanel(
              condition = paste0("input.dose_type_",i," == 'infusion'"),
              textInput(paste0("inf_dur_",i),
                        paste0("Infusion duration(s):"),
                        value="")
              )))
              
          #   if(input$dose_type=="bolus") textInput(paste0("amt33_",i),paste0("dooo amount"))
          #}
          #if(get_dose_type()=="bolus") out <- c(out,list(h3(paste0("Group ", i))))
          
        }
        
        if(input$struct_PK_model!="NULL"){
          out <- c(out,list(textInput(paste0("xt_pk_",i), paste0("PK Sample times:"))))
        }
        if(input$struct_PD_model!="NULL"){
          out <- c(out,list(textInput(paste0("xt_pd_",i), paste0("PD Sample times:"))))
        }
        if(any(DF$covariate)){
          cov_names <- DF[DF["covariate"]==T,"name"]
          names_par <- cov_names[!cov_names %in% c("Dose","DOSE","dose","tau","TAU","Tau")]
          for(j in names_par){
            out <- c(out,list(textInput(paste0(j,"_",i),paste0(j,":"))))
          }
        }
        #       if(num_groups > 1){
        #         out <- c(out,list(actionButton(paste0("remove_group_",i),paste0("Remove Group ",i)))) 
        #         #out <- c(out,list(renderPrint({ input[[paste0("remove_group_",i)]] })))
        #       }
      }
    #out <- c(out,list(renderPrint({ input$new_group })))
    #out <- c(out,list(actionButton("new_group","Add a new group")))  
    return(as.list(out))
  })
  
  # result <- list()
  # for(i in 1:input$num_groups){
  #   test <- renderUI({
  #     out <- list()
  #     if(input$dose_type=="bolus") out <- c(out,list(h3(paste0("Group "))))
  #     return(as.list(out))
  #   })
  #   
  #   test2 <- renderUI({
  #     out <- list()
  #     if(input$dose_type=="bolus") out <- c(out,list(h3(paste0("Group "))))
  #     return(as.list(out))
  #   })
  #   result <- c(result, test,test2)
  # }
  # output$test <- result
  # 
  
  output$parameter_vales <- renderUI({
    out <- list()
    parameter_names <- codetools::findGlobals(eval(parse(text="ff.PK.1.comp.oral.sd.CL")),merge=F)$variables  
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
    # df_2 <- df[df$covariate==F,]
    # bpop <- df_2[["pop_val"]]
    # names(bpop) <- df_2[["name"]]
    # bpop_notfixed <- !df_2[["pop_fixed"]]
    # names(bpop_notfixed) <- df_2[["name"]]
    # par_names <- df_2[["name"]]

    
    sfg <- build_sfg(model=NULL,
                     par_names = df[["name"]],
                     covariates = df[["name"]][df[["covariate"]]],
                     etas=df[["bsv_model"]],
                     #covariates=c("DOSE"),
                     #etas="exp", # can be exp, prop, add
                     no_etas=NULL,
                     env = parent.frame())
    
    
    #bpop=c(CL=0.15, V=8, KA=1.0, Favail=1)
    bpop=df[df$covariate==F,"pop_val"]
    names(bpop) <- df[df$covariate==F,"name"]
    #bpop_notfixed <- c(CL=1, V=1, KA=1, Favail=0) 
    bpop_notfixed <- !df[df$covariate==F,"pop_fixed"]
    names(bpop_notfixed) <- df[df$covariate==F,"name"]
    #d_vec <- c(CL=0.07, V=0.02, KA=0.6)
    d_vec <- df[df$bsv_model!="none","variance"]
    names(d_vec) <- df[df$bsv_model!="none","name"]
    d_notfixed <- !df[df$bsv_model!="none","var_fixed"]
    names(d_notfixed) <- df[df$bsv_model!="none","name"]
    #sigma <- c(prop=0.01,add=0.1)
    sigma <- c(prop=0.01,add=0.1)
    
    # parameter_names <- codetools::findGlobals(eval(parse(text=input$struct_PK_model)),merge=F)$variables  
    # new_bpop <- c()
    # new_bpop_notfixed <- c()
    # new_d_vec <- c()
    # for(var in parameter_names){
    #   if(var %in% names(bpop)) new_bpop <- c(new_bpop,bpop[var])
    #   if(var %in% names(bpop_notfixed))  new_bpop_notfixed <- c(new_bpop_notfixed,bpop_notfixed[var])
    #   if(var %in% names(d_vec)) new_d_vec <- c(new_d_vec,d_vec[var])
    # }
    
    # new_sigma <- sigma
    
    design <- updateDesign()
    
    ## -- Define initial design  and design space
    
    poped.db <- create.poped.database(ff_file=input$struct_PK_model,
                                      #ff_file="ff",
                                      #fg_fun=model$sfg,
                                      #fg_fun=sfg,
                                      fg_fun=sfg,
                                      #fError_file="feps",
                                      fError_file=input$ruv_pk_model,
                                      #bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                      bpop=bpop,  
                                      notfixed_bpop=bpop_notfixed,
                                      #d=c(CL=0.07, V=0.02, KA=0.6), 
                                      d=d_vec, 
                                      notfixed_d = d_notfixed,
                                      sigma=sigma,
                                      groupsize=design$groupsize,
                                      #xt=c(0.5,1,2,6,24,36,72,120),
                                      xt=design$xt,
                                      #xt=eval(parse(text=paste("c(",input$xt,")"))),
                                      #xt=design$xt[[1]],
                                      minxt = 0,
                                      maxxt = 120,
                                      a=design$a)
    #plot_model_prediction(poped.db)
    #print(plot_model_prediction(poped.db))
    
    plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups)
    #print(plot_model_prediction(poped.db.1,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #print(plot_model_prediction(poped.db.2,IPRED=TRUE,DV=TRUE))
  })
})