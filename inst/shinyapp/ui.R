library(shiny)
library(rhandsontable)

# Define UI for miles per gallon application
shinyUI(
  navbarPage(title="PopED - Population Experimental Design",
             collapsible = TRUE,
             tabPanel("Model Definition",fluidPage(
               title = "Hello Shiny!",
               h1("Model Definition"),
               'hh',
               h2("PK Model"),
               checkboxInput("pk_mod", label = "PK model", value = TRUE),  
               conditionalPanel(
                 condition = "input.pk_mod == true",
                 
                 ## TODO: add more structural models, create function like in PKPDmodels or with PKPDsim
                 ## allow for different parameterizations as in PKPDsim (list to change values in equations)
                 ## print out equation so that users can sse what they are using
                 selectInput("struct_PK_model", "Structural model:",
                             list(
                               "1-cpt, 1st order abs., single dose, CL param." = "ff.PK.1.comp.oral.sd.CL",
                               "1-cpt, 1st order abs., single dose, KE param." = "ff.PK.1.comp.oral.sd.KE",
                               "1-cpt, 1st order abs., multi. dose, CL param." = "ff.PK.1.comp.oral.md.CL",
                               "1-cpt, 1st order abs., multi. dose, KE param" = "ff.PK.1.comp.oral.md.KE"
                             )),
                 
                 selectInput("ruv_pk_model", "Residual error model",
                             list(
                               "Additive + Proportional" = "feps.add.prop",
                               "Proportional" = "feps.prop",
                               "Additive" = "feps.add"
                             ))
               ),
               br(),
               br(),
               h2("PD Model"),
               checkboxInput("pd_mod", label = "PD model", value = FALSE), 
               conditionalPanel(
                 condition = "input.pd_mod == true && input.pk_mod == true",
                 checkboxInput("effect_compartment", label = "Use effect compartment", value = FALSE),
                 checkboxInput("turnover", label = "Use turnover model", value = FALSE),
                 conditionalPanel(
                   condition = "input.turnover == true && input.pd_mod == true",
                   selectInput("turnover_model", label = "Choose turnover model", list(
                     "Inhibit Kin"="inhib_kin",
                     "Stimulate Kin"="stim_kin",
                     "Inhibit Kout"="inhib_kout",
                     "Stimulate Kout"="stim_kout"))
                 )
               ),
               
               conditionalPanel(
                 condition = "input.pd_mod == true",
                 
                 ## TODO: add more structural models, create function like in PKPDmodels or with PKPDsim
                 ## allow for different parameterizations as in PKPDsim (list to change values in equations)
                 ## print out equation so that users can sse what they are using
                 selectInput("struct_PD_model", "Structural PD Model:",
                             list(
                               "Emax" = "emax",
                               "Linear" = "linear"
                             )),
                 selectInput("ruv_pd_model", "Residual error model",
                             list(
                               "Additive + Proportional" = "feps.add.prop",
                               "Proportional" = "feps.prop",
                               "Additive" = "feps.add"
                             ))
               ),
               br(),
               
               # h2("Between Subject Variability Model"),
               # conditionalPanel(
               #   condition = "input.pk_mod == true",
               #   
               #   selectInput("bsv_pk_model","PK BSV",
               #               list(
               #                 "Exponential" = "exp",
               #                 "Proportional" = "prop",
               #                 "Additive" = "add",
               #                 "None" = "none"
               #               )),
               #   checkboxInput("per_pd_param", label = "Choose per PK parameter", value = FALSE)  
               # ),
               # conditionalPanel(
               #   condition = "input.pd_mod == true",
               #   
               #   selectInput("bsv_pd_model","PD BSV",
               #               list(
               #                 "Exponential" = "exp",
               #                 "Proportional" = "prop",
               #                 "Additive" = "add",
               #                 "None" = "none"
               #               )),
               #   checkboxInput("per_pd_param", label = "Choose per PD parameter", value = FALSE)  
               # ),
               
               
               #textInput("param"),
               br()
               
               
             )),
             tabPanel("Parameter Definition",
                      #uiOutput("parameter_vales"),
                      #radioButtons("useType", "Use Data Types", c("TRUE", "FALSE")),
                      rHandsontableOutput("hot"),
                      #helpText(paste0("value is ")),
                      #h3("Residual Unexplained Variability Model"),
                      # conditionalPanel(
                      #   condition = "input.bsv_per_param == false",
                      #   selectInput("bsv_model","",
                      #               list(
                      #                 "Exponential" = "exp",
                      #                 "Proportional" = "prop",
                      #                 "Additive" = "add",
                      #                 "None" = "none"
                      #               ))
                      # ),
                      #checkboxInput("bsv_per_param", label = "Choose BSV model per parameter", value = FALSE),
                      br()
             ),
             tabPanel("Design Definition",
                      textInput("num_groups", "Nunber of design groups", "1"),
                      uiOutput("group_designs"),
                      hr(),
                      #actionButton("new_group","Add a new group"),
                      ##  create plot of model 
                      ## plot_model_prediction(poped.db,IPRED=T,DV=T)
                      
                      #     checkboxInput("IPRED", "Show IPRED", FALSE),
                      #     
                      #     checkboxInput("DV", "Show DV", FALSE),
                      #     checkboxInput("separate.groups", "Separate Groups", FALSE)
                      
                      ## identifiers
                      br()
             ),
             
             #mainPanel(
             #h3(textOutput("caption")),
             #     checkboxInput("smooth", "Smooth"),
             #     conditionalPanel(
             #       condition = "input.smooth == true",
             #       selectInput("smoothMethod", "Method",
             #                   list("lm", "glm", "gam", "loess", "rlm"))
             #     ),
             
             tabPanel("Plot model/design", 
                      plotOutput("modelPlot"),
                      #submitButton("Update View"),
                      checkboxInput("IPRED", "Show IPRED", FALSE),
                      checkboxInput("DV", "Show DV", FALSE),
                      checkboxInput("separate.groups", "Separate Groups", FALSE),
                      br()
             ), 
             tabPanel("Evaluate design", verbatimTextOutput("summary"),
                      br()
             ),
             tabPanel("Optimize design", tableOutput("table"),
                      br()
             ),
             footer = list(
               hr(),
               img(src = "poped_splash.png", height = 72, width = 72),
               a(paste("PopED for R (", packageVersion("PopED"),")",sep=""), 
                 href = "http://poped.sf.net"),
               div(tags$small(
                 img(src = "", height = 2, width = 3),
                 "(c) 2014-2016, Andrew C. Hooker,",
                 tags$br(),
                 img(src = "", height = 2, width = 3),
                 "Pharmacometrics Research Group,",
                 tags$br(),
                 img(src = "", height = 2, width = 3),
                 "Uppsala University,",
                 tags$br(),
                 img(src = "", height = 2, width = 3),
                 "Sweden"
               ))
             )
             )
  )
  #plotOutput("modelPlot")
  #   footer=(
  #     br(),
  #     img(src = "poped_splash.png", height = 72, width = 72), 
  #     a(paste("PopED for R (", packageVersion("PopED"),")",sep=""), 
  #       href = "http://poped.sf.net"),
  #     h6("(c) 2014, Andrew C. Hooker, Pharmacometrics Research Group, Uppsala University, Sweden")
  #     
  #   )
  