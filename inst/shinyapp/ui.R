library(shiny)
library(rhandsontable)

offset <- 0

# Define UI for miles per gallon application
shinyUI(
  navbarPage(title="PopED - Population Experimental Design",
             collapsible = TRUE,
             tabPanel("Model Definition",fluidPage(
               fluidRow(
                 column(9, offset = offset,
                        #title = "Model Definition",
                        #titlePanel("Model Definition"),
                        h1("Structural Model"),
                        helpText("Define the structural pharmacokinetic (PK) and/or the pharmacodynamic (PD) model.")
                 )
               ),
               fluidRow(
                 column(3, offset = offset,
                        #checkboxInput("pk_mod", label = "PK model", value = TRUE),  
                        #h4("PK Model"),
                        selectInput("struct_PK_model", "Structural PK model:",
                                    choices = list(
                                      "None"="NULL",
                                      "1-compartment" = "ff.PK.1.comp.iv",
                                      "1-compartment, 1st-order absorption" = "ff.PK.1.comp.oral",
                                      "2-compartment" = "ff.PK.2.comp.iv",
                                      "2-compartment, 1st-order absorption" = "ff.PK.2.comp.oral",
                                      "3-compartment" = "ff.PK.3.comp.iv",
                                      "3-compartment, 1st-order absorption" = "ff.PK.3.comp.oral"
                                    ),
                                    selected = "ff.PK.1.comp.oral")
                        
                 ),
                 #),
                 #fluidRow(
                 column(3, offset = 0,#offset+1,
                        #checkboxInput("pk_mod", label = "PK model", value = TRUE),  
                        #h4("PK Model"),
                        selectInput("param_PK_model", "PK model parameterization:",
                                    choices = list(
                                      "Clearace and Volume"="CL",
                                      "Rate constants" = "KE"
                                    ),
                                    selected = "cl")
                        
                 )
               ),
               #br(),
               fluidRow(
                 column(3, offset = offset,
                        #h4("PD Model"),
                        #checkboxInput("pd_mod", label = "PD model", value = FALSE),
                        selectInput("struct_PD_model", "Structural PD Model:",
                                    list(
                                      "None" = "NULL",
                                      "Linear" = "linear",
                                      "Emax" = "emax",
                                      "Emax with hill coefficient" = "hill"
                                    ))
                 ),
               
                 column(3, offset = 0,
                        conditionalPanel(
                          condition = "input.struct_PD_model != 'NULL' && input.struct_PK_model != 'NULL'",
                          selectInput("link_fcn", label = "Link Function", 
                                      list(
                                        "Direct Eeffect" = "direct",
                                        "Effect compartment" = "effect",
                                        "Turnover inhibit Kin"="inhib_kin",
                                        "Turnover stimulate Kin"="stim_kin",
                                        "Turnover inhibit Kout"="inhib_kout",
                                        "Turnover stimulate Kout"="stim_kout"))
                        )
                 )
               ),
               fluidRow(
                 column(9, offset = offset,
                        h1("Between Subject Variability Model"),
                        helpText("Define the between subject variability (BSV) for the PK and/or PD model.")
                 )
               ),
               rHandsontableOutput("hot2"),
               
               #uiOutput("parameter_vales"),
               fluidRow(
                 column(9, offset = offset,
                        h1("Residual Unexplained Variability Model"),
                        helpText("Define the residual unexplained variability (RUV) for the PK and/or PD model.")
                 )
               ),
               conditionalPanel(
                 condition = "input.struct_PK_model != 'NULL'",
                 ## TODO: add more structural models, create function like in PKPDmodels or with PKPDsim
                 ## allow for different parameterizations as in PKPDsim (list to change values in equations)
                 ## print out equation so that users can sse what they are using
                 
                 selectInput("ruv_pk_model", "PK RUV model",
                             list(
                               "Additive + Proportional" = "feps.add.prop",
                               "Proportional" = "feps.prop",
                               "Additive" = "feps.add"
                             ))
               ),
               conditionalPanel(
                 condition = "input.struct_PD_model != 'NULL'",
                 
                 ## TODO: add more structural models, create function like in PKPDmodels or with PKPDsim
                 ## allow for different parameterizations as in PKPDsim (list to change values in equations)
                 ## print out equation so that users can sse what they are using
                 
                 selectInput("ruv_pd_model", "PK RUV model",
                             list(
                               "Additive + Proportional" = "feps.add.prop",
                               "Proportional" = "feps.prop",
                               "Additive" = "feps.add"
                             ))
               ),
               #hr(),
               
               
               #fluidRow(box(rHandsontableOutput("hot2", height = 400))),
               
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
               
               br()
               
               
             )),
             tabPanel("Parameter Definition",
                      #uiOutput("parameter_vales"),
                      #radioButtons("useType", "Use Data Types", c("TRUE", "FALSE")),
                      
                      fluidRow(
                        column(9, offset = offset,
                               h1("Parameter definition"),
                               helpText("Define the parameter values for the PK and/or PD model.")
                        )
                      ),
                      
                      h5("Fixed effects"),
                      rHandsontableOutput("hot3"),
                      hr(),
                      h5("BSV parameters (variance units)"),
                      rHandsontableOutput("hot4"),
                      h5("Fixed BSV parameters"),
                      rHandsontableOutput("hot5"),
                      hr(),
                      
                      h5("RUV parameters (variance units)"),
                      rHandsontableOutput("hot6"),
                      h5("Fixed RUV parameters"),
                      rHandsontableOutput("hot7"),
                      
                      #rHandsontableOutput("hot"),
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
                      uiOutput("test"),
                      rHandsontableOutput("hot8"),
                      
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
               fluidRow(
                 column(6, offset = offset,
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
