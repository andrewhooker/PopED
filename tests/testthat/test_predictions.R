context("Model Predictions")

test_that("model_prediction works", {
  
  source("examples_fcn_doc/examples_model_prediction.R")
  
  expect_equal(length(unique(df_2$ID)),32)
  
  expect_null(df_3$DV)
  
  expect_null(df_4$a_i)
  
  expect_equal(length(unique(df_5$Group)),2)
  expect_equal(length(unique(df_5$a_i)),2)
  expect_equal(length(unique(df_5$ID)),6)
  
  expect_equal(length(unique(df_6$Group)),2)
  expect_true(all(is.na(df_6$PRED)))
  
  expect_true(all(c("WT","AGE") %in% names(df_7)))
  
  expect_equal(length(unique(df_8$WT)),2)
  expect_equal(length(unique(df_8$AGE)),2)
  
  expect_equal(length(unique(df_9$WT)),2)
  expect_equal(length(unique(df_9$AGE)),2)
  expect_equal(length(unique(df_9$ID)),6)
  
  expect_equal(length(unique(df_10$WT)),6)
  
  expect_equal(length(unique(df_11$AGE)),6)
  
  expect_equal(length(unique(df_12$AMT)),3)
  
  expect_equal(length(unique(df_13$AMT)),2)
  
  expect_equal(length(unique(df_15$AMT[df_15$ID==1])),3)
  
  
  df_16 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_4,filename="test.csv")
  expect_true("test.csv" %in% list.files())
  unlink("test.csv")
  
  dosing_2 <- list(list(AMT=1000,RATE=NA,Time=0.5),list(AMT=3000,RATE=NA,Time=0.5),list(AMT=6000,RATE=NA,Time=0.5))
  
  expect_error(model_prediction(design=design_3,DV=T,dosing=dosing_2))

  sfg <- function(x,a,bpop,b,bocc){
    parameters=c(CL=bpop[1]*exp(b[1]),
                 V=bpop[2]*exp(b[2]),
                 KA=bpop[3]*exp(b[3]),
                 Favail=bpop[4],
                 DOSE=a[1])
    return(parameters) 
  }
  
  ## -- Define initial design  and design space
  poped.db.2 <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                    fg_fun=sfg,
                                    fError_fun=feps.add.prop,
                                    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6), 
                                    sigma=c(prop=0.01,add=1),
                                    groupsize=32,
                                    xt=c( 0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=120,
                                    a=70)
  plot_model_prediction(poped.db.2,PI=T,DV=T)#,groupsize_sim = 500)
  df_20 <- model_prediction(poped.db.2,PI=TRUE)
  expect_true(all(c("PI_l","PI_u") %in% names(df_20)))
  
  sfg.3 <- function(x,a,bpop,b,bocc){
    parameters=c(CL=bpop[1]*exp(b[1]),
                 V=bpop[2]*exp(b[2]),
                 KA=bpop[3]*exp(b[3]),
                 Favail=bpop[4],
                 DOSE=a[1],
                 TAU=a[2])
    return(parameters) 
  }
  poped.db.3 <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                      fg_fun=sfg.3,
                                      fError_fun=feps.add.prop,
                                      bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(CL=0.07, V=0.02, KA=0.6), 
                                      sigma=c(prop=0.01,add=1),
                                      groupsize=32,
                                      xt=c( 0.5,1,2,6,24,36,72,120),
                                      minxt=0,
                                      maxxt=120,
                                      a=c(DOSE=70,TAU=200))
  plot_model_prediction(poped.db.3,PI=T,DV=T)#,groupsize_sim = 500)

})

test_that("plot_model_prediction works", {

  source("examples_fcn_doc/examples_plot_model_prediction.R")
  expect_output(str(plot_model_prediction(poped.db)))
  
})
