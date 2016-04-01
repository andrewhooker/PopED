context("E-family")

test_that("ed_laplace_ofv works", {
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_ed_laplace_ofv.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")  
  sink()
  file.remove("tmp.txt","temp.R")
  
  expect_gt(output$E_ofv,1.0e+24)
  

  ## -- parameter definition function 
  sfg <- function(x,a,bpop,b,bocc){
    parameters=c( V=bpop[1]*exp(b[1]),
                  KA=bpop[2]*exp(b[2]),
                  CL=bpop[3]*exp(b[3]),
                  Favail=bpop[4],
                  DOSE=a[1],
                  TAU=a[2],
                  E0=bpop[5]*exp(b[4]),
                  IMAX=bpop[6],
                  IC50=bpop[7]
    )
    return( parameters ) 
  }
  
  ##-- Model: One comp first order absorption + inhibitory imax
  ## -- works for both mutiple and single dosing
  ff <- function(model_switch,xt,parameters,poped.db){
    with(as.list(parameters),{
      y=xt
      MS <- model_switch
      
      # PK model
      N = floor(xt/TAU)+1
      CONC = (DOSE*Favail/V)*(KA/(KA - CL/V)) * 
        (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
           exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
      
      # PD model
      EFF = E0*(1 - CONC*IMAX/(IC50 + CONC))
      
      y[MS==1] = CONC[MS==1]
      y[MS==2] = EFF[MS==2]
      
      return(list( y= y,poped.db=poped.db))
    })
  }
  
  ## -- Residual unexplained variablity (RUV) function
  ## -- Additive + Proportional  on both PK and PD measurements
  feps <- function(model_switch,xt,parameters,epsi,poped.db){
    returnArgs <- ff(model_switch,xt,parameters,poped.db) 
    y <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    MS <- model_switch
    
    pk.dv <- y*(1+epsi[,1])+epsi[,2]
    pd.dv <-  y*(1+epsi[,3])+epsi[,4]
    
    y[MS==1] = pk.dv[MS==1]
    y[MS==2] = pd.dv[MS==2]
    
    return(list( y= y,poped.db =poped.db )) 
  }
  
  bpop_vals=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,
              E0=1120,IMAX=0.807,IC50=0.0993)
  bpop_vals_ed <- cbind(zeros(length(bpop_vals),1),
                        bpop_vals,
                        zeros(length(bpop_vals),1)) 
  bpop_vals_ed["IC50",1] <- 4
  bpop_vals_ed["IC50",3] <- (bpop_vals_ed["IC50",2]*0.1)^2
  
  poped_db <- create.poped.database(ff_fun=ff,
                                      fError_fun=feps,
                                      fg_fun = sfg,
                                      groupsize=20,
                                      m=3,
                                      sigma=c(0.04,5e-6,0.09,100),
                                      bpop=bpop_vals_ed,  
                                      notfixed_bpop=c(1,1,1,0,1,1,1),
                                      notfixed_sigma=c(0,0,0,0),
                                      d=c(V=0.09,KA=0.09,CL=0.25^2,E0=0.09), 
                                      xt=c(0,0,0,240,240,0,0,0,240,240),
                                      minxt=c(0,0,0,240,240,0,0,0,240,240),
                                      maxxt=c(10,10,10,248,248,10,10,10,248,248),
                                      model_switch=c(1,1,1,1,1,2,2,2,2,2),
                                      bUseGrouped_xt=TRUE,
                                      a=list(c(DOSE=20,TAU=24),
                                             c(DOSE=40,TAU=24),
                                             c(DOSE=0,TAU=24)),
                                      maxa = c(DOSE=200,TAU=48),
                                      mina = c(DOSE=0,TAU=8),
                                      G_xt = c(1,2,3,4,5,1,2,3,4,5),
                                      ourzero = 0,
                                      discrete_xt =list(0:10,0:10,0:10,240:248,240:248,
                                                        0:10,0:10,0:10,240:248,240:248),
                                      discrete_a = list(0:50,8:48),
                                      ED_samp_size=20) 
  
  
  
  expect_warning(evaluate.e.ofv.fim(poped_db,use_laplace = T))
  
  expect_warning(expect_error(poped_optim(poped_db,opt_xt=T,parallel = T, d_switch=0,use_laplace=T)))
  
  
    
})