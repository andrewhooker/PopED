context("Design Evaluation")

test_that("Design evaluation works", {
  
  source("examples_fcn_doc/warfarin_basic.R")
  source("examples_fcn_doc/examples_evaluate_design.R")

  expect_output(str(evaluate_design(poped.db)),"List of 3")
})

test_that("Power evaluation works", {
  
  source("examples_fcn_doc/examples_evaluate_power.R")
  

  expect_error(evaluate_power(poped.db,bpop_idx = 4))

  expect_equivalent(evaluate_power(poped.db_2,bpop_idx = 4)$power$power_pred,
                    91.52704,
                    tolerance=1e-5)
  
  # does combining the fim from groups work?
  n_per_group = poped.db_2$design$groupsize
  n_tot <- sum(n_per_group)
  props = c(n_per_group/n_tot)
  norm_group_fim <- extract_norm_group_fim(poped.db_2)
  expect_equivalent(evaluate_design(poped.db_2)$ofv,
                    log(det(combine_norm_group_fim(norm_group_fim,props*n_tot))))
  
  # How about with a prior
  ff <- function(model_switch,xt,parameters,poped.db){
    with(as.list(parameters),{
      y=xt
      N = floor(xt/TAU)+1
      y=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
        (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
           exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
      return(list( y=y,poped.db=poped.db))
    })
  }
  
  sfg <- function(x,a,bpop,b,bocc){
    parameters=c( V=bpop[1]*exp(b[1]),
                  KA=bpop[2]*exp(b[2]),
                  CL=bpop[3]*exp(b[3])*bpop[5]^x[1], # add covariate for pediatrics
                  Favail=bpop[4],
                  isPediatric = x[1],
                  DOSE=a[1],
                  TAU=a[2])
    return( parameters ) 
  }
  
  ## -- Define design and design space for adults (isPediatric = 0)
  ## Two arms, 5 time points
  poped.db <- create.poped.database(ff_fun=ff,
                                    fg_fun=sfg,
                                    fError_fun=feps.add.prop,
                                    bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,pedCL=0.8), 
                                    notfixed_bpop=c(1,1,1,0,1),
                                    d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                    sigma=c(0.04,5e-6),
                                    notfixed_sigma=c(0,0),
                                    m=2,
                                    groupsize=20,
                                    xt=c( 1,8,10,240,245),
                                    bUseGrouped_xt=1,
                                    a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                    x=list(isPediatric = 0)
  )

  ## Define pediatric model/design (isPediatric = 1)
  ## One arm, 4 time points only
  poped.db.ped <- create.poped.database(ff_fun=ff,
                                        fg_fun=sfg,
                                        fError_fun=feps.add.prop,
                                        bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,pedCL=0.8), 
                                        notfixed_bpop=c(1,1,1,0,1),
                                        d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                        sigma=c(0.04,5e-6),
                                        notfixed_sigma=c(0,0),
                                        m=1,
                                        groupsize=6,
                                        xt=c( 1,2,6,240),
                                        bUseGrouped_xt=1,
                                        a=list(c(DOSE=40,TAU=24)),
                                        x=list(isPediatric = 1)
  )
  
  poped.db.all <- create.poped.database(
    poped.db.ped,
    prior_fim = evaluate.fim(poped.db)
  )
  
  expect_equivalent(evaluate_design(poped.db.all)$ofv,
                    log(det(evaluate_design(poped.db.all)$fim+evaluate.fim(poped.db))))
  
  evaluate_power(poped.db.all, bpop_idx=5, h0 = 1)

  
  n_per_group = poped.db.all$design$groupsize
  n_tot <- sum(n_per_group)
  props = c(n_per_group/n_tot)
  norm_group_fim <- extract_norm_group_fim(poped.db.all)
  expect_equivalent(props*n_tot*norm_group_fim[[1]],evaluate_design(poped.db.all)$fim)
  expect_equivalent(evaluate_design(poped.db.all)$ofv,
                    log(det(combine_norm_group_fim(norm_group_fim,props*n_tot)+evaluate.fim(poped.db))))
  
})

test_that("Shrinkage evaluation works", {
  
  ## following example 1.4.2 in PFIM
  sfg <- function(x,a,bpop,b,bocc){
    ## -- parameter definition function 
    parameters=c(KA=bpop[1]*exp(b[1]),
                 K=bpop[2]*exp(b[2]),
                 V=bpop[3]*exp(b[3]),
                 DOSE=a[1])
    return(parameters) 
  }
  
  ff <- function(model_switch,xt,parameters,poped.db){
    ##-- Model: One comp first order absorption
    with(as.list(parameters),{
      y<-(DOSE/V*KA/(KA-K)*(exp(-K*xt)-exp(-KA*xt)))
      return(list(y=y,poped.db=poped.db))
    })
  }
  
  ## -- Residual unexplained variablity (RUV) function
  ## -- Additive + Proportional  
  feps <- function(model_switch,xt,parameters,epsi,poped.db){
    returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
    f <- returnArgs[[1]]
    y = f + (0.5 + 0.15*f)*epsi[,1]
    
    return(list( y= y,poped.db =poped.db )) 
  }
  
  ## -- Define initial design  and design space
  poped.db <- create.poped.database(ff_fun=ff,
                                    fg_fun=sfg,
                                    fError_fun =feps,
                                    bpop=c(KA=2, K=0.25, V=15), 
                                    d=c(KA=1, V=0.25,0.1), 
                                    sigma=c(1),
                                    notfixed_sigma = c(0),
                                    groupsize=500,
                                    xt=c( 1,3,8),
                                    minxt=0,
                                    maxxt=10,
                                    a=100)
  
  
  shr_out <- shrinkage(poped.db)
  expect_equivalent(c(shr_out[1,1:3]),c(.5035130,.3673427,.4244184),tolerance=1e-5)
  
  ff <- function(model_switch,xt,parameters,poped.db){
    with(as.list(parameters),{
      y=xt
      MS <- model_switch
      
      # PK model
      CONC = DOSE/V*exp(-CL/V*xt) 
      
      # PD model
      EFF = E0 + CONC*EMAX/(EC50 + CONC)
      
      y[MS==1] = CONC[MS==1]
      y[MS==2] = EFF[MS==2]
      
      return(list( y= y,poped.db=poped.db))
    })
  }
  
  ## -- parameter definition function
  sfg <- function(x,a,bpop,b,bocc){
    parameters=c( 
      CL=bpop[1]*exp(b[1])  ,
      V=bpop[2]*exp(b[2])	,
      E0=bpop[3]*exp(b[3])	,
      EMAX=bpop[4],
      EC50=bpop[5]*exp(b[4])	,
      DOSE=a[1]
    )
    return( parameters ) 
  }
  
  ## -- Residual Error function
  ## -- Proportional PK + additive PD
  feps <- function(model_switch,xt,parameters,epsi,poped.db){
    returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
    y <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    MS <- model_switch
    
    prop.err <- y*(1+epsi[,1])
    add.err <- y+epsi[,2]
    
    y[MS==1] = prop.err[MS==1]
    y[MS==2] = add.err[MS==2]
    
    return(list( y= y,poped.db =poped.db )) 
  }
  
  poped.db <- create.poped.database(ff_fun=ff,
                                    fError_fun=feps,
                                    fg_fun=sfg,
                                    groupsize=20,
                                    m=3,
                                    sigma=diag(c(0.15,0.015)),
                                    bpop=c(CL=0.5,V=0.2,E0=1,EMAX=1,EC50=1),  
                                    d=c(CL=0.09,V=0.09,E0=0.04,EC50=0.09), 
                                    xt=c( 0.33,0.66,0.9,5,0.1,1,2,5),
                                    bUseGrouped_xt=1,
                                    model_switch=c( 1,1,1,1,2,2,2,2),
                                    minxt=0,
                                    maxxt=5,
                                    ourzero = 0,
                                    a=list(c(DOSE=0),c(DOSE=1),c(DOSE=2)),
                                    maxa=c(DOSE=10),
                                    mina=c(DOSE=0))
  
  

  shrk <- shrinkage(poped.db)

  
})
