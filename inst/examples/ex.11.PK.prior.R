library(PopED)

##-- Model: One comp first order absorption
## -- Analytic solution for both mutiple and single dosing
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

## -- parameter definition function 
## -- names match parameters in function ff
## -- note, covariate on clearance for pediatrics (using isPediatric x[1])
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

## -- Residual unexplained variablity (RUV) function
## -- Additive + Proportional  
feps <- function(model_switch,xt,parameters,epsi,poped.db){
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list( y= y,poped.db =poped.db )) 
}

## -- Define design and design space for adults (isPediatric = 0)
## Two arms, 5 time points
poped.db <- create.poped.database(ff_fun="ff",
                                  fg_fun="sfg",
                                  fError_fun="feps",
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,pedCL=0.8), 
                                  notfixed_bpop=c(1,1,1,1,1),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(0,0),
                                  m=2,
                                  groupsize=20,
                                  xt=c( 1,2,8,240,245),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                  x=list(isPediatric = 0)
)

##  create plot of model without variability
plot_model_prediction(poped.db, model_num_points = 300)

## evaluate initial design
evaluate_design(poped.db)

## To store FIM from adult design - need intermediate result
outAdult <- calc_ofv_and_fim(poped.db)
get_rse(outAdult$fim, poped.db)

## Define pediatric model/design
## One arm, 4 time points only - insufficient on its own
poped.db.ped <- create.poped.database(ff_fun="ff",
                                  fg_fun="sfg",
                                  fError_fun="feps",
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,pedCL=0.8), 
                                  notfixed_bpop=c(1,1,1,1,1),
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

##  create plot of model without variability 
plot_model_prediction(poped.db.ped, model_num_points = 300)

## evaluate design of pediatrics only - insufficient
evaluate_design(poped.db.ped)

## Add adult prior
poped.db.all <- create.poped.database(
  poped.db.ped,
  prior_fim = outAdult$fim
)

## evaluate design using prior FIM from adults
evaluate_design(poped.db.all)
