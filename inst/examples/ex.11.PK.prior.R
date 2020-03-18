library(PopED)
# This example shows how to include a prior FIM into the design evaluation.
# We look at PK assessment in pediatrics, where we are mainly interested 
# in assessing if there is a substantial difference of more than 20% in 
# clearance (CL) between children and adults.

# First we setup the general model, then the designs for adults and pediatrics, 
# and finally we can evaluate the separate designs and the pooled data.

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
                CL=bpop[3]*exp(b[3])*bpop[5]^a[3], # add covariate for pediatrics
                Favail=bpop[4],
                isPediatric = a[3],
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
                                  notfixed_bpop=c(1,1,1,0,1),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(0,0),
                                  m=2,
                                  groupsize=20,
                                  xt=c( 1,8,10,240,245),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=20,TAU=24,isPediatric = 0),
                                         c(DOSE=40, TAU=24,isPediatric = 0)))

# Note, to be able to use the adults FIM to combine with the pediatrics, 
# both have to have the parameter "pedCL" defined and set notfixed_bpop to 1.

## Define pediatric model/design (isPediatric = 1)
## One arm, 4 time points only
poped.db.ped <- create.poped.database(ff_fun="ff",
                                  fg_fun="sfg",
                                  fError_fun="feps",
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,pedCL=0.8), 
                                  notfixed_bpop=c(1,1,1,0,1),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(0,0),
                                  m=1,
                                  groupsize=6,
                                  xt=c( 1,2,6,240),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=40,TAU=24,isPediatric = 1)))

##  Create plot of model of adult data without variability
plot_model_prediction(poped.db, model_num_points = 300)

##  Create plot of model of pediatric (single dose level) data without variability
plot_model_prediction(poped.db.ped, model_num_points = 300)

## To store FIM from adult design - need FIM from evaluate_design
(outAdult <- evaluate_design(poped.db))
# It is obvious that we cannot estimate the pediatric covariate from adult 
# data only - therefore the message from the calculation.
# You can also note the zeros in the 4th column and 4th row of the FIM.

# We can evaluate the adult design without warning, by setting the pedCL 
# parameter to be fixed (i.e., not estimated)
evaluate_design(create.poped.database(poped.db, notfixed_bpop=c(1,1,1,0,0)))
# One obtains good estimates for all parameters for adults (<60% RSE for all).

## evaluate design of pediatrics only - insufficient
# Similarly as before with only pediatrics we cannot estimate the covariate effect, 
# so we fix it.
evaluate_design(create.poped.database(poped.db.ped, notfixed_bpop=c(1,1,1,0,0)))
# Due to having less subjects, less samples per subject, and only one dose level 
# the variability in pediatrics cannot be estimated well.

## Add adult prior
# Now we combined the two sutdies, where we assume that we only assess a difference 
# between adults and pediatrics in CL.
# We can set the prior FIM to the adult one:
poped.db.all <- create.poped.database(
  poped.db.ped,
  prior_fim = outAdult$fim
)

## evaluate design using prior FIM from adults
(out.all <- evaluate_design(poped.db.all))
# Obviously, the pooled data leads to much higher precision in parameter estimates 
# compared to the pediatrics only.

# One can also obtain the power for estimating the covariate to be different from 1.
round(evaluate_power(poped.db.all, bpop_idx=5, h0=1,out=out.all)$power,1)

