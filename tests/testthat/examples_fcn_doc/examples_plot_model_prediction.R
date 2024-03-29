## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.md.CL

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
    return(parameters) 
}

## -- Define initial design  and design space
poped.db <- create.poped.database(
  ff_fun=ff.PK.1.comp.oral.sd.CL,
  fg_fun=sfg,
  fError_fun=feps.prop,
  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
  notfixed_bpop=c(1,1,1,0),
  d=c(CL=0.07, V=0.02, KA=0.6), 
  sigma=0.01,
  groupsize=32,
  xt=c( 0.5,1,2,6,24,36,72,120),
  minxt=0,
  maxxt=120,
  a=70)

##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability by simulating from OMEGA and SIGMA
plot_model_prediction(poped.db,IPRED=TRUE,DV=TRUE)

##  create plot of model with variability by 
##  computing the expected variance (using an FO approximation) 
##  and then computing a prediction interval 
##  based on an assumption of normality
##  computation is faster but less accurate 
##  compared to using DV=TRUE (and groupsize_sim = 500)
plot_model_prediction(poped.db,PI=TRUE)

##-- Model: One comp first order absorption + inhibitory imax
## -- works for both mutiple and single dosing  
ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    
    y=xt
    MS <- model_switch
    
    # PK model
    N = floor(xt/TAU)+1
    CONC=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
      (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    
    # PD model
    EFF = E0*(1 - CONC*IMAX/(IC50 + CONC))
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

## -- parameter definition function 
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU = a[2],
                E0=bpop[5]*exp(b[4]),
                IMAX=bpop[6],
                IC50=bpop[7])
  return( parameters ) 
}


## -- Residual Error function
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

poped.db <- 
  create.poped.database(
    ff_fun=ff,
    fError_fun=feps,
    fg_fun=sfg,
    groupsize=20,
    m=3,
    bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,
           E0=1120,IMAX=0.807,IC50=0.0993),  
    notfixed_bpop=c(1,1,1,0,1,1,1),
    d=c(V=0.09,KA=0.09,CL=0.25^2,E0=0.09), 
    sigma=c(0.04,5e-6,0.09,100),
    notfixed_sigma=c(0,0,0,0),
    xt=c( 1,2,8,240,240,1,2,8,240,240),
    minxt=c(0,0,0,240,240,0,0,0,240,240),
    maxxt=c(10,10,10,248,248,10,10,10,248,248),
    discrete_xt = list(0:248),
    G_xt=c(1,2,3,4,5,1,2,3,4,5),
    bUseGrouped_xt=1,
    model_switch=c(1,1,1,1,1,2,2,2,2,2),
    a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24),c(DOSE=0, TAU=24)),
    maxa=c(DOSE=200,TAU=40),
    mina=c(DOSE=0,TAU=2),
    ourzero=0)

##  create plot of model and design 
plot_model_prediction(poped.db,facet_scales="free",
                      model.names = c("PK","PD"))

##  create plot of model with variability by  
##  computing the expected variance (using an FO approximation) 
##  and then computing a prediction interval 
##  based on an assumption of normality
##  computation is faster but less accurate 
##  compared to using DV=TRUE (and groupsize_sim = 500)
plot_model_prediction(poped.db,facet_scales="free",
                      model.names = c("PK","PD"),
                      PI=TRUE,
                      separate.groups = TRUE)


