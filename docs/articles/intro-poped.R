## ----eval=TRUE-----------------------------------------------------------
library(PopED)

## ----include = FALSE-----------------------------------------------------
set.seed(1234)
knitr::opts_chunk$set(cache = FALSE)

## ----struct_model--------------------------------------------------------
ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    N = floor(xt/TAU)+1
    y=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
      (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    return(list( y=y,poped.db=poped.db))
  })
}

## ------------------------------------------------------------------------
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
}

## ------------------------------------------------------------------------
feps <- function(model_switch,xt,parameters,epsi,poped.db){
  returnArgs <- ff(model_switch,xt,parameters,poped.db) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
 
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list(y=y,poped.db=poped.db)) 
}

## ------------------------------------------------------------------------
poped.db <- create.poped.database(ff_fun=ff,
                                  fg_fun=sfg,
                                  fError_fun=feps,
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(1,0),
                                  m=2,
                                  groupsize=20,
                                  a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                  maxa=c(DOSE=200,TAU=24),
                                  mina=c(DOSE=0,TAU=24),
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  bUseGrouped_xt=TRUE)

## ---- fig.width=6--------------------------------------------------------
plot_model_prediction(poped.db, model_num_points = 500)

## ----simulate_with_BSV, fig.width=6--------------------------------------
plot_model_prediction(poped.db, model_num_points=500, IPRED=T)

## ------------------------------------------------------------------------
dat <- model_prediction(poped.db,DV=T)
head(dat,n=5);tail(dat,n=5)

## ------------------------------------------------------------------------
evaluate_design(poped.db)

## ---- fig.width=6--------------------------------------------------------
summary(output)
plot_model_prediction(output$poped.db)

## ---- message = FALSE,results='hide'-------------------------------------
poped.db.discrete <- create.poped.database(poped.db,discrete_xt = list(0:248))
                                          
output_discrete <- poped_optim(poped.db.discrete, opt_xt=T)


## ----fig.width=6---------------------------------------------------------
summary(output_discrete)
plot_model_prediction(output_discrete$poped.db)

## ------------------------------------------------------------------------
crit_fcn <- function(poped.db,...){
  pred_df <- model_prediction(poped.db)
  sum((pred_df[pred_df["Time"]==240,"PRED"] - c(0.2,0.35))^2)
}
crit_fcn(output$poped.db)

## ---- fig.width=6--------------------------------------------------------
summary(output_cost)
get_rse(output_cost$FIM, output_cost$poped.db)
plot_model_prediction(output_cost$poped.db)

