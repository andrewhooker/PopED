
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
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
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

## -- Define design and design space
poped.db <- create.poped.database(ff_fun=ff,
                                  fg_fun=sfg,
                                  fError_fun=feps,
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(0,0),
                                  m=2,
                                  groupsize=20,
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                  maxa=c(DOSE=200,TAU=24),
                                  mina=c(DOSE=0,TAU=24))

##  create plot of model without variability 
plot_model_prediction(poped.db, model_num_points = 300)


\dontrun{
  
  ##  create plot of model with variability 
  plot_model_prediction(poped.db, IPRED=T, DV=T, separate.groups=T, model_num_points = 300)
  
}

## evaluate initial design
evaluate_design(poped.db)

\dontrun{
  
  # Optimization of sample times
  output <- poped_optim(poped.db, opt_xt=TRUE, parallel = TRUE)
  summary(output)
  get_rse(output$FIM, output$poped.db)
  plot_model_prediction(output$poped.db)
  
  # Optimization of sample times and doses
  output_2 <- poped_optim(poped.db, opt_xt=TRUE, opt_a=TRUE, parallel = TRUE)
  summary(output_2)
  get_rse(output_2$FIM,output_2$poped.db)
  plot_model_prediction(output_2$poped.db)
  
  # Optimization of sample times with only integer time points in design space
  # faster than continuous optimization in this case
  poped.db.discrete <- create.poped.database(ff_fun=ff,
                                             fg_fun=sfg,
                                             fError_fun=feps,
                                             bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                             notfixed_bpop=c(1,1,1,0),
                                             d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                             sigma=c(0.04,5e-6),
                                             notfixed_sigma=c(0,0),
                                             m=2,
                                             groupsize=20,
                                             xt=c( 1,2,8,240,245),
                                             minxt=c(0,0,0,240,240),
                                             maxxt=c(10,10,10,248,248),
                                             discrete_xt = list(0:248),
                                             bUseGrouped_xt=1,
                                             a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                             maxa=c(DOSE=200,TAU=24),
                                             mina=c(DOSE=0,TAU=24),
                                             ourzero = 0)
  
  output_discrete <- poped_optim(poped.db.discrete, opt_xt=T, parallel = TRUE)
  
  
  summary(output_discrete)
  get_rse(output_discrete$FIM,output_discrete$poped.db)
  plot_model_prediction(output_discrete$poped.db)
  
  # Efficiency of sampling windows
  plot_efficiency_of_windows(output_discrete$poped.db,xt_windows=0.5)
  plot_efficiency_of_windows(output_discrete$poped.db,xt_windows=1)
  
}
