## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error to 
##   avoid sample times at very low concentrations (time 0 or very late samoples).
library(PopED)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list(y=y,poped.db=poped.db))
  })
}

feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional + additive
  y <- ff(model_switch,xt,parameters,poped.db)[[1]] 
  y = y*(1+epsi[,1]) + epsi[,2]  
  return(list(y=y,poped.db=poped.db)) 
}


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff",
                                  fg_file="sfg",
                                  fError_file="feps",
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70,
                                  mina=0,
                                  maxa=100)

##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability 
plot_model_prediction(poped.db,PI=TRUE)

## evaluate initial design
evaluate_design(poped.db)

##############
# Optimization
##############

# below are a number of ways to optimize the problem

# Optimization of sample times
output <- poped_optim(poped.db, opt_xt = T, parallel = T)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db)


# Examine efficiency of sampling windows
plot_efficiency_of_windows(output$poped.db,xt_windows=0.5)


# Optimization of DOSE and sampling times
output_D_T <- poped_optim(poped.db, opt_xt = T, opt_a = T, parallel = T)

summary(output_D_T)
plot_model_prediction(output_D_T$poped.db)

# Discrete optimization with only integer times allowed
# and Dose in units of 10
poped.db.discrete <- create.poped.database(poped.db,
                                           xt=c( 1,2,3,6,24,36,72,120),
                                           discrete_xt = list(0:120),
                                           discrete_a = list(seq(10,100,10)))

output_discrete <- poped_optim(poped.db.discrete, opt_xt = T, opt_a = T, parallel = T)

get_rse(output_discrete$FIM,output_discrete$poped.db)
plot_model_prediction(output_discrete$poped.db)


# Optimization using a genetic algorithm
output_ga <- poped_optim(poped.db, opt_xt = T, parallel = T, method = c("GA"))

summary(output_ga)
plot_model_prediction(output_ga$poped.db)
