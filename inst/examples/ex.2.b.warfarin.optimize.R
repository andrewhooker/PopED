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
plot_model_prediction(poped.db,IPRED=T,DV=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

##############
# Optimization
##############

# below are a number of ways to optimize the problem

# RS+SG+LS optimization of sample times
# optimization with just a few iterations
output <- poped_optimize(poped.db,opt_xt=T,
                         rsit=5,sgit=5,ls_step_size=5)

# RS+SG+LS optimization of sample times (longer optimization times)
output <- poped_optimize(poped.db,opt_xt=T)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db)

# MFEA optimization with only integer times allowed
mfea.output <- poped_optimize(poped.db,opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=1)
get_rse(mfea.output$fmf,mfea.output$poped.db)
plot_model_prediction(mfea.output$poped.db)

# Examine efficiency of sampling windows
plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5)
plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1)

# Random search optimization (300 itertations) of DOSE and sampling times
rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=300)
rs.output$xt
names(rs.output)
result.db <- rs.output$poped.db
plot_model_prediction(result.db)

# RS using the full FIM
rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,fim.calc.type=0) 

# RS within poped_optimize (just a few samples here)
rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                            bUseRandomSearch= 1,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0)
names(rs.output)
get_rse(rs.output$fmf,rs.output$poped.db)

# line search, DOSE and sample time optimization
ls.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 1,
                            ls_step_size=10)

# Stochastic gradient search, DOSE and sample time optimization
sg.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1, 
                            bUseRandomSearch= 0,bUseStochasticGradient = 1,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            sgit=20)

# BFGS search, DOSE and sample time optimization
bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                              bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 1,bUseLineSearch = 0)

