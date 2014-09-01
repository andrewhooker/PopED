
#' An implementation of a two compartment model with oral absorption using ODEs.
library(PopED)
library(deSolve)

#' Define the ODE system
PK.2.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {    
    dA1 <- -KA*A1 
    dA2 <- KA*A1 + A3* Q/V2 -A2*(CL/V1+Q/V1)
    dA3 <- A2* Q/V1-A3* Q/V2
    return(list(c(dA1, dA2, dA3)))
  })
}

#' define the initial conditions and the dosing
ff.PK.2.comp.oral.md.ode <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0, A3=0)
    times_xt <- drop(xt)
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, PK.2.comp.oral.ode, parameters, events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V1/Favail)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

#' parameter definition function 
#' names match parameters in function ff
fg <- function(x,a,bpop,b,bocc){
  parameters=c( CL=bpop[1]*exp(b[1]),
                V1=bpop[2],
                KA=bpop[3]*exp(b[2]),
                Q=bpop[4],
                V2=bpop[5],
                Favail=bpop[6],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}


#' create poped database
poped.db <- create.poped.database(ff_file="ff.PK.2.comp.oral.md.ode",
                                  fError_file="feps.add.prop",
                                  fg_file="fg",
                                  groupsize=20,
                                  m=1,      #number of groups
                                  sigma=c(prop=0.1^2,add=0.05^2),
                                  bpop=c(CL=10,V1=100,KA=1,Q= 3.0, V2= 40.0, Favail=1),
                                  d=c(CL=0.15^2,KA=0.25^2),
                                  notfixed_bpop=c(1,1,1,1,1,0),
                                  xt=c( 48,50,55,65,70,85,90,120), 
                                  minxt=0,
                                  maxxt=144,
                                  a=c(100,24),
                                  maxa=c(1000,24),
                                  mina=c(0,8))

#' plot intial design just PRED
plot_model_prediction(poped.db)

#' plot intial design with BSV and RUV in model
plot_model_prediction(poped.db,IPRED=T,DV=T)

#' how long does one evaluation of the FIM take? 
tic()
FIM <- evaluate.fim(poped.db) 
toc()

#' and the results
FIM
det(FIM)
get_rse(FIM,poped.db)

#' optimize using line search using the following code:
# ls.output <- poped_optimize(poped.db,opt_xt=1,
#                             bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 1,
#                             ls_step_size=1)

#' Plot final design
# plot_model_prediction(ls.output$poped.db)

