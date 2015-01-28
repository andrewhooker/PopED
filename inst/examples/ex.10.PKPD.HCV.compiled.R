## HCV example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)
library(deSolve)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c(p=bpop[1],
               d=bpop[2],
               e=bpop[3],
               s=bpop[4],
               KA=bpop[5] + b[1],
               KE=bpop[6] + b[2],
               VD=bpop[7] + b[3],
               EC50=bpop[8] + b[4],
               n=bpop[9] + b[5],
               delta=bpop[10] + b[6],
               c=bpop[11] + b[7],
               DOSE=a[1],
               TINF=a[2],
               TAU=a[3])
  return(parameters)
}

#' To make evaluation time more reasonable we use compiled code
#' To set this up see the 
#' "R Package deSolve, Writing Code in Compiled Languages" 
#' vingette in the deSolve documentation
#' make sure your working directory is where this file is located
system("R CMD SHLIB HCV_ode.c")
dyn.load(paste("HCV_ode", .Platform$dynlib.ext, sep = ""))

ff.ODE.compiled <- function(model_switch,xt,parameters,poped.db){
  parameters[5:11] <- exp(parameters[5:11])
  with(as.list(parameters),{
    A_ini  <- c(A1 = 0, A2 = 0, A3=c*delta/(p*e), 
                A4=(s*e*p-d*c*delta)/(p*delta*e),
                A5=(s*e*p-d*c*delta)/(c*delta*e))
    
    #Set up time points for the ODE
    times_xt <- drop(xt)
    times <- c(0,times_xt) ## add extra time for start of integration
    times <- sort(times) 
    times <- unique(times) # remove duplicates
    
    # compute values from ODEs
    #out   <- ode(A_ini, times, model.ode, parameters,atol=1e-10,rtol=1e-10)
    out <- ode(A_ini, times, func = "derivs", parms = parameters,
               #jacfunc = "jac", # not really needed, speed up is minimal if this is defined or not.
               dllname = "HCV_ode",
               initfunc = "initmod", #nout = 1, outnames = "Sum",
               atol=1e-10,rtol=1e-10)
    
    # grab timepoint values
    out = out[match(times_xt,out[,"time"]),]
    
    y <- xt*0
    pk <- out[,"A2"]/VD
    pd <- out[,"A5"]
    y[model_switch==1] <- pk[model_switch==1]
    y[model_switch==2] <- log10(pd[model_switch==2])
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}

feps.ODE.compiled <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  MS<-model_switch
  y <- ff.ODE.compiled(model_switch,xt,parameters,poped.db)[[1]]
  
  y[MS==1] <- y[MS==1]+epsi[,1]
  y[MS==2] <- y[MS==2]+epsi[,2]
  
  return(list(y=y,poped.db=poped.db))
}

## -- Define initial design  and design space
poped.db.compiled <- create.poped.database(ff_file="ff.ODE.compiled",
                                           fg_file="sfg",
                                           fError_file="feps.ODE.compiled",
                                           bpop=c(p=100,
                                                  d=0.001,
                                                  e=1e-7,
                                                  s=20000,
                                                  KA=log(0.8),
                                                  KE=log(0.15),
                                                  VD=log(100), 
                                                  #VD=100000,
                                                  EC50=log(0.12), 
                                                  #EC50=0.00012,
                                                  n=log(2),
                                                  delta=log(0.2),
                                                  c=log(7)),
                                           notfixed_bpop=c(0,0,0,0,1,1,1,1,1,1,1),
                                           d=c(KA=0.25,
                                               KE=0.25,
                                               VD=0.25,
                                               EC50=0.25,
                                               n=0.25,
                                               delta=0.25,
                                               c=0.25),
                                           sigma=c(0.04,0.04),
                                           #sigma=c(4e-8,0.04),
                                           groupsize=30,
                                           xt=c(0.00001,0.25,0.5,1,2,3,4,7,10,14,21,28,
                                                0,0.25,0.5,1,2,3,4,7,10,14,21,28),
                                           model_switch=c(rep(1,12),rep(2,12)),
                                           a=c(180,1,7))

##  create plot of model without variability 
plot_model_prediction(poped.db.compiled, facet_scales = "free")


#' how long does one evaluation of the FIM take? 
tic()
FIM.compiled <- evaluate.fim(poped.db.compiled) 
toc()

#' design evaluation compared to 
#' Nyberg et al., "Methods and software tools for design evaluation 
#' for population pharmacokinetics-pharmacodynamics studies", 
#' Br. J. Clin. Pharm., 2014. 
crit_reference <- 248.8
rse_reference <- c(12.1,10.5,10.0,15.8,10.4,9.4,11.0,40.0,30.8,28.8,60.4,28.8,27.2,32.8,8.5,9.0)

crit <- det(FIM.compiled)^(1/length(get_unfixed_params(poped.db)[["all"]]))
crit

# difference between crit and reference
round(crit,1) - crit_reference 

rse <- get_rse(FIM.compiled,poped.db.compiled) # this is for the log of the fixed effect parameters
rse_norm <- sqrt(diag(inv(FIM.compiled)))*100 # this is approximately the RSE for the normal scale of the fixed effects 
rse[1:7] <- rse_norm[1:7]
rse

# difference between %rse and reference %rse
round(rse,1) - rse_reference

