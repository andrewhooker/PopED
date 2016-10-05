library(PopED)
library(deSolve)

##-- Model: One comp first order absorption, analytic solution
PK.1.comp.oral.sd.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
    return(list( y= y,poped.db=poped.db))
  })
}

# ODE solution
PK.1.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dA1  <- -KA*A1
    dA2  <- KA*A1 - KE*A2
    
    return(list(c(dA1, dA2)))
  })
}

PK.1.comp.oral.sd.ff.ode <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE, A2 = 0)
    times <- drop(xt)
    times <- sort(times) 
    times <- c(0,times) # add extra time for start of integration
    out   <- ode(A_ini, times, PK.1.comp.oral.ode, parameters) #,atol=1e-13,rtol=1e-13)
    y = out[,"A2"]/(V/Favail)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) ## must be a row vector
    return(list( y= y,poped.db=poped.db)) 
  })
}

## -- parameter definition function 
PK.1.comp.oral.sd.fg.param.1 <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}

PK.1.comp.oral.sd.fg.param.2 <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]), 
                KE=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  return( parameters ) 
}

poped.db.1 <- create.poped.database(ff_fun="PK.1.comp.oral.sd.ff",
                                    fg_fun="PK.1.comp.oral.sd.fg.param.1",
                                    fError_fun="feps.add.prop",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,CL=0.15,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=c(1,1,1,0),
                                    xt=c(1,2,3,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    discrete_xt = list(1:120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))

poped.db.2 <- create.poped.database(poped.db.1,
                                    fg_fun="PK.1.comp.oral.sd.fg.param.2",
                                    bpop=c(V=8,KA=1,KE=0.15/8,Favail=1),
                                    d=c(V=0.02,KA=0.6,KE=0.07)) 
                                    

poped.db.3 <- create.poped.database(poped.db.1,
                                    ff_fun="PK.1.comp.oral.sd.ff.ode")
                                    

##  create plot of models
(plot1 <- plot_model_prediction(poped.db.1,IPRED=T))
(plot2 <- plot_model_prediction(poped.db.2,IPRED=T))
(plot3 <- plot_model_prediction(poped.db.3,IPRED=T))

# different results for different parameterizations
library(ggplot2)
library(gridExtra)
plot1 <- plot1 + ggtitle("CL Parameterization")
plot2 <- plot2 + ggtitle("KE Parameterization")
grid.arrange(plot1,plot2)


## evaluate initial designs
# different results for different parameterizations
evaluate_design(poped.db.1)
evaluate_design(poped.db.2)
evaluate_design(poped.db.3)
  

# Optimization of sample times
# different results for different parameterizations
output.1 <- poped_optim(poped.db.1,opt_xt=T,parallel=T)
output.2 <- poped_optim(poped.db.2,opt_xt=T,parallel=T)

summary(output.1)
summary(output.2)
