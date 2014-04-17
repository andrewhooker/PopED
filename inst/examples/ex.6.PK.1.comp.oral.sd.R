library(PopED)
library(deSolve)

PK.1.comp.oral.sd.fg.param.1 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}

bpop.vals.param.1 <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
d.vals.param.1 <- c(V=0.09,KA=0.09,CL=0.25^2)
notfixed_bpop <- c(1,1,1,0)

PK.1.comp.oral.sd.fg.param.2 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]), 
                KE=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  return( parameters ) 
}

bpop.vals.param.2 <- c(V=72.8,KA=0.25,KE=3.75/72.8,Favail=0.9)
d.vals.param.2 <- c(V=0.09,KA=0.09,KE=0.25^2)
notfixed_bpop <- c(1,1,1,0)


PK.1.comp.oral.sd.discrete.fg.param.1 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=x[1])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}


PK.1.comp.oral.sd.ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
    return(list( y= y,poped.db=poped.db))
  })
}

PK.1.comp.oral.sd.ff.ode <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
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

PK.1.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
        
    dA1  <- -KA*A1
    dA2  <- KA*A1 - KE*A2
    
    return(list(c(dA1, dA2)))
  })
}

poped.db.1 <- create.poped.database(ff_file="PK.1.comp.oral.sd.ff",
                                    fg_file="PK.1.comp.oral.sd.fg.param.1",
                                    fError_file="feps.add.prop",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,CL=0.15,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=notfixed_bpop,
                                    xt=c(0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))

poped.db.2 <- create.poped.database(ff_file="PK.1.comp.oral.sd.ff",
                                    fg_file="PK.1.comp.oral.sd.fg.param.2",
                                    fError_file="feps.add.prop",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,KE=0.15/8,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=notfixed_bpop,
                                    xt=c(0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))

poped.db.3 <- create.poped.database(ff_file="PK.1.comp.oral.sd.ff.ode",
                                    fError_file="feps.add.prop",
                                    fg_file="PK.1.comp.oral.sd.fg.param.1",
                                    groupsize=32,
                                    m=1,
                                    sigma=diag(c(0.01,0.25)),
                                    bpop=c(V=8,KA=1,CL=0.15,Favail=1), 
                                    d=c(V=0.02,KA=0.6,CL=0.07), 
                                    notfixed_bpop=notfixed_bpop,
                                    xt=c(0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=c(25,25,25,120,120,120,120,120),
                                    a=cbind(c(70)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200),
                                    mina=c(0))


for(i in 1:3){
  poped.db <- eval(parse(text=paste("poped.db.",i,sep="")))
  ##  create plot of model 
  if(i==3) print(plot_model_prediction(poped.db))
  # plot_model_prediction(poped.db,IPRED=T)
  #print(plot_model_prediction(poped.db,IPRED=T,DV=T))
  if(i!=3) print(plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T))
  
  ## evaluate initial design
  FIM <- evaluate.fim(poped.db) 
  print(FIM)
  print(det(FIM))
  print(get_rse(FIM,poped.db))
  
  # RS+SG+LS optimization of sample times
  # output <- poped_optimize(poped.db,opt_xt=T)
  # 
  # get_rse(output$fmf,output$poped.db)
  # result.db <- output$poped.db
  # plot_model_prediction(result.db,IPRED=F,DV=F)
  
  # MFEA optimization with only integer times allowed
  # mfea.output <- poped_optimize(poped.db,opt_xt=1,
  #                               bUseExchangeAlgorithm=1,
  #                               EAStepSize=1)
  # 
  # get_rse(mfea.output$fmf,mfea.output$poped.db)
  # result.db <- mfea.output$poped.db
  # plot_model_prediction(result.db,IPRED=T,DV=T)
  
  # efficiency of windows
  # plot_efficiency_of_windows(result.db,xt_windows=0.5)
  # plot_efficiency_of_windows(result.db,xt_windows=1)
  
}
