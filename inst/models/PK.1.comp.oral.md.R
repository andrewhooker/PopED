library(PopED)
library(deSolve)

PK.1.comp.oral.md.fg.param.1 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}

bpop.vals.param.1 <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
d.vals.param.1 <- c(V=0.09,KA=0.09,CL=0.25^2)
notfixed_bpop <- c(1,1,1,0)

PK.1.comp.oral.md.fg.param.2 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]), 
                KE=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}

bpop.vals.param.2 <- c(V=72.8,KA=0.25,KE=3.75/72.8,Favail=0.9)
d.vals.param.2 <- c(V=0.09,KA=0.09,KE=0.25^2)
notfixed_bpop <- c(1,1,1,0)


PK.1.comp.oral.md.ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  with(as.list(parameters),{
    y=xt
    N = floor(xt/TAU)+1
    y=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * TAU)) * (1 - exp(-N * KE * TAU))/(1 - exp(-KE * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    return(list( y=y,poped.db=poped.db))
  })
}

PK.1.comp.oral.md.ff.ode <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt) #xt[,,drop=T] 
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, PK.1.comp.oral.ode, parameters, events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V/Favail)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

PK.1.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {    
    dA1 <- -KA*A1
    dA2 <- KA*A1 - KE*A2
    return(list(c(dA1, dA2)))
  })
}


poped.db.1 <- create.poped.database(ff_file="PK.1.comp.oral.md.ff",
                                    fg_file="PK.1.comp.oral.md.fg.param.1",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    m=2,
                                    sigma=diag(c(0.04,5e-6)),
                                    bpop=bpop.vals.param.1, 
                                    d=d.vals.param.1, 
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=c(0,0),
                                    xt=c( 1,2,8,240,245),
                                    minxt=c(0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248),
                                    a=cbind(c(20,40),c(24,24)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200,40),
                                    mina=c(0,2))

poped.db.2 <- create.poped.database(ff_file="PK.1.comp.oral.md.ff",
                                    fg_file="PK.1.comp.oral.md.fg.param.2",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    m=2,
                                    sigma=diag(c(0.04,5e-6)),
                                    bpop=bpop.vals.param.2, 
                                    d=d.vals.param.2, 
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=c(0,0),
                                    xt=c( 1,2,8,240,245),
                                    minxt=c(0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248),
                                    a=cbind(c(20,40),c(24,24)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200,40),
                                    mina=c(0,2))

poped.db.3 <- create.poped.database(ff_file="PK.1.comp.oral.md.ff.ode",
                                    fError_file="feps.add.prop",
                                    fg_file="PK.1.comp.oral.md.fg.param.1",
                                    groupsize=20,
                                    m=2,      #number of groups
                                    sigma=diag(c(0.04,5e-6)),
                                    bpop=bpop.vals.param.1,  
                                    d=d.vals.param.1, 
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=c(0,0),
                                    xt=c( 1,2,8,240,245),
                                    minxt=c(0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248),
                                    a=cbind(c(20,40),c(24,24)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200,40),
                                    mina=c(0,2))


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