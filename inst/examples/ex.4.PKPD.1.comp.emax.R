library(PopED)

ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    MS <- model_switch
    
    # PK model
    CONC = DOSE/V*exp(-CL/V*xt) 
    
    # PD model
    EFF = E0 + CONC*EMAX/(EC50 + CONC)
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

## -- parameter definition function
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( 
    CL=bpop[1]*exp(b[1])  ,
    V=bpop[2]*exp(b[2])	,
    E0=bpop[3]*exp(b[3])	,
    EMAX=bpop[4],
    EC50=bpop[5]*exp(b[4])	,
    DOSE=a[1]
  )
  return( parameters ) 
}

## -- Residual Error function
## -- Proportional PK + additive PD
feps <- function(model_switch,xt,parameters,epsi,poped.db){
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  MS <- model_switch
  
  prop.err <- y*(1+epsi[,1])
  add.err <- y+epsi[,2]
  
  y[MS==1] = prop.err[MS==1]
  y[MS==2] = add.err[MS==2]
  
  return(list( y= y,poped.db =poped.db )) 
}

poped.db <- create.poped.database(ff_fun=ff,
                                  fError_fun=feps,
                                  fg_fun=sfg,
                                  groupsize=20,
                                  m=3,
                                  sigma=diag(c(0.15,0.015)),
                                  bpop=c(CL=0.5,V=0.2,E0=1,EMAX=1,EC50=1),  
                                  d=c(CL=0.09,V=0.09,E0=0.04,EC50=0.09), 
                                  xt=c( 0.33,0.66,0.9,5,0.1,1,2,5),
                                  bUseGrouped_xt=1,
                                  model_switch=c( 1,1,1,1,2,2,2,2),
                                  minxt=0,
                                  maxxt=5,
                                  ourzero = 0,
                                  a=list(c(DOSE=0),c(DOSE=1),c(DOSE=2)),
                                  maxa=c(DOSE=10),
                                  mina=c(DOSE=0))


plot_model_prediction(poped.db,facet_scales="free")
plot_model_prediction(poped.db,IPRED=T,DV=T,facet_scales="free",separate.groups=T)

## evaluate initial design
evaluate_design(poped.db)
shrinkage(poped.db)

# Optimization of sample times and doses
output <- poped_optim(poped.db, opt_xt = T, opt_a = T, parallel = T,method = c("LS"))

get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db,facet_scales="free")


