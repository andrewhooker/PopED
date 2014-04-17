sfg.emax.hill <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( EMAX=bpop[1]*exp(b[1]),
                ED50=bpop[2]*exp(b[2]),
                GAMMA=bpop[3],
                BASE=bpop[4]+b[3])
  return( parameters ) 
}


ff.emax.hill <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    DOSE = xt
    y=BASE + EMAX*DOSE^(GAMMA)/(ED50^(GAMMA) + DOSE^(GAMMA))
    return(list( y= y,poped.db=poped.db))
  })
}

poped.db <- create.poped.database(ff_file="ff.emax.hill",
                                  fError_file="feps.add.prop",
                                  fg_file="sfg.emax.hill",
                                  groupsize=100,
                                  m=1,
                                  bpop=c(EMAX=100,ED50=20,GAMMA=4.5,BASE=1),  
                                  d=c(EMAX=0.0625,ED50=0.0625,BASE=0.0625), 
                                  sigma=diag(c(0.01,.1)),
                                  xt=seq(0,50,length.out=8),
                                  minxt=0,
                                  maxxt=50,
                                  ourzero=0)

plot_model_prediction(poped.db,IPRED=T,DV=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of doses 
output <- poped_optimize(poped.db,opt_xt=T)

# plot output of optimization
plot_model_prediction(output$poped.db,IPRED=T,DV=T)









