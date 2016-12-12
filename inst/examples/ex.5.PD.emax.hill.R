library(PopED)


ff.emax.hill <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    DOSE = xt
    y=BASE + EMAX*DOSE^(GAMMA)/(ED50^(GAMMA) + DOSE^(GAMMA))
    return(list( y= y,poped.db=poped.db))
  })
}

## -- parameter definition function
sfg.emax.hill <- function(x,a,bpop,b,bocc){
  parameters=c( EMAX=bpop[1]*exp(b[1]),
                ED50=bpop[2]*exp(b[2]),
                GAMMA=bpop[3],
                BASE=bpop[4]+b[3])
  return( parameters ) 
}

poped.db <- create.poped.database(ff_fun=ff.emax.hill,
                                  fError_fun=feps.add.prop,
                                  fg_fun=sfg.emax.hill,
                                  groupsize=100,
                                  m=1,
                                  bpop=c(EMAX=100,ED50=20,GAMMA=4.5,BASE=1),  
                                  d=c(EMAX=0.0625,ED50=0.0625,BASE=0.0625), 
                                  sigma=diag(c(0.01,.1)),
                                  xt=seq(0,50,length.out=8),
                                  minxt=0,
                                  maxxt=50,
                                  ourzero=0)

library(ggplot2)
plot1 <- plot_model_prediction(poped.db,IPRED=T,DV=T)
plot1 + xlab("Dose")

## evaluate initial design
evaluate_design(poped.db)

# Optimization of doses
output <- poped_optim(poped.db, opt_xt = T, parallel = T)

summary(output)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db)









