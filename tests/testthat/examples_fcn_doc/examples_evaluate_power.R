# Folowing the examples presented in Retout, 2007

ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    
    lambda1 <- lam1a
    if(TREAT==2) lambda1 <- lam1b
    
    y=log10(P1*exp(-lambda1*xt)+P2*exp(-lam2*xt))
    
    return(list(y=y,poped.db=poped.db))
  })
}

sfg <- function(x,a,bpop,b,bocc){
  parameters=c(P1=exp(bpop[1]+b[1]),
               P2=exp(bpop[2]+b[2]),
               lam1a=exp(bpop[3]+b[3]),
               lam1b=exp(bpop[3]+bpop[4]+b[3]),
               lam2=exp(bpop[5]+b[4]),
               TREAT=a[1])
  return(parameters) 
}
  

poped.db <- create.poped.database(ff_fun = ff,
                                  fg_fun = sfg,
                                  fError_fun = feps.add,
                                  bpop=c(P1=12, P2=8,
                                         lam1=-0.7,beta=0,lam2=-3.0),
                                  d=c(P1=0.3, P2=0.3,
                                      lam1=0.3,lam2=0.3), 
                                  sigma=c(0.065^2),
                                  groupsize=100,
                                  m=2,
                                  xt=c(1, 3, 7, 14, 28, 56),
                                  minxt=0,
                                  maxxt=100,
                                  a=list(c(TREAT=1),c(TREAT=2)))

plot_model_prediction(poped.db)
evaluate_design(poped.db)

poped.db_2 <- create.poped.database(poped.db,bpop=c(P1=12, P2=8,
                                      lam1=-0.7,beta=0.262,lam2=-3.0))

plot_model_prediction(poped.db_2)
evaluate_design(poped.db_2)

evaluate_power(poped.db_2,bpop_idx = 4)
