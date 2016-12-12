## using libary models and reparameterizing the problen to KA, KE and V 
## optimization of dose and dose interval
library(PopED)


## -- names match parameters in function defined in ff_file
fg.PK.1.comp.oral.md.param.2 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]), 
                KE=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}

## -- Define design and design space
poped.db <- create.poped.database(ff_file="ff.PK.1.comp.oral.md.KE",
                                  fg_file="fg.PK.1.comp.oral.md.param.2",
                                  fError_file="feps.add.prop",
                                  groupsize=20,
                                  m=2,
                                  sigma=c(0.04,5e-6),
                                  bpop=c(V=72.8,KA=0.25,KE=3.75/72.8,Favail=0.9), 
                                  d=c(V=0.09,KA=0.09,KE=0.25^2), 
                                  notfixed_bpop=c(1,1,1,0),
                                  notfixed_sigma=c(0,0),
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                  maxa=c(DOSE=200,TAU=40),
                                  mina=c(DOSE=0,TAU=2))


##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability 
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)

## evaluate initial design
evaluate_design(poped.db)

# Optimization of sample times
output <- poped_optim(poped.db, opt_xt =TRUE, parallel=TRUE)

# Evaluate optimization results
summary(output)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db)

# Optimization of sample times, doses and dose intervals
output_2 <- poped_optim(output$poped.db, opt_xt =TRUE, opt_a = TRUE, parallel = TRUE)

summary(output_2)
get_rse(output_2$FIM,output_2$poped.db)
plot_model_prediction(output_2$poped.db)

# Optimization of sample times with only integer time points in design space
# faster than continuous optimization in this case
poped.db.discrete <- create.poped.database(poped.db,discrete_xt = list(0:248))

output_discrete <- poped_optim(poped.db.discrete, opt_xt=T, parallel = TRUE)

summary(output_discrete)
get_rse(output_discrete$FIM,output_discrete$poped.db)
plot_model_prediction(output_discrete$poped.db)


# Efficiency of sampling windows
plot_efficiency_of_windows(output_discrete$poped.db, xt_windows=1)

