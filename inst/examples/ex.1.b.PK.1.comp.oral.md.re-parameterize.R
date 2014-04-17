library(PopED)

## using libary models and reparameterizing the problen to KA, KE and V 

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
                                    a=cbind(c(20,40),c(24,24)),
                                    bUseGrouped_xt=1,
                                    maxa=c(200,40),
                                    mina=c(0,2))


##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability 
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of sample times
output <- poped_optimize(poped.db,opt_xt=T)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db,IPRED=F,DV=F)

# RS+SG+LS optimization of sample times and doses
output <- poped_optimize(poped.db,opt_xt=T,opt_a=T)
get_rse(output$fmf,output$poped.db)
plot_model_prediction(output$poped.db,IPRED=F,DV=F)

# MFEA optimization with only integer times allowed
mfea.output <- poped_optimize(poped.db,opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=1)
get_rse(mfea.output$fmf,mfea.output$poped.db)
plot_model_prediction(mfea.output$poped.db)

# Efficiency of sampling windows
plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5)
plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1)

