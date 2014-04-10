# load the model and database
source("MODEL.1.PK.one.comp.oral.SOLUTION.R")

## Choose database
poped.db <- poped.db.5 # Design space reduced according to clinical constraints

##  create plot of model 
plot_model_prediction(poped.db)
plot_model_prediction(poped.db,IPRED=T)
plot_model_prediction(poped.db,IPRED=T,DV=T)
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of sample times
output <- poped_optimize(poped.db,opt_xt=T)

get_rse(output$fmf,output$poped.db)
result.db <- output$poped.db
plot_model_prediction(result.db,IPRED=F,DV=F)

# MFEA optimization with only integer times allowed
mfea.output <- poped_optimize(poped.db,opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=1)

get_rse(mfea.output$fmf,mfea.output$poped.db)
result.db <- mfea.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)

# efficiency of windows
plot_efficiency_of_windows(result.db,xt_windows=0.5)
plot_efficiency_of_windows(result.db,xt_windows=1)


