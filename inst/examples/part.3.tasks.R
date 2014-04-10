# load the model and database
source("MODEL.2.PK.one.comp.oral.PD.imax.SOLUTION.R")

## Choose database
poped.db <- poped.db.1 # original model and design
poped.db <- poped.db.2  # reducing design space
poped.db <- poped.db.3 # Add placebo group

##  create plot of model 
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T,facet_scales="free")

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# MFEA optimization with only integer times allowed
mfea.output <- poped_optimize(poped.db,opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=1)

get_rse(mfea.output$fmf,mfea.output$poped.db)
result.db <- mfea.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T,separate.groups=T,facet_scales="free")



