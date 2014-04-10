# load the model and database
source("MODEL.1.PK.one.comp.oral.SOLUTION.R")

## Choose database
poped.db <- poped.db.1 # original model and design
poped.db <- poped.db.2 # updated with correct residual error
poped.db <- poped.db.3 # parameters fixed

##  create plot of model 
plot_model_prediction(poped.db)
plot_model_prediction(poped.db,IPRED=T)
plot_model_prediction(poped.db,IPRED=T,DV=T)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of sample times
output <- poped_optimize(poped.db,opt_xt=T)

get_rse(output$fmf,output$poped.db)
result.db <- output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)

# Dose and time optimization
output <- poped_optimize(poped.db.4,opt_xt=T,opt_a=T)

get_rse(output$fmf,output$poped.db)
result.db <- output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)
