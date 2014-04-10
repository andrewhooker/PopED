# load the model and database
source("MODEL.2.PK.one.comp.oral.PD.imax.SOLUTION.R")

## Choose database
poped.db <- poped.db.4 # With uncertainty around IC50

## ED evaluate using API. 
output <- evaluate.e.ofv.fim(poped.db,ofv_calc_type=4)
output$E_ofv

## Increase sample size
output <- evaluate.e.ofv.fim(poped.db,ofv_calc_type=4,ED_samp_size=100)
output$E_ofv

# MFEA for API
mfea.output <- poped_optimize(poped.db,
                              opt_xt=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=1,
                              ofv_calc_type=4,
                              d_switch=0)

get_rse(mfea.output$fmf,mfea.output$poped.db)
result.db <- mfea.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T,separate.groups=T,facet_scales="free")



