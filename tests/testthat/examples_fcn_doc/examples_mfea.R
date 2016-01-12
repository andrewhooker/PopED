##############
# typically one will use poped_optimize 
# This then calls mfea 
##############

# optimization of covariate, with coarse grid
out_1 <- poped_optimize(poped.db,opt_a=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=25)


\dontrun{
  
  
  
  # MFEA optimization with only integer times allowed
  out_2 <- poped_optimize(poped.db,opt_xt=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(out_2$fmf,out_2$poped.db)
  plot_model_prediction(out_2$poped.db)
  
  
  ##############
  # If you really want to you can use mfea dirtectly
  ##############
  dsl <- downsizing_general_design(poped.db)
  
  output <- mfea(poped.db,
                 model_switch=dsl$model_switch,
                 ni=dsl$ni,
                 xt=dsl$xt,
                 x=dsl$x,
                 a=dsl$a,
                 bpopdescr=dsl$bpop,
                 ddescr=dsl$d,
                 maxxt=dsl$maxxt,
                 minxt=dsl$minxt,
                 maxa=dsl$maxa,
                 mina=dsl$mina,
                 fmf=0,dmf=0,
                 EAStepSize=1,
                 opt_xt=1)
  
  
}


