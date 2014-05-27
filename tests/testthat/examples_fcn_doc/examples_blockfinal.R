# warfarin optimization model 


FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)


blockfinal(fn="",fmf=FIM,
             dmf=dmf,
             groupsize=poped.db$groupsize,
             ni=poped.db$gni,
             xt=poped.db$gxt,
             x=poped.db$gx,a=poped.db$ga,
             model_switch=poped.db$global_model_switch,
             poped.db$param.pt.val$bpop,
             poped.db$param.pt.val$d,
             poped.db$docc,
             poped.db$param.pt.val$sigma,
             poped.db,
             opt_xt=TRUE,
             fmf_init=FIM,
             dmf_init=dmf,
             param_cvs_init=rbind(get_rse(FIM,poped.db,use_percent=FALSE)))
  
  
