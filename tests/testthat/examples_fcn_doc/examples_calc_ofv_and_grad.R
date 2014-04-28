# warfarin optimization model 

opta=TRUE
aa=opta*poped.db$cfaa*matrix(1,poped.db$m,poped.db$na)
aa

optxt=TRUE
axt=optxt*poped.db$cfaxt*matrix(1,poped.db$m,poped.db$maxni)
axt

calc_ofv_and_grad(x=c(poped.db$gxt,poped.db$ga),
                  optxt=optxt, opta=opta, 
                  model_switch=poped.db$global_model_switch,
                  aa=aa,
                  axt=axt,
                  groupsize=poped.db$groupsize,
                  ni=poped.db$gni,
                  xtopto=poped.db$gxt,
                  xopto=poped.db$gx,
                  aopto=poped.db$ga,
                  bpop=poped.db$param.pt.val$bpop,
                  d=poped.db$param.pt.val$d,
                  sigma=poped.db$param.pt.val$sigma,
                  docc_full=poped.db$param.pt.val$docc,
                  poped.db,
                  only_fim=FALSE)

\dontrun{
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)
  
}




