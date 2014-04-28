# warfarin optimization model 


\dontrun{  
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=0,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)
  
  f_name <- 'calc_ofv_and_grad' 
  gen_des <- downsizing_general_design(poped.db)
  
  aa <- 0*poped.db$cfaa*matrix(1,poped.db$m,poped.db$na)
  axt=1*poped.db$cfaxt*matrix(1,poped.db$m,poped.db$maxni)
  
  f_options_1 <- list(gen_des$x,1, 0, gen_des$model_switch,
                    aa=aa,axt=axt,poped.db$groupsize,
                    gen_des$ni,
                    gen_des$xt,gen_des$x,gen_des$a,gen_des$bpop[,2,drop=F],
                    getfulld(gen_des$d[,2,drop=F],poped.db$covd),
                    poped.db$sigma,
                    getfulld(poped.db$docc[,2,drop=F],poped.db$covdocc),poped.db)
  
  options=list('factr'=poped.db$BFGSConvergenceCriteriaMinStep,
               #'factr'=0.01,
               'pgtol'=poped.db$BFGSProjectedGradientTol,
               'ftol'=poped.db$BFGSTolerancef,
               'gtol'=poped.db$BFGSToleranceg,
               'xtol'=poped.db$BFGSTolerancex)
  
  x_k=t(gen_des$xt)
  lb=t(gen_des$minxt)
  ub=t(gen_des$maxxt)
  
  output <- bfgsb_min(f_name,f_options, x_k,lb,ub,options) 
  
}

