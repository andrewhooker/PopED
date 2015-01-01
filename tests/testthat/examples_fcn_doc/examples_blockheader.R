# warfarin optimization model 


FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)

blockheader(poped.db,name="")

blockheader(name="",iter=1,poped.db)


blockheader(name='',
              iter=1,
              poped.db,
              e_flag=FALSE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$settings$optsw[4],
              opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$parameters$param.pt.val$bpop,
              d=poped.db$parameters$param.pt.val$d,
              docc=poped.db$parameters$docc,sigma=poped.db$parameters$param.pt.val$sigma)



blockheader(name='',
              iter=1,
              poped.db,
              e_flag=TRUE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$settings$optsw[4],
              opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$parameters$param.pt.val$bpop,
              d=poped.db$parameters$param.pt.val$d,
              docc=poped.db$parameters$docc,sigma=poped.db$parameters$param.pt.val$sigma)
  
  


