# warfarin optimization model 


FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)

blockheader_2(name="",iter=1,poped.db)


blockheader_2(name='',
              iter=1,
              poped.db,
              e_flag=FALSE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$optsw[4],
              opt_samps=poped.db$optsw[1],opt_inds=poped.db$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$param.pt.val$bpop,
              d=poped.db$param.pt.val$d,
              docc=poped.db$docc,sigma=poped.db$param.pt.val$sigma)




  
  


