# warfarin optimization model

#for the FO approximation
ind=1
LinMatrixL(model_switch=t(poped.db$global_model_switch[ind,,drop=FALSE]),
          xt_ind=t(poped.db$gxt[ind,,drop=FALSE]),
          x=zeros(0,1),
          a=t(poped.db$ga[ind,,drop=FALSE]),
          bpop=poped.db$gbpop[,2,drop=FALSE],
          b_ind=zeros(poped.db$NumRanEff,1),
          bocc_ind=zeros(poped.db$NumDocc,1),
          poped.db)["y"]

