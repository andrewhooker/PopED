# warfarin optimization model

#for the FO approximation
ind=1

# no occasion defined in this example, so result is zero
output <- mf5(model_switch=t(poped.db$global_model_switch[ind,,drop=FALSE]),
   xt=t(poped.db$gxt[ind,,drop=FALSE]),
   x=zeros(0,1),
   a=t(poped.db$ga[ind,,drop=FALSE]),
   bpop=poped.db$gbpop[,2,drop=FALSE],
   d=poped.db$param.pt.val$d,
   sigma=poped.db$sigma,
   docc=poped.db$param.pt.val$docc,
   poped.db)

# in this simple case the full FIM is just the sum of the individual FIMs
# and all the individual FIMs are the same
det(output$ret*32) == det(evaluate.fim(poped.db,fim.calc.type=4))  
