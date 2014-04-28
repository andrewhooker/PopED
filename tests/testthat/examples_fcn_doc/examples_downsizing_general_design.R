# warfarin basic model

dgd <- downsizing_general_design(poped.db)

output = mftot(dgd$model_switch,poped.db$groupsize,dgd$ni,dgd$xt,dgd$x,dgd$a,
               poped.db$param.pt.val$bpop,poped.db$param.pt.val$d,
               poped.db$param.pt.val$sigma,poped.db$param.pt.val$docc,poped.db)
FIM <- output$ret
det(FIM)
