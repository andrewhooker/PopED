# warfarin ed model

# very few samples
poped.db$ED_samp_size=10
ed_mftot(model_switch=poped.db$global_model_switch,
         groupsize=poped.db$groupsize,
         ni=poped.db$gni,
         xtoptn=poped.db$gxt,
         xoptn=poped.db$gx,
         aoptn=poped.db$ga,
         bpopdescr=poped.db$gbpop,
         ddescr=poped.db$gd,
         covd=poped.db$covd,
         sigma=poped.db$sigma,
         docc=poped.db$docc, 
         poped.db)["ED_ofv"]

  
