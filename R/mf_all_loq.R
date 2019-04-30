mf_all_loq <- function(model_switch_i,xt_i,x_i,a_i,bpop_val,d_full,sigma_full,docc_full,
                       poped.db,
                       loq, # vector of length number of models
                       loq_method=1,#poped.db$settings$loq_method,
                       loq_PI_conf_level = 0.95,#poped.db$settings$loq_PI_conf_level,
                       loq_prob_limit = 0.001,#poped.db$settings$loq_prob_limit,
                       uloq = NULL,
                       ...){

  verbose=FALSE
  
  # PRED calculations based on FO
  b_ind = poped.db$parameters$b_global[,1,drop=F]*0
  bocc_ind = poped.db$parameters$bocc_global[[1]]*0
  g0 = feval(poped.db$model$fg_pointer,x_i,a_i,bpop_val,b_ind,bocc_ind)
  pred <- feval(poped.db$model$ff_pointer,model_switch_i,xt_i,g0,poped.db)
  pred <- drop(pred[[1]])
  
  fim_size <- get_fim_size(poped.db)
  
  if(loq_method==1){ # D6 method
    
    # COV calculations based on FO
    cov <- v(model_switch_i,xt_i,
             x_i,a_i,bpop_val,b_ind,bocc_ind,d_full,
             sigma_full,docc_full,poped.db)[[1]]
    
    # compute points that have PI that overlaps LOQ
    PI_alpha <- 1-loq_PI_conf_level
    z_val <- qnorm(1-PI_alpha/2)
    se_val <- sqrt(diag(cov))
    ci_u <- pred + z_val*se_val 
    ci_l <- pred - z_val*se_val 
    
    # df <- tibble::tibble(pred=c(pred),ci_l=c(ci_l),ci_u=c(ci_u),loq=loq)
    # df <- df %>% dplyr::mutate(above=dplyr::if_else(ci_l>loq,1,0)) %>% 
    #   dplyr::mutate(below=dplyr::if_else(ci_u<loq,1,0)) %>% 
    #   dplyr::mutate(overlap=dplyr::if_else(ci_u>loq & ci_l<loq,1,0)) %>% 
    #   dplyr::mutate(bloq_obs=dplyr::if_else(below==1 & overlap==0,1,dplyr::if_else(overlap==1,2,0)))
    
    loq_full <- rep(0,length(pred))
    for(k in 1:length(loq)){
      loq_full[model_switch_i==k] <- loq[k]
    }
    
    above <- rep(0,length(pred))
    above[ci_l>loq_full] <- 1
    
    below <- 0*above
    below[ci_u<loq_full] <- 1
    
    overlap <- 0*above
    overlap[ci_u>loq_full&ci_l<loq_full] <- 1
    
    bloq_obs_master <- 0*above + 2
    bloq_obs_master[below==1 & overlap==0] <- 1
    bloq_obs_master[above==1 & overlap==0] <- 0
    
    
    #bloq_obs_master <- df$bloq_obs
    #bloq_obs_master <- bloq_obs_master*0+2
    
    
    # number of potential bloq_obs
    n_pot_blq <- sum(bloq_obs_master==2)
    
    if(n_pot_blq>0){
      
      # combination of potential bloq_obs with datapoints above LOQ
      bloq_obs <- gtools::permutations(2,n_pot_blq,v=c(0,1),repeats.allowed=TRUE)
      loq_pot_blq <- matrix(rep(loq_full[bloq_obs_master==2],nrow(bloq_obs)),byrow=T)
      bloq_comb_l <- zeros(size(bloq_obs))
      bloq_comb_l[bloq_obs==1] <- -Inf
      bloq_comb_l[bloq_obs==0] <- loq_pot_blq[bloq_obs==0]
      bloq_comb_u <- zeros(size(bloq_obs))
      bloq_comb_u[bloq_obs==1] <- loq_pot_blq[bloq_obs==1]
      bloq_comb_u[bloq_obs==0] <- Inf
    
      # compute all probabilities
      pred_pot_blq <- pred[bloq_obs_master==2]
      cov_pot_blq <- cov[bloq_obs_master==2,bloq_obs_master==2]
      p_bloq_comb <- rep(0,nrow(bloq_obs))
      p_bloq_comb_full <- rep(0,nrow(bloq_obs)) # for diagnostics
      for(j in 1:nrow(bloq_obs)){
        p_bloq_comb_tmp <- 
          mvtnorm::pmvnorm(bloq_comb_l[j,],
                           bloq_comb_u[j,], 
                           mean=pred_pot_blq, 
                           sigma=cov_pot_blq)
        #p_bloq_comb_tmp <- mnormt::sadmvn(bloq_comb_l[j,],bloq_comb_u[j,], pred, cov)
        
        p_bloq_comb_full[j] <- p_bloq_comb_tmp # for diagnostics
        
        # filter out low probability values
        if(p_bloq_comb_tmp < loq_prob_limit) p_bloq_comb_tmp <- 0
        p_bloq_comb[j] <- p_bloq_comb_tmp
      }
      
      # rescale probabilities
      p_bloq_comb <- p_bloq_comb/sum(p_bloq_comb)
      
      if(verbose){
        bloq_obs_tmp <- bloq_obs_master 
        model <- xt <- NULL
        for(j in 1:nrow(bloq_obs)){
          bloq_obs_tmp[bloq_obs_master==2] <- bloq_obs[j,] 
          df_p <- tibble::tibble(model=c(model_switch_i),xt=c(xt_i),pred=c(pred),BLQ=bloq_obs_tmp)
          df_p <- df_p %>% dplyr::arrange(model,xt)
          # print(df_p)
          cat("Time: ",sprintf("%1.f",df_p$xt),
              "\nBLQ: ",sprintf("%1.f",df_p$BLQ), 
              "\np_initial: ", sprintf("%8.4g",p_bloq_comb_full[j]),
              " p_final: ",sprintf("%8.4g",p_bloq_comb[j]),
              "\n\n")
        }
      }
      
      # compute FIM for each case and combine
      fim <- zeros(fim_size)
      bloq_obs_tmp <- bloq_obs_master 
      for(j in 1:nrow(bloq_obs)){
        #j=2
        bloq_obs_tmp[bloq_obs_master==2] <- bloq_obs[j,] 
        if(any(bloq_obs_tmp==0) & p_bloq_comb[j]!=0){
          fim_tmp <- mf_all(model_switch_i[bloq_obs_tmp==0,1,drop=F],
                            xt_i[bloq_obs_tmp==0,1,drop=F],
                            x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped.db)$ret
          fim <- fim + p_bloq_comb[j]*fim_tmp
        }
        
      }
    } else {
      bloq_obs <- bloq_obs_master
      fim <- zeros(fim_size)
      if(any(bloq_obs==0)){
        fim <- mf_all(model_switch_i[bloq_obs==0,1,drop=F],
                      xt_i[bloq_obs==0,1,drop=F],
                      x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped.db)$ret
      }
    }
  }
  if(loq_method==2){ #D2 method
    loq_full <- rep(0,length(pred))
    for(k in 1:length(loq)){
      loq_full[model_switch_i==k] <- loq[k]
    }
    
    if(!is.null(uloq)){
      uloq_full <- rep(0,length(pred))
      for(k in 1:length(uloq)){
        uloq_full[model_switch_i==k] <- uloq[k]
      }
    }
    
    bloq_obs <- pred<loq_full
    if(!is.null(uloq)){
      uloq_obs <- pred>uloq_full
      bloq_obs <- bloq_obs | uloq_obs
    } 
    fim <- zeros(fim_size)
    if(any(bloq_obs==0)){
      fim <- mf_all(model_switch_i[bloq_obs==0,1,drop=F],
                    xt_i[bloq_obs==0,1,drop=F],
                    x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped.db)$ret
    }
  }
  
  return(list(fim=fim,poped.db=poped.db)) 
}