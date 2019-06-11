mf_all_loq <- function(model_switch_i,xt_i,x_i,a_i,bpop_val,d_full,sigma_full,docc_full,
                       poped.db,
                       loq = -Inf, # vector of length number of models
                       loq_method=1,#poped.db$settings$loq_method,
                       loq_PI_conf_level = 0.95,#poped.db$settings$loq_PI_conf_level,
                       loq_prob_limit = 0.001,#poped.db$settings$loq_prob_limit,
                       loq_start_time = NULL,
                       uloq = Inf,
                       uloq_method=1,
                       uloq_start_time = NULL,
                       verbose=FALSE,
                       ...){

  
  # TODO: add to poped.db
  
  # PRED calculations based on FO
  b_ind = poped.db$parameters$b_global[,1,drop=F]*0
  bocc_ind = poped.db$parameters$bocc_global[[1]]*0
  g0 = feval(poped.db$model$fg_pointer,x_i,a_i,bpop_val,b_ind,bocc_ind)
  pred <- feval(poped.db$model$ff_pointer,model_switch_i,xt_i,g0,poped.db)
  pred <- drop(pred[[1]])
  
  fim_size <- get_fim_size(poped.db)
  
  n_mod <- unique(c(model_switch_i))
  
  loq_full <- rep(NA,length(pred))
  uloq_full <- rep(NA,length(pred))
  
  if(length(loq)==1) loq_full <- rep(loq,length(pred))
  if(length(uloq)==1) uloq_full <- rep(uloq,length(pred))
  
  if(length(loq) == n_mod){
    for(k in unique(c(model_switch_i))){
      loq_full[model_switch_i==k] <- loq[k]
    }
  }
  
  if(length(uloq) == n_mod){
    for(k in unique(c(model_switch_i))){
      uloq_full[model_switch_i==k] <- uloq[k]
    }
  }
  
  if(!is.null(loq_start_time)) loq_full[xt_i<loq_start_time] <- -Inf
  if(!is.null(uloq_start_time)) uloq_full[xt_i<uloq_start_time] <- Inf
  
  if(any(is.na(loq_full)) | any(is.na(uloq_full))) stop("loq or uloq not specified properly") 
  
  # D2 method
  if(uloq_method==2) uloq_obs_master <- pred>uloq_full
  if(loq_method==2) bloq_obs_master <- pred<loq_full
  
  # D6 method
  if(loq_method==1 | uloq_method==1){ 
    
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
    
    if(loq_method==1){
      above <- below <- overlap <- loq_full*0
      
      above[ci_l>loq_full] <- 1
      below[ci_u<loq_full] <- 1
      overlap[ci_u>loq_full&ci_l<loq_full] <- 1
      
      bloq_obs_master <- 0*above + 2
      bloq_obs_master[below==1 & overlap==0] <- 1
      bloq_obs_master[above==1 & overlap==0] <- 0
    }

    if(uloq_method==1){
      above_u <- below_u <- overlap_u <- uloq_full*0
      above_u[ci_l>uloq_full] <- 1
      below_u[ci_u<uloq_full] <- 1
      overlap_u[ci_u>uloq_full&ci_l<uloq_full] <- 1
      
      uloq_obs_master <- 0*above_u + 2
      uloq_obs_master[below_u==1 & overlap_u==0] <- 0
      uloq_obs_master[above_u==1 & overlap_u==0] <- 1
    } 
  }
  
  #bloq_obs_master <- df$bloq_obs
  #bloq_obs_master <- bloq_obs_master*0+2
  
  loq_obs_master <- bloq_obs_master
  loq_obs_master[uloq_obs_master==1] <- 1
  loq_obs_master[uloq_obs_master==2 & bloq_obs_master!=1] <- 2
  
  # number of potential loq_obs
  n_pot_loq <- sum(loq_obs_master==2)
  
  if(n_pot_loq>0){ # D6 Method
    
    # combination of potential obs with datapoints above LOQ or below ULOQ
    loq_obs_init <- gtools::permutations(2,n_pot_loq,v=c(0,1),repeats.allowed=TRUE)
    
    # map for type of observation
    # 0 = normal observations
    # 1 = BLOQ
    # 2 = ULOQ
    # 3 = Could be either BLOQ or ULOQ (and need expanding)
    loq_obs_map <- loq_obs_master[loq_obs_master==2]*0 + 1
    uloq_obs_map <- uloq_obs_master[loq_obs_master==2]
    bloq_obs_map <- bloq_obs_master[loq_obs_master==2]
    loq_obs_map[uloq_obs_map==2 & bloq_obs_map!=2] <- 2 
    loq_obs_map[uloq_obs_map==2 & bloq_obs_map==2] <- 3 
    
    loq_obs_short <- matrix(loq_obs_map,ncol=length(loq_obs_map),nrow=nrow(loq_obs_init),byrow=T)
    loq_obs_short[loq_obs_init==0] <- 0
    
    
    # expand rows that could be BLOQ or ULOQ 
    exp_rows <- apply(loq_obs_short==3,1,any)
    loq_obs <- loq_obs_short[!exp_rows,,drop=F]
    if(any(exp_rows)){
      loq_obs_tmp <- loq_obs_short[exp_rows,]
      
      # expand rows
      for(i in 1:nrow(loq_obs_tmp)){
        #i <- 1
        obs_tmp <- loq_obs_tmp[i,]
        perm_tmp <- gtools::permutations(2,sum(obs_tmp==3),v=c(1,2),repeats.allowed=TRUE)
        
        obs_tmp_exp <- matrix(obs_tmp,ncol=length(obs_tmp),nrow=nrow(perm_tmp),byrow=T)
        obs_tmp_exp[obs_tmp_exp==3] <- perm_tmp[,]
        loq_obs <- rbind(loq_obs,obs_tmp_exp)
      }
    }
    # make sure that mapped values are all accounted for 
    if(any(loq_obs==3)) stop("Combinations not fully expanded")
    
    # cat(loq_obs,"\n")
    # cat(loq_obs_master,"\n")
    # if(sum(loq_obs_master==2)==1) browser()
    lloq_mat <- matrix(rep(loq_full[loq_obs_master==2],nrow(loq_obs)),nrow=nrow(loq_obs),byrow=T)
    uloq_mat <- matrix(rep(uloq_full[loq_obs_master==2],nrow(loq_obs)),nrow=nrow(loq_obs),byrow=T)
    
    
    loq_comb_l <- loq_obs*NA
    loq_comb_u <- loq_obs*NA
    
    # BLOQ
    loq_comb_l[loq_obs==1] <- -Inf
    loq_comb_u[loq_obs==1] <- lloq_mat[loq_obs==1]
    
    # ULOQ
    loq_comb_l[loq_obs==2] <- uloq_mat[loq_obs==2]
    loq_comb_u[loq_obs==2] <- Inf
    
    # normal observations
    loq_comb_l[loq_obs==0] <- lloq_mat[loq_obs==0]
    loq_comb_u[loq_obs==0] <- uloq_mat[loq_obs==0]
    
    # compute all probabilities
    pred_pot_loq <- pred[loq_obs_master==2]
    cov_pot_loq <- cov[loq_obs_master==2,loq_obs_master==2]
    p_loq_comb <- rep(0,nrow(loq_obs))
    p_loq_comb_full <- rep(0,nrow(loq_obs)) # for diagnostics
    for(j in 1:nrow(loq_obs)){
      p_loq_comb_tmp <- 
        mvtnorm::pmvnorm(loq_comb_l[j,],
                         loq_comb_u[j,], 
                         mean=pred_pot_loq, 
                         sigma=cov_pot_loq)
      #p_bloq_comb_tmp <- mnormt::sadmvn(bloq_comb_l[j,],bloq_comb_u[j,], pred, cov)
      
      # filter out low probability values
      p_loq_comb_full[j] <- p_loq_comb_tmp # save initial probs for diagnostics
      if(p_loq_comb_tmp < loq_prob_limit) p_loq_comb_tmp <- 0
      p_loq_comb[j] <- p_loq_comb_tmp
    }
    
    # sum of probabilities
    tot_p <- sum(p_loq_comb_full)
    max_diff <- PI_alpha/2*length(loq_obs_master==2) # max p missed if all points are truncated with PI 
    if(tot_p > 1.01 | tot_p < (1-max_diff)) stop("Sum of initial probabilities: ", sprintf("%6.5g",tot_p),"\n",
                                                 "Probabilities do not add up to one!")
    
    # rescale probabilities
    p_loq_comb <- p_loq_comb/sum(p_loq_comb)
    
    if(verbose){
      loq_obs_tmp <- loq_obs_master 
      model <- xt <- NULL
      for(j in 1:nrow(loq_obs)){
        loq_obs_tmp[loq_obs_master==2] <- loq_obs[j,] 
        df_p <- tibble::tibble(model=c(model_switch_i),xt=c(xt_i),pred=c(pred),LOQ=loq_obs_tmp)
        df_p <- df_p %>% dplyr::arrange(model,xt)
        #print(df_p)
        cat("Time: ",sprintf("%1.f",df_p$xt),
            "\nLOQ: ",sprintf("%1.f",df_p$LOQ), 
            "\np_initial: ", sprintf("%8.4g",p_loq_comb_full[j]),
            " p_final: ",sprintf("%8.4g",p_loq_comb[j]),
            "\n\n")
      }
      cat("sum of initial probabilities: ", sprintf("%6.5g",tot_p),"\n")
      cat("sum of final probabilities: ", sprintf("%6.5g",sum(p_loq_comb)),"\n")
      cat("\n")
    }
    
    # compute FIM for each case and combine
    fim <- zeros(fim_size)
    loq_obs_tmp <- loq_obs_master 
    for(j in 1:nrow(loq_obs)){
      #j=2
      loq_obs_tmp[loq_obs_master==2] <- loq_obs[j,] 
      if(any(loq_obs_tmp==0) & p_loq_comb[j]!=0){
        fim_tmp <- mf_all(model_switch_i[loq_obs_tmp==0,1,drop=F],
                          xt_i[loq_obs_tmp==0,1,drop=F],
                          x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped.db)$ret
        fim <- fim + p_loq_comb[j]*fim_tmp
      }
      
    }
  } else { # D2 method for BLOQ
    fim <- zeros(fim_size)
    if(any(loq_obs_master==0)){
      fim <- mf_all(model_switch_i[loq_obs_master==0,1,drop=F],
                    xt_i[loq_obs_master==0,1,drop=F],
                    x_i,a_i,bpop_val,d_full,sigma_full,docc_full,poped.db)$ret
    }
  }
  
  return(list(fim=fim,poped.db=poped.db)) 
}