
create_ofv <- function(poped.db,
                       opt_xt=poped.db$settings$optsw[2],
                       opt_a=poped.db$settings$optsw[4],
                       opt_x=poped.db$settings$optsw[3],
                       opt_samps=poped.db$settings$optsw[1],
                       opt_inds=poped.db$settings$optsw[5],
                       fim.calc.type=poped.db$settings$iFIMCalculationType,
                       ofv_calc_type=poped.db$settings$ofv_calc_type,
                       approx_type=poped.db$settings$iApproximationMethod,
                       d_switch=poped.db$settings$d_switch,
                       ED_samp_size = poped.db$settings$ED_samp_size,
                       bLHS=poped.db$settings$bLHS,
                       use_laplace=poped.db$settings$iEDCalculationType,
                       ofv_fun = poped.db$settings$ofv_fun,
                       transform_parameters=T,
                       ...){
  
  #------------ update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
    }
  }
  
  #----------- checks
  if((sum(poped.db$settings$optsw)==0)){
    stop('No optimization parameter is set.')
  }
  
  if(is.null(ofv_fun) || is.function(ofv_fun)){
    ofv_fun_user <- ofv_fun 
  } else {
    # source explicit file
    # here I assume that function in file has same name as filename minus .txt and pathnames
    if(file.exists(as.character(ofv_fun))){
      source(as.character(ofv_fun))
      ofv_fun_user <- eval(parse(text=fileparts(ofv_fun)[["filename"]]))
    } else {
      stop("ofv_fun is not a function or NULL, and no file with that name was found")
    }
    
  }
  if(!is.null(ofv_fun)){
    poped.db$settings$ofv_calc_type = 0
  }
  
  # Collect the parameters to optimize
  par <- c()
  upper <- c()
  lower <- c()
  par_grouping <- c()
  par_type <- c()
  par_dim <- list()
  allowed_values <- NULL
  build_allowed_values <- FALSE
  if(!is.null(poped.db$design_space$xt_space) ||
     !is.null(poped.db$design_space$a_space)) build_allowed_values <- TRUE
  if(opt_samps) stop('Sample number optimization is not yet implemented in the R-version of PopED.')
  if(opt_inds) stop('Optimization  of number of individuals in different groups is not yet implemented in the R-version of PopED.')
  if(opt_xt){ 
    #par <- c(par,poped.db$design$xt)
    # upper <- c(upper,poped.db$design_space$maxxt)
    # lower <- c(lower,poped.db$design_space$minxt)
    # par_grouping <- c(par_grouping,poped.db$design_space$G_xt)
    # par_type <- c(par_type,rep("xt",length(poped.db$design$xt)))
    par_dim$xt <- dim(poped.db$design$xt)
    if(is.null(poped.db$design_space$xt_space) && build_allowed_values){ 
      poped.db$design_space$xt_space <- cell(par_dim$xt)
    }
    # allowed_values <- c(allowed_values,poped.db$design_space$xt_space)
    
    for(i in 1:poped.db$design$m){
      if((poped.db$design$ni[i]!=0 && poped.db$design$groupsize[i]!=0)){
        par <- c(par,poped.db$design$xt[i,1:poped.db$design$ni[i]])
        upper <- c(upper,poped.db$design_space$maxxt[i,1:poped.db$design$ni[i]])
        lower <- c(lower,poped.db$design_space$minxt[i,1:poped.db$design$ni[i]])
        par_grouping <- c(par_grouping,poped.db$design_space$G_xt[i,1:poped.db$design$ni[i]])
        par_type <- c(par_type,rep("xt",length(poped.db$design$xt[i,1:poped.db$design$ni[i]])))
        allowed_values <- c(allowed_values,poped.db$design_space$xt_space[i,1:poped.db$design$ni[i]])
      }
    }
  }
  if(opt_a) { 
    par <- c(par,poped.db$design$a)
    upper <- c(upper,poped.db$design_space$maxa)
    lower <- c(lower,poped.db$design_space$mina)
    if(opt_xt){
      par_grouping <- c(par_grouping,poped.db$design_space$G_a + max(par_grouping)) 
    } else {
      par_grouping <- c(par_grouping,poped.db$design_space$G_a) 
    }
    par_type <- c(par_type,rep("a",length(poped.db$design$a)))
    par_dim$a <- dim(poped.db$design$a)
    if(is.null(poped.db$design_space$a_space) && build_allowed_values){ 
      poped.db$design_space$a_space <- cell(par_dim$a)
    }
    allowed_values <- c(allowed_values,poped.db$design_space$a_space)
    
  }
  if(opt_x) NULL # par <- c(par,poped.db$design$x)
  
  # continuous and discrete parameters
  npar <- max(c(length(lower),length(upper),length(allowed_values),length(par)))
  par_cat_cont <- rep("cont",npar)
  if(!is.null(allowed_values)){
    for(k in 1:npar){
      if(!is.na(allowed_values[[k]]) && length(allowed_values[[k]]>0)){
        par_cat_cont[k] <- "cat"          
      }
    }
  }
  
  # Parameter grouping
  par_df <- data.frame(par,par_grouping,upper,lower,par_type,par_cat_cont)
  par_df_unique <- NULL
  #allowed_values_full <- allowed_values 
  if(!all(!duplicated(par_df$par_grouping))){
    par_df_unique <- par_df[!duplicated(par_df$par_grouping),]
    par <- par_df_unique$par
    lower <- par_df_unique$lower
    upper <- par_df_unique$upper
    par_cat_cont <- par_df_unique$par_cat_cont
    allowed_values <- allowed_values[!duplicated(par_df$par_grouping)]
  }
  
  par_df_2 <- data.frame(par,upper,lower,par_cat_cont)
  #par_fixed_index <- which(upper==lower)
  par_fixed_index <- which(upper==lower & par_cat_cont=="cont")
  for(npar in 1:length(par)){
    if(par_cat_cont[npar]=="cont") next
    if(all(par[npar]==allowed_values[[npar]])){
      par_fixed_index <- c(par_fixed_index,npar)
    } 
  }
  par_fixed_index <- sort(unique(par_fixed_index))
  par_not_fixed_index <- 1:length(par)
  par_not_fixed_index <- par_not_fixed_index[!(par_not_fixed_index %in% par_fixed_index)]
  
  if(length(par_fixed_index)!=0){
    par <- par[-c(par_fixed_index)]
    npar <- length(par)
    lower <- lower[-c(par_fixed_index)]
    upper <- upper[-c(par_fixed_index)]
    par_cat_cont <- par_cat_cont[-c(par_fixed_index)]
    allowed_values <- allowed_values[-c(par_fixed_index)]
  }
  
  if(length(par)==0) stop("No design parameters have a design space to optimize")
  
  # if(length(par)==0){
  #   message("No design parameters have a design space to optimize")
  #   return(invisible(list( ofv= output$ofv, FIM=fmf, poped.db = poped.db )))
  # } 
  
  if(!is.null(allowed_values)){
    for(k in 1:npar){
      if(!all(is.na(allowed_values[[k]]))){
        if(length(upper)>0) allowed_values[[k]] <- allowed_values[[k]][allowed_values[[k]]<=upper[k]]
        if(length(lower)>0) allowed_values[[k]] <- allowed_values[[k]][allowed_values[[k]]>=lower[k]]
      }
    }
  }
  
  lower_opt <- lower
  upper_opt <- upper
  if(transform_parameters){
    for(i in 1:length(par)){
      if(par_cat_cont[i]=="cont"){
        # par[i] <- FastImputation::NormalizeBoundedVariable(par[i],
        #                                                    constraints = 
        #                                                      list(lower=lower[i],
        #                                                           upper=upper[i]))
        par[i] <- unbound_par(par[i],lower=lower[i],upper=upper[i])
        lower_opt[i] = -Inf
        upper_opt[i] = Inf
        
        #     if(is.finite(lower[i]) && is.finite(upper[i])){
        #       par[i] <- (par[i] - lower[i])/(upper[i]-lower[i])
        #       par[i] <- stats::qlogis(par[i])
        #     }
        #     
      }
    }
  }
  
  #------- create optimization function with optimization parameters first
  ofv_fun <- function(par,only_cont=F,...){
    
    if(transform_parameters){
      for(i in 1:length(par)){
        if(par_cat_cont[i]=="cont"){
          # par[i] <- FastImputation::BoundNormalizedVariable(par[i],
          #                                                   constraints = 
          #                                                     list(lower=lower[i],
          #                                                          upper=upper[i]))
          par[i] <- bound_par(par[i],lower=lower[i],upper=upper[i])
          # if(is.finite(lower[i]) && is.finite(upper[i])){
          #   par[i] <- stats::plogis(par[i])
          #   par[i] <- (par[i] )*(upper[i]-lower[i])+ lower[i]
          # }
          
        }
      }
    }
    
    if(length(par_fixed_index)!=0){
      par_df_2[par_not_fixed_index,"par"] <- par
      par <- par_df_2$par
    }
    
    if(!is.null(par_df_unique)){
      if(only_cont){ 
        par_df_unique[par_df_unique$par_cat_cont=="cont","par"] <- par
      } else {
        par_df_unique$par <- par
      }
      for(j in par_df_unique$par_grouping){
        par_df[par_df$par_grouping==j,"par"] <- par_df_unique[par_df_unique$par_grouping==j,"par"]
      }  
      
      #par_full[par_cat_cont=="cont"] <- par 
      par <- par_df$par
    } else if (only_cont){ 
      par_df[par_df$par_cat_cont=="cont","par"] <- par
      par <- par_df$par
    }
    xt <- NULL
    #if(opt_xt) xt <- matrix(par[par_type=="xt"],par_dim$xt)
    if(opt_xt){
      xt <- zeros(par_dim$xt)
      par_xt <- par[par_type=="xt"]
      for(i in 1:poped.db$design$m){
        if((poped.db$design$ni[i]!=0 && poped.db$design$groupsize[i]!=0)){
          xt[i,1:poped.db$design$ni[i]] <- par_xt[1:poped.db$design$ni[i]]
          par_xt <- par_xt[-c(1:poped.db$design$ni[i])]
        }
      }
    } 
    
    a <- NULL
    if(opt_a) a <- matrix(par[par_type=="a"],par_dim$a)
    
    # if(d_switch){
    #   FIM <- evaluate.fim(poped.db,xt=xt,a=a,...)
    #   ofv <- ofv_fim(FIM,poped.db,...)
    # } else{
    #   output <-calc_ofv_and_fim(poped.db,d_switch=d_switch,
    #                             ED_samp_size=ED_samp_size,
    #                             bLHS=bLHS,
    #                             use_laplace=use_laplace,
    #                             ofv_calc_type=ofv_calc_type,
    #                             fim.calc.type=fim.calc.type,
    #                             xt=xt,
    #                             a=a,
    #                             ...)
    #   
    #   FIM <- output$fim
    #   ofv <- output$ofv
    # }
    
    
    extra_args <- dots(...)
    extra_args$evaluate_fim <- FALSE
    
    output <- do.call(calc_ofv_and_fim,
                      c(list(
                        poped.db,
                        d_switch=d_switch,
                        ED_samp_size=ED_samp_size,
                        bLHS=bLHS,
                        use_laplace=use_laplace,
                        ofv_calc_type=ofv_calc_type,
                        fim.calc.type=fim.calc.type,
                        xt=xt,
                        a=a,
                        ofv_fun = ofv_fun_user
                      ),
                      extra_args))
    
    
    # output <-calc_ofv_and_fim(poped.db,d_switch=d_switch,
    #                           ED_samp_size=ED_samp_size,
    #                           bLHS=bLHS,
    #                           use_laplace=use_laplace,
    #                           ofv_calc_type=ofv_calc_type,
    #                           fim.calc.type=fim.calc.type,
    #                           xt=xt,
    #                           a=a,
    #                           evaluate_fim = F,
    #                           ...)
    #FIM <- output$fim
    ofv <- output$ofv
    
    
    #ofv <- tryCatch(ofv_fim(FIM,poped.db,...), error = function(e) e)
    if(!is.finite(ofv) && ofv_calc_type==4){
      #ofv <- -Inf 
      ofv <- NA
    } else {
      #if(!is.finite(ofv)) ofv <- 1e-15
      if(!is.finite(ofv)) ofv <- NA
      #if(!is.finite(ofv)) ofv <- -Inf
    }
    
    #cat(ofv,"\n")
    return(ofv)
  }
  
  back_transform_par <- function(par,only_cont=F,...){
    if(transform_parameters){
      for(i in 1:length(par)){
        if(par_cat_cont[i]=="cont"){
          # par[i] <- FastImputation::BoundNormalizedVariable(par[i],
          #                                                   constraints = 
          #                                                     list(lower=lower[i],
          #                                                          upper=upper[i]))
          par[i] <- bound_par(par[i],lower=lower[i],upper=upper[i])
          # if(is.finite(lower[i]) && is.finite(upper[i])){
          #   par[i] <- stats::plogis(par[i])
          #   par[i] <- (par[i] )*(upper[i]-lower[i])+ lower[i]
          # }
          
        }
      }
    }
    return(par)
  }
  
  
  return(list(fun=ofv_fun,
              par=par,
              back_transform_par=back_transform_par, 
              space=list(lower=lower_opt,upper=upper_opt,par_cat_cont=par_cat_cont,
                         par_fixed_index=par_fixed_index,par_df_unique=par_df_unique,
                         par_df=par_df, par_dim=par_dim,
                         par_type=par_type, par_not_fixed_index = par_not_fixed_index)) )
}

