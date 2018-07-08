par_and_space_tbl <- function(poped.db) {
  df <- NULL
  allowed_values <- allowed_values_2 <- cont <- lower <- upper <- par <- NULL
  
  design=poped.db$design
  design_space=poped.db$design_space
  
  ############################
  # Collect the xt parameters 
  ############################
  df_xt <- NULL
  sel_mat_xt <- zeros(size(design$xt))
  for(i in 1:length(design$ni)){
    if((design$ni[i]!=0 && design$groupsize[i]!=0)) sel_mat_xt[i,1:(design$ni[i])] <- 1
  }
  sel_mat_xt <- sel_mat_xt==1
  
  df_xt <- tibble::tibble(par=t(design$xt)[t(sel_mat_xt)])
  df_xt <- tibble::add_column(df_xt,lower=t(design_space$minxt)[t(sel_mat_xt)])
  df_xt <- tibble::add_column(df_xt,upper=t(design_space$maxxt)[t(sel_mat_xt)])
  df_xt <- tibble::add_column(df_xt,grouping=t(design_space$G_xt)[t(sel_mat_xt)])
  df_xt <- tibble::add_column(df_xt,type="xt")
  if(is.null(design_space$xt_space)){
    df_xt <- tibble::add_column(df_xt,allowed_values=list(NA))
  } else{
    df_xt <- tibble::add_column(df_xt,allowed_values=t(design_space$xt_space)[t(sel_mat_xt)])
  }
  
  # group column
  group_mat <- ones(size(design$xt))*c(1:nrow(design$xt))
  df_xt <- tibble::add_column(df_xt,group=t(group_mat)[t(sel_mat_xt)])
  
  # name column
  name_mat <- matrix(colnames(design$xt),nrow=nrow(design$xt),ncol = ncol(design$xt),byrow=T)
  df_xt <- tibble::add_column(df_xt,name=t(name_mat)[t(sel_mat_xt)])
  
  # model_switch column
  df_xt <- tibble::add_column(df_xt,model=t(design$model_switch)[t(sel_mat_xt)])
  
  # Matrix index
  df_xt <- tibble::add_column(df_xt,index=which(t(sel_mat_xt)))
  
  df <- dplyr::bind_rows(df,df_xt)
  
  ############################
  # collect the a parameters
  ############################
  df_a <- NULL
  sel_mat_a <- ones(size(design$a))
  for(i in 1:length(design$ni)){
    if((design$ni[i]!=0 && design$groupsize[i]!=0)) sel_mat_a[i,] <- 1
  }
  sel_mat_a <- sel_mat_a==1
  
  df_a <- tibble::tibble(par=t(design$a)[t(sel_mat_a)])
  df_a <- tibble::add_column(df_a,lower=t(design_space$mina)[t(sel_mat_a)])
  df_a <- tibble::add_column(df_a,upper=t(design_space$maxa)[t(sel_mat_a)])
  df_a <- tibble::add_column(df_a,grouping=t(design_space$G_a)[t(sel_mat_a)])
  
  df_a <- tibble::add_column(df_a,type="a")
  if(is.null(design_space$a_space)){
    df_a <- tibble::add_column(df_a,allowed_values=list(NA))
  } else{
    df_a <- tibble::add_column(df_a,allowed_values=t(design_space$a_space)[t(sel_mat_a)])
  }
  
  # group column
  group_mat <- ones(size(design$a))*c(1:nrow(design$a))
  df_a <- tibble::add_column(df_a,group=t(group_mat)[t(sel_mat_a)])
  
  # name column
  if(!any(dim(design$a)==0) & !all(is.null(colnames(design$a)))){
    name_mat <- matrix(colnames(design$a),nrow=nrow(design$a),ncol = ncol(design$a),byrow=T)
    df_a <- tibble::add_column(df_a,name=t(name_mat)[t(sel_mat_a)])
  }
  
  # Matrix index
  df_a <- tibble::add_column(df_a,index=which(t(sel_mat_a)))
  
  df <- dplyr::bind_rows(df,df_a)
  
  ############################
  # collect the x parameters
  ############################
  df_x <- NULL
  sel_mat_x <- ones(size(design$x))
  for(i in 1:length(design$ni)){
    if((design$ni[i]!=0 && design$groupsize[i]!=0)) sel_mat_x[i,] <- 1
  }
  sel_mat_x <- sel_mat_x==1
  
  df_x <- tibble::tibble(par=t(design$x)[t(sel_mat_x)])
  df_x <- tibble::add_column(df_x,lower=t(design_space$minx)[t(sel_mat_x)])
  df_x <- tibble::add_column(df_x,upper=t(design_space$maxx)[t(sel_mat_x)])
  df_x <- tibble::add_column(df_x,grouping=t(design_space$G_x)[t(sel_mat_x)])
  
  df_x <- tibble::add_column(df_x,type="x")
  if(is.null(design_space$x_space)){
    df_x <- tibble::add_column(df_x,allowed_values=list(NA))
  } else{
    df_x <- tibble::add_column(df_x,allowed_values=t(design_space$x_space)[t(sel_mat_x)])
  }
  
  # group column
  group_mat <- ones(size(design$x))*c(1:nrow(design$x))
  df_x <- tibble::add_column(df_x,group=t(group_mat)[t(sel_mat_x)])
  
  # name column
  if(!any(dim(design$x)==0)){
    name_mat <- matrix(colnames(design$x),nrow=nrow(design$x),ncol = ncol(design$x),byrow=T)
    df_x <- tibble::add_column(df_x,name=t(name_mat)[t(sel_mat_x)])
  } 
  
  df <- dplyr::bind_rows(df,df_x)
  
  
  ############################
  # number of groups
  ############################
  df <- tibble::add_row(df,par=design$m, lower=design$m, upper=design$m, type="n_grp",name="n_grp",allowed_values=list(design$m))
  
  
  ############################
  # group-size
  ############################
  df_g <- NULL
  df_g <- tibble::tibble(par=design$groupsize[,])
  df_g <- tibble::add_column(df_g,lower=design_space$mingroupsize[,])
  df_g <- tibble::add_column(df_g,upper=design_space$maxgroupsize[,])
  if(!is.null(design_space$G_groupsize)) df_g <- tibble::add_column(df_g,grouping=design_space$G_groupsize[,])
  
  df_g <- tibble::add_column(df_g,type="n_id_grp")
  if(is.null(design_space$groupsize_space)){
    space <- c()
    for(i in 1:nrow(df_g)){
      space <- c(space,list(df_g$lower[i]:df_g$upper[i]))
    }
    df_g <- tibble::add_column(df_g,allowed_values=space)
  } else{
    df_g <- tibble::add_column(df_g,allowed_values=design_space$groupsize_space[,])
  }
  
  # group column
  group_mat <- ones(size(design$groupsize))*c(1:nrow(design$groupsize))
  df_g <- tibble::add_column(df_g,group=group_mat[,])
  
  
  # name column
  if(!any(dim(design$groupsize)==0)){
    name_mat <- matrix(colnames(design$groupsize),nrow=nrow(design$groupsize),ncol = ncol(design$groupsize),byrow=T)
    df_g <- tibble::add_column(df_g,name=name_mat[,])
  }
  
  df <- dplyr::bind_rows(df,df_g)
  
  
  df <- tibble::add_row(df,par=sum(design$groupsize), 
                        lower=design_space$mintotgroupsize, 
                        upper=design_space$maxtotgroupsize, 
                        type="n_id",
                        name="n_id",
                        allowed_values=list(design_space$mintotgroupsize:design_space$maxtotgroupsize))
  
  
  ############################
  # number of samples
  ############################
  df_g <- NULL
  df_g <- tibble::tibble(par=design$ni[,])
  df_g <- tibble::add_column(df_g,lower=design_space$minni[,])
  df_g <- tibble::add_column(df_g,upper=design_space$maxni[,])
  if(!is.null(design_space$G_ni)) df_g <- tibble::add_column(df_g,grouping=design_space$G_ni[,])
  
  df_g <- tibble::add_column(df_g,type="n_obs_grp")
  if(is.null(design_space$ni_space)){
    space <- c()
    for(i in 1:nrow(df_g)){
      space <- c(space,list(df_g$lower[i]:df_g$upper[i]))
    }
    df_g <- tibble::add_column(df_g,allowed_values=space)
  } else{
    df_g <- tibble::add_column(df_g,allowed_values=design_space$ni_space[,])
  }
  
  # group column
  group_mat <- ones(size(design$ni))*c(1:nrow(design$ni))
  df_g <- tibble::add_column(df_g,group=group_mat[,])
  
  # name column
  if(!any(dim(design$ni)==0)){
    name_mat <- matrix(colnames(design$ni),nrow=nrow(design$ni),ncol = ncol(design$ni),byrow=T)
    df_g <- tibble::add_column(df_g,name=name_mat[,])
  }
  
  df <- dplyr::bind_rows(df,df_g)
  
  
  df <- tibble::add_row(df,par=sum(design$ni), 
                        lower=design_space$mintotni, 
                        upper=design_space$maxtotni, 
                        type="n_obs",
                        name="n_obs",
                        allowed_values=list(design_space$mintotni:design_space$maxtotni))
  
  
  ############################
  # Extra stuff
  ############################
  
  # find continuous and discrete parameters
  df <- dplyr::mutate(df,allowed_values_2=ifelse(is.na(allowed_values) | sapply(allowed_values,length)==0,NA,allowed_values))
  df <- dplyr::mutate(df,allowed_values=NULL)
  df <- dplyr::rename(df,allowed_values=allowed_values_2)
  df <- dplyr::mutate(df,cont=ifelse(is.na(allowed_values),TRUE,FALSE))
  
  df <- dplyr::mutate(df,fixed=ifelse(cont,
                                      ifelse(lower==upper & upper==par,
                                             TRUE,
                                             FALSE),
                                      ifelse(sapply(allowed_values,length)==1,
                                             ifelse(mapply("%in%",par,allowed_values),
                                                    TRUE,
                                                    FALSE),
                                             FALSE)))
  
  
  
  
  # # update allowed values to be within limits
  df <- dplyr::mutate(df,allowed_values_2 = mapply("[",
                                                   allowed_values,
                                                   mapply("&",
                                                          mapply(">=",upper,allowed_values), 
                                                          mapply("<=",lower,allowed_values))))
  df <- dplyr::mutate(df,allowed_values=NULL)
  df <- dplyr::rename(df,allowed_values=allowed_values_2)
  
  # #-----------old 
  # if(!is.null(allowed_values)){
  #   for(k in 1:npar){
  #     if(!all(is.na(allowed_values[[k]]))){
  #       if(length(upper)>0) allowed_values[[k]] <- allowed_values[[k]][allowed_values[[k]]<=upper[k]]
  #       if(length(lower)>0) allowed_values[[k]] <- allowed_values[[k]][allowed_values[[k]]>=lower[k]]
  #     }
  #   }
  # }
  # #-----------old_end
  # 
  #df <- dplyr::mutate(df,allowed_values_2=allowed_values[allowed_values[])
  
  
  
  return(df)
}

get_par_and_space_optim <- function(poped.db,
                                    #df,
                                    opt_xt=poped.db$settings$optsw[2],
                                    opt_a=poped.db$settings$optsw[4],
                                    opt_x=poped.db$settings$optsw[3],
                                    opt_samps=poped.db$settings$optsw[1],
                                    opt_inds=poped.db$settings$optsw[5],
                                    transform_parameters=TRUE,
                                    cont_cat = "both",
                                    warn_when_none=TRUE) 
{
    
  type <- index <- fixed <- cont <- par <- lower <- upper <- NULL
  
  #----------- checks
  if(!any(as.logical(opt_xt),as.logical(opt_a),as.logical(opt_samps),as.logical(opt_inds))){
    stop('No optimization parameter is set.')
  }
  
  if(opt_samps) stop('Sample number optimization is not yet implemented in the R-version of PopED.')
  if(opt_inds) stop('Optimization  of number of individuals in different groups is not yet implemented in the R-version of PopED.')
  
  df <- par_and_space_tbl(poped.db)
  
  # Parameter grouping
  tmp <- df %>% dplyr::group_by(type,grouping) %>% dplyr::summarize(index=list(unique(index)))
  df <- dplyr::distinct(df,grouping,type,.keep_all=TRUE) 
  df <- dplyr::select(df,-index) 
  df <- dplyr::left_join(df,tmp,by=c("grouping","type"))
  
  
  # removed fixed parameters
  df <- dplyr::filter(df,fixed==FALSE)
  
  # remove cont or cat values
  if(cont_cat=="cont")   df <- dplyr::filter(df,cont==TRUE)
  if(cont_cat=="cat")   df <- dplyr::filter(df,cont==FALSE)
  
  if(nrow(df)==0){
    if(warn_when_none) warning("No parameters to optimize")
    return(df)
  }
  
  # filter out non-optimized parameters
  filter_var <- c()
  if(opt_xt) filter_var <- c(filter_var,"'xt'")
  if(opt_a) filter_var <- c(filter_var,"'a'")
  exp <- parse(text=paste0("type==",filter_var,collapse = " | "))
  df <- dplyr::filter(df,eval(exp))
  
  if(nrow(df)==0){
    if(warn_when_none) warning("No parameters to optimize")
    return(df)
  }
  
  if(transform_parameters){
    df <- df %>% dplyr::rowwise() %>% 
      dplyr::mutate(par=dplyr::if_else(cont==TRUE,
                                   # FastImputation::NormalizeBoundedVariable(par,
                                   #                                          constraints =
                                   #                                            list(lower=lower,
                                   #                                                 upper=upper)),
                                   unbound_par(par,lower=lower,upper=upper),
                                   par))
    
    df <- df %>% dplyr::mutate(transformed=dplyr::if_else(cont==TRUE,TRUE,FALSE),
                               lower_orig=lower, upper_orig=upper)
    df <- df %>% dplyr::mutate(
                               upper=dplyr::if_else(cont==TRUE,Inf,upper),
                               lower=dplyr::if_else(cont==TRUE,-Inf,lower)
                               )
    
    #     if(is.finite(lower[i]) && is.finite(upper[i])){
    #       par[i] <- (par[i] - lower[i])/(upper[i]-lower[i])
    #       par[i] <- stats::qlogis(par[i])
    #     }
    #
  } else {
    df <- df %>% dplyr::mutate(transformed=FALSE,
                               lower_orig=lower, 
                               upper_orig=upper)
  }
    
  
  attr(df,"opt_xt") <- ifelse(opt_xt,TRUE,FALSE)
  attr(df,"opt_a") <- ifelse(opt_a,TRUE,FALSE)
  attr(df,"opt_samps") <- ifelse(opt_samps,TRUE,FALSE)
  attr(df,"opt_inds") <- ifelse(opt_inds,TRUE,FALSE)
  attr(df,"transformed") <- ifelse(any(df$transformed),TRUE,FALSE)
  
  return(df)
  
}

# put_par_optim <- function(tbl_opt,tbl_full) {
#   
#   transformed <- par <- lower_orig <- upper_orig <- type <- par_opt <- NULL
#   
#   #transformed back to normal scale
#   if("transformed" %in% names(tbl_opt)){
#     tbl_opt <- 
#       dplyr::mutate(tbl_opt,
#                     par=
#                       dplyr::if_else(
#                         transformed==TRUE,
#                         # FastImputation::BoundNormalizedVariable(
#                         #   par,
#                         #   constraints =
#                         #     list(lower=lower_orig,
#                         #          upper=upper_orig)
#                         # ),
#                         bound_par(par,lower=lower_orig,upper=upper_orig),
#                         par
#                       )
#       )
#   }
#     
#   # tbl_opt_small <- tbl_opt %>% dplyr::select(par,type,grouping) %>% dplyr::rename(par_opt=par)
#   tbl_opt_small <- dplyr::select(tbl_opt, par,type,grouping) 
#   tbl_opt_small <-  dplyr::rename(tbl_opt_small, par_opt=par)
# 
#   # tbl_full <- dplyr::left_join(tbl_full,tbl_opt_small,by=c("grouping","type")) %>%
#     # dplyr::mutate(par=dplyr::if_else(is.na(par_opt),par,par_opt)) %>%
#     # dplyr::select(-par_opt)
#   tbl_full <- dplyr::left_join(tbl_full,tbl_opt_small,by=c("grouping","type"))
#   tbl_full <- dplyr::mutate(tbl_full,par=dplyr::if_else(is.na(par_opt),par,par_opt))
#   tbl_full <- dplyr::select(tbl_full,-par_opt)
#   
#   return(tbl_full)
# }

# 
# get_xt_from_tbl <- function(tbl_full,add_names=FALSE){
#   ret <- tbl_full %>% dplyr::filter(type=="xt") %>% dplyr::select(par,group,name) %>% 
#     tidyr::spread(key=name,value = par) #%>% 
#   if(add_names){ 
#     ret <- ret %>% dplyr::mutate(group=paste0("grp_",group)) %>% as.data.frame() %>% 
#       tibble::column_to_rownames(var="group") %>% data.matrix()  
#   } else {
#     ret <- ret %>% dplyr::select(-group)
#   }
#   return(ret)
# }

# get_a_from_tbl <- function(tbl_full){
#   tbl_full %>% dplyr::filter(type=="a") %>% dplyr::select(par,group,name) %>% 
#     tidyr::spread(key=name,value = par) %>% 
#     #dplyr::select(-group)
#     dplyr::mutate(group=paste0("grp_",group)) %>% as.data.frame() %>% 
#     tibble::column_to_rownames(var="group") %>% data.matrix()  
# }

# get_type_from_tbl <- function(type_str,tbl_full,add_row_names=TRUE){
#   type <- par <- group <- name <- NULL
#   ret <- dplyr::filter(tbl_full,type==!!type_str) 
#   ret <- dplyr::select(ret,par,group,name) 
#   ret <- tidyr::spread(ret,key=name,value = par) 
#   if(add_row_names){ 
#     ret <- dplyr::mutate(ret,group=paste0("grp_",group)) 
#     ret <- data.matrix(tibble::column_to_rownames(as.data.frame(ret),var="group"))
#   } else {
#     ret <- dplyr::select(ret,-group)
#   }
#   return(ret)
# }



