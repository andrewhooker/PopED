
#' Get parameter list to optimize and the space of those parameters.
#'
#' @param poped.db A poped database
#' @param opt_xt Optimize on time values?
#' @param opt_a Optimize on non-time values?
#' @param opt_samps Optimize on number of samples?
#' @param opt_inds Optimize on number of individuals?
#' @param transform_parameters Transform parameters that have bounds to be boundless before optimization?
#' @param ... Arguments passed to other functions.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' get_par_and_space(poped.db,opt_xt=T)
#' get_par_and_space(poped.db,opt_xt = T,opt_a = T)
#' get_par_and_space(poped.db,opt_a = T)
#' 
get_par_and_space <- function(poped.db,
                              opt_xt=as.logical(poped.db$settings$optsw[2]),
                              opt_a=as.logical(poped.db$settings$optsw[4]),
                              opt_samps=as.logical(poped.db$settings$optsw[1]),
                              opt_inds=as.logical(poped.db$settings$optsw[5]),
                              #transform_parameters=T,
                              ...) {
  
  #----------- checks
  if(!any(opt_xt,opt_a,opt_samps,opt_inds)){
    stop('No optimization parameter is set.')
  }
  
  df_xt <- NULL
  df_a <- NULL
  df <- NULL
  
  design=poped.db$design
  design_space=poped.db$design_space
  
  # Collect the parameters to optimize
  if(opt_samps) stop('Sample number optimization is not yet implemented in the R-version of PopED.')
  if(opt_inds) stop('Optimization  of number of individuals in different groups is not yet implemented in the R-version of PopED.')
  if(opt_xt){ 
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
      df_xt <- tibble::add_column(df_xt,allowed_values=NA)
    } else{
      df_xt <- tibble::add_column(df_xt,allowed_values=t(design_space$xt_space)[t(sel_mat_xt)])
    }
    
    # group column
    group_mat <- ones(size(design$xt))*c(1:nrow(design$xt))
    df_xt <- tibble::add_column(df_xt,group=t(group_mat)[t(sel_mat_xt)])
    
    # model_switch column
    df_xt <- tibble::add_column(df_xt,model=t(design$model_switch)[t(sel_mat_xt)])

    df <- dplyr::bind_rows(df,df_xt)
  }
  if(opt_a) { 
    
    sel_mat_a <- ones(size(design$a))
    for(i in 1:length(design$ni)){
      if((design$ni[i]!=0 && design$groupsize[i]!=0)) sel_mat_a[i,] <- 1
    }
    sel_mat_a <- sel_mat_a==1
    
    #grouping_max <- 0
    #if(!is.null(df)) grouping_max <- max(df$grouping)

    df_a <- tibble::tibble(par=t(design$a)[t(sel_mat_a)])
    df_a <- tibble::add_column(df_a,lower=t(design_space$mina)[t(sel_mat_a)])
    df_a <- tibble::add_column(df_a,upper=t(design_space$maxa)[t(sel_mat_a)])
    #df_a <- tibble::add_column(df_a,grouping=grouping_max+t(design_space$G_a)[t(sel_mat_a)])
    df_a <- tibble::add_column(df_a,grouping=t(design_space$G_a)[t(sel_mat_a)])
    
    df_a <- tibble::add_column(df_a,type="a")
    if(is.null(design_space$a_space)){
      df_a <- tibble::add_column(df_a,allowed_values=NA)
    } else{
      df_a <- tibble::add_column(df_a,allowed_values=t(design_space$a_space)[t(sel_mat_a)])
    }
    
    # group column
    group_mat <- ones(size(design$a))*c(1:nrow(design$a))
    df_a <- tibble::add_column(df_a,group=t(group_mat)[t(sel_mat_a)])
    
    df <- dplyr::bind_rows(df,df_a)
    
  }

  # continuous and discrete parameters
  df <- dplyr::mutate(df,allowed_values_2=ifelse(is.na(allowed_values) | length(allowed_values)==0,NA,allowed_values))
  df <- dplyr::mutate(df,allowed_values=NULL)
  df <- dplyr::rename(df,allowed_values=allowed_values_2)
  df <- dplyr::mutate(df,cat_cont=ifelse(is.na(allowed_values),"cont","cat"))
  
  
  # Parameter grouping
  df <- dplyr::distinct(df,grouping,type,.keep_all=TRUE) 

  
  # find and filter fixed parameters 
  df <- dplyr::mutate(df,fixed=ifelse(lower==upper | 
                                        (length(allowed_values==1) 
                                         & !is.na(allowed_values[[1]]) 
                                         & allowed_values[[1]]==par),
                                      TRUE,FALSE)) 
  
  df <- dplyr::filter(df,fixed==FALSE)
  if(nrow(df)==0) stop("No design parameters have a design space to optimize")
  
  
  # # update allowed values to be within limits
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
  
  
  # if(transform_parameters){
  #   for(i in 1:length(par)){
  #     if(par_cat_cont[i]=="cont"){
  #       par[i] <- FastImputation::NormalizeBoundedVariable(par[i],
  #                                                          constraints = 
  #                                                            list(lower=lower[i],
  #                                                                 upper=upper[i]))
  #       
  #       #     if(is.finite(lower[i]) && is.finite(upper[i])){
  #       #       par[i] <- (par[i] - lower[i])/(upper[i]-lower[i])
  #       #       par[i] <- stats::qlogis(par[i])
  #       #     }
  #       #     
  #     }
  #   }
  # }
  return(df)
}


