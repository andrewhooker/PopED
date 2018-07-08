add_to_poped_db <- function(
  poped_db,
  ps_tbl,
  par=ps_tbl$par,
  ps_transformed=ps_tbl$transformed
){
  
  # TODO: matches with ofv_optim should join this code  
  # variables needed
  xt <- poped_db$design$xt
  a <- poped_db$design$a
  ps_type <- ps_tbl$type
  ps_index <- ps_tbl$index
  ps_lower_orig=ps_tbl$lower_orig
  ps_upper_orig=ps_tbl$upper_orig

  # transform parameters if needed
  par <- transform_back_par(ps_tbl,
                            par=par,
                            ps_transformed=ps_transformed,
                            ps_lower_orig=ps_lower_orig,
                            ps_upper_orig=ps_upper_orig)
  
  # put xt parameters in right places
  if(attr(ps_tbl,"opt_xt")==TRUE){
    par_xt <- par[ps_type=="xt"]
    par_xt_index <- ps_index[ps_type=="xt"]
    xt_tmp <- t(xt)
    for(i in 1:length(par_xt)) xt_tmp[par_xt_index[i][[1]]] = par_xt[i]
    xt <- t(xt_tmp)
    poped_db$design$xt[,]=xt[,]
  }
  
  # put a parameters in right places
  if(attr(ps_tbl,"opt_a")==TRUE){
    par_a <- par[ps_type=="a"]
    par_a_index <- ps_index[ps_type=="a"]
    a_tmp <- t(a)
    for(i in 1:length(par_a)) a_tmp[par_a_index[i][[1]]] = par_a[i]
    a <- t(a_tmp)
    poped_db$design$a[,]=a[,]
  }
  
  return(poped_db)
}