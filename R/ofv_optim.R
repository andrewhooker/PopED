#------- create optimization function with optimization parameters first
ofv_optim <- function(par,ps_tbl,poped_db,...){
  
  # TODO: matches with add_to_poped_db should join this code  
  
  # variables needed
  xt <- poped_db$design$xt
  a <- poped_db$design$a
  ps_type <- ps_tbl$type
  ps_index <- ps_tbl$index
  ps_transformed <- ps_tbl$transformed
  ps_lower_orig <- ps_tbl$lower_orig
  ps_upper_orig <- ps_tbl$upper_orig
  
  # transform parameters
  if(attr(ps_tbl,"transformed")==TRUE)
    par[ps_transformed] <-
    mapply(bound_par,par[ps_transformed],
           ps_lower_orig[ps_transformed],
           ps_upper_orig[ps_transformed])

  
  # put xt parameters in right places
  if(attr(ps_tbl,"opt_xt")==TRUE){
    par_xt <- par[ps_type=="xt"]
    par_xt_index <- ps_index[ps_type=="xt"]
    xt_tmp <- t(xt)
    for(i in 1:length(par_xt)) xt_tmp[par_xt_index[i][[1]]] = par_xt[i]
    xt <- t(xt_tmp)
  }
  
  # put a parameters in right places
  if(attr(ps_tbl,"opt_a")==TRUE){
    par_a <- par[ps_type=="a"]
    par_a_index <- ps_index[ps_type=="a"]
    a_tmp <- t(a)
    for(i in 1:length(par_a)) a_tmp[par_a_index[i][[1]]] = par_a[i]
    a <- t(a_tmp)
  }

  # output <- calc_ofv_and_fim(poped_db,xt=xt,a=a,...,evaluate_fim = FALSE)
  # extra_args <- dots(...)
  extra_args <- list(...) 
  extra_args$evaluate_fim <- FALSE
  output <- do.call(calc_ofv_and_fim,
                    c(list(
                      poped_db,
                      xt=xt,
                      a=a
                    ),
                    extra_args))
  
  ofv <- output$ofv
  
  if(!is.finite(ofv)) ofv <- NA
  
  return(ofv)
}
