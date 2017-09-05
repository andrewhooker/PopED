#------- create optimization function with optimization parameters first
ofv_optim <- function(par,ps_tbl_opt,ps_tbl_full,poped_db,...){
  
  # add par to ps_tbl_opt
  ps_tbl_opt$par <- par
  
  # add ps_tbl_opt to ps_tbl_full
  ps_tbl_full <- put_par_optim(ps_tbl_opt,ps_tbl_full)
  
  xt <- NULL
  if(ps_tbl_opt %>% dplyr::filter(type=="xt") %>% nrow()>0) xt <-  get_xt_from_tbl(ps_tbl_full)
  
  a <- NULL
  if(ps_tbl_opt %>% dplyr::filter(type=="a") %>% nrow()>0) a <- get_a_from_tbl(ps_tbl_full)
  
  extra_args <- dots(...)
  extra_args$evaluate_fim <- FALSE
  
  output <- do.call(calc_ofv_and_fim,
                    c(list(
                      poped.db,
                      xt=xt,
                      a=a
                    ),
                    extra_args))
  
  ofv <- output$ofv
  
  if(!is.finite(ofv)) ofv <- NA
  
  return(ofv)
}