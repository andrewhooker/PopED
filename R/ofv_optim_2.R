#------- create optimization function with optimization parameters first
ofv_optim_2 <- function(par,ps_tbl,poped_db,...){
  
  browser()
  library(microbenchmark)
  compare <- microbenchmark(
    dplyr::mutate(ps_tbl,
                  par=
                    dplyr::if_else(
                      transformed==TRUE,
                      FastImputation::BoundNormalizedVariable(
                        par,
                        constraints =
                          list(lower=lower_orig,
                               upper=upper_orig)
                      ),
                      par
                    )
    ) 
    ,mapply(transform_back,par,ps_tbl$lower_orig,ps_tbl$upper_orig) 
    ,times = 100L
  )
  compare
  library(ggplot2)
  autoplot(compare)
  
  # transform parameters
  par[ps_tbl$transformed] <- 
    mapply(transform_back,par,
           ps_tbl$lower_orig[ps_tbl$transformed],
           ps_tbl$upper_orig[ps_tbl$transformed])
  
  # put parameters in right places
  if(attr(ps_tbl,"opt_xt")==TRUE){
    par_xt <- par[ps_tbl$type=="xt"]
  }
  
  xt <- NULL
  if(any(ps_tbl_opt$type=="xt")) xt <-  get_type_from_tbl("xt",ps_tbl_full, add_row_names = FALSE)
  #if(ps_tbl_opt %>% dplyr::filter(type=="xt") %>% nrow()>0) xt <-  get_xt_from_tbl(ps_tbl_full)
  
  a <- NULL
  if(any(ps_tbl_opt$type=="a")) a <-  get_type_from_tbl("a",ps_tbl_full, add_row_names = FALSE)
  #if(ps_tbl_opt %>% dplyr::filter(type=="a") %>% nrow()>0) a <- get_a_from_tbl(ps_tbl_full)
  
  #extra_args <- dots(...)
  #extra_args$evaluate_fim <- FALSE
  
  output <- calc_ofv_and_fim(poped.db,xt=xt,a=a,evaluate_fim = FALSE,...)
  # output <- do.call(calc_ofv_and_fim,
  #                   c(list(
  #                     poped.db,
  #                     xt=xt,
  #                     a=a
  #                   ),
  #                   extra_args))
  
  ofv <- output$ofv
  
  if(!is.finite(ofv)) ofv <- NA
  
  return(ofv)
}
