
#' Display a summary of output from poped_optim
#'
#' @param object An object returned from \code{\link{poped_optim}} to summarize. 
# @param fn A file handle to write to.  Default is to the R console.
#' @param ... Additional arguments. Passed to \code{\link{blockfinal}}.
#'
#' @return NULL
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_summary.poped_optim.R
#' @export

summary.poped_optim <- function(object,
                              #fn="",
                              ...){
  
  # fn=blockheader(object$initial$poped.db,
  #                name=fn,
  #                #e_flag=!object$poped.db$d_switch,
  #                fmf=object$initial$FIM,
  #                dmf=object$initial$ofv,
  #                #out_file=out_file,
  #                #trflag=trace,
  #                ...)
  # 
  # blockexp(fn=fn,object$initial$poped.db,
  #          e_flag=!(object$initial$poped.db$settings$d_switch),...)
  
  #object <- x
  
  blockfinal(fn="",fmf=object$FIM,
             dmf=object$ofv,
             groupsize=object$poped.db$design$groupsize,
             ni=object$poped.db$design$ni,
             xt=object$poped.db$design$xt,
             x=object$poped.db$design$x,
             a=object$poped.db$design$a,
             model_switch=object$poped.db$design$model_switch,
             object$poped.db$parameters$param.pt.val$bpop,
             object$poped.db$parameters$param.pt.val$d,
             object$poped.db$parameters$docc,
             object$poped.db$parameters$param.pt.val$sigma,
             object$poped.db,
             fmf_init=object$initial$FIM,
             dmf_init=object$initial$ofv,
             run_time = object$run_time,
             ...)
  #NextMethod("print")
  return(invisible())
  
}