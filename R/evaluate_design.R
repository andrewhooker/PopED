#' Evaluate a design
#' 
#' This function evaluates the design defined in a poped database.
#' 
#' @param poped.db A poped database
#' @param ... Extra parameters passed to \code{\link{calc_ofv_and_fim}} and \code{\link{get_rse}}
#' @return A list of elements evaluating the current design.
#' @export
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_evaluate_design.R
#' @family evaluate_design

evaluate_design <- function(poped.db, ...) {
  out <- calc_ofv_and_fim(poped.db,...)
  if(is.null(out$fim)){
    out$rse <- NULL
  } else{
    out$rse <- get_rse(out$fim,poped.db,...)
  }
  return(out)
}