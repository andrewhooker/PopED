#' Evaluate power of a design for detecting a parameter to be non-zero
#' 
#' This tunction evaluates the design defined in a poped database.
#' 
#' @param poped.db A poped database
#' @param bpopIdx Indices for unfixed parameters for which power should be evaluated for being non-zero
#' @param alpha Type 1 error (default = 0.05)
#' @param power Targeted power (default = 80%)
#' @param twoSided Is two-sided test (default = TRUE)
#' @param fim Optional to provide FIM from a previous calculation
#' @param out Optional to provide output from a previous calculation (e.g., calc_ofv_and_fim, ...)
#' @param ... Extra parameters passed to \code{\link{calc_ofv_and_fim}} and \code{\link{get_rse}}
#' @return A list of elements evaluating the current design including the power.
#' @export
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_evaluate_design.R
#' @family evaluate_design

evaluate_power <- function(poped.db, bpopIdx=NULL, fim=NULL, out=NULL, alpha=0.05, power=80, twoSided=TRUE, ...) {
  # If two-sided then halve the alpha
  if (twSided == TRUE) alpha = alpha/2
  
  # Check if bpopIdx is given and within the non-fixed parameters
  if (is.null(bpopIdx)) stop("Population parameter index must be given in bpopIdx")
  if (!all(bpopIdx %in% which(poped.db$parameters$notfixed_bpop==1))) stop("bpopIdx can only include non-fixed population parameters bpop")
  
  # Prepare output structure with at least out$fim available
  if (is.null(fim) & is.null(out$fim)) {
    out <- calc_ofv_and_fim(poped.db,...)
  } else if (!is.null(fim)) {
    out = list(fim = fim)
  }
  # Add out$rse
  out$rse <- get_rse(out$fim,poped.db,...)

  # Derive power and RSE needed for the selected parameter(s)
  norm.val = abs(qnorm(alpha, mean=0, sd=1))
  val = poped.db$parameters$param.pt.val$bpop[bpopIdx]
  rse = out$rse[which(poped.db$parameters$notfixed_bpop==1)[bpopIdx]] # in percent!!

  # Following the paper of Retout et al., 2007 for the Wald-test:
  powPred = round(100*(1 - pnorm(norm.val-(100/rse)) + pnorm(-norm.val-(100/rse))), digits=1)
  needRSE = 100/(norm.val-qnorm(1-power/100))

  out$power = data.frame(Value=val, RSE=rse, predPower=powPred, wantPower=power, needRSE=needRSE)
  return(out)
}