#' Power of a design to estimate a parameter. 
#' 
#' Evaluate the power of a design to estimate a parameter value different than
#' some assumed value (often the assumed value is zero). The power is calculated 
#' using the linear Wald test and the the design is defined in a poped database.
#'
#' @param poped.db A poped database
#' @param bpop_idx Index for an unfixed population parameter (bpop) for 
#'   which the power should be
#'   evaluated for being different than the null hypothesis (h0).
#' @param alpha Type 1 error.
#' @param h0 The null hypothesized value for the parameter.
#' @param power Targeted power.
#' @param twoSided Is this a two-sided test.
#' @param fim Provide the FIM from a previous calculation
#' @param out provide output from a previous calculation (e.g.,
#'   calc_ofv_and_fim, ...)
#' @param ... Extra parameters passed to \code{\link{calc_ofv_and_fim}} and
#'   \code{\link{get_rse}}
#' @param find_min_n Should the function compute the minimum n needed (given the
#'   current design) to achieve the desired power?
#' @return A list of elements evaluating the current design including the power.
#' @references \enumerate{ \item Retout, S., Comets, E., Samson, A., and Mentre,
#'   F. (2007). Design in nonlinear mixed effects models: Optimization using the
#'   Fedorov-Wynn algorithm and power of the Wald test for binary covariates.
#'   Statistics in Medicine, 26(28), 5162-5179.
#'   \doi{10.1002/sim.2910}. \item Ueckert, S., Hennig, S.,
#'   Nyberg, J., Karlsson, M. O., and Hooker, A. C. (2013). Optimizing disease
#'   progression study designs for drug effect discrimination. Journal of
#'   Pharmacokinetics and Pharmacodynamics, 40(5), 587-596.
#'   \doi{10.1007/s10928-013-9331-3}. }
#'
#' @example tests/testthat/examples_fcn_doc/examples_evaluate_power.R
#'
#' @family evaluate_design
#' @export

evaluate_power <- function(poped.db, bpop_idx, h0=0, alpha=0.05, power=0.80, twoSided=TRUE, 
                           find_min_n=TRUE,
                           fim=NULL, out=NULL,...) {
  # If two-sided then halve the alpha
  if (twoSided == TRUE) alpha = alpha/2
  
  # Check if bpop_idx is given and within the non-fixed parameters
  if (!all(bpop_idx %in% which(poped.db$parameters$notfixed_bpop==1))) 
    stop("bpop_idx can only include non-fixed population parameters bpop")
  if (poped.db$parameters$param.pt.val$bpop[bpop_idx]==0) 
    stop("
  Population parameter is assumed to be zero, 
  there is 0% power in identifying this parameter 
  as non-zero assuming no bias in parameter estimation")
  
  # Prepare output structure with at least out$fim available
  if (is.null(fim) & is.null(out$fim)) {
    out <- calc_ofv_and_fim(poped.db,...)
  } else if (!is.null(fim)) {
    out = list(fim = fim)
  }
  # Add out$rse
  out$rse <- get_rse(out$fim,poped.db,...) # in percent!!

  # Derive power and RSE needed for the selected parameter(s)
  norm.val = abs(qnorm(alpha, mean=0, sd=1))
  val = poped.db$parameters$param.pt.val$bpop[bpop_idx]
  rse = out$rse[cumsum(poped.db$parameters$notfixed_bpop==1)[bpop_idx]]/100 
  se <- rse*val
  wald_stat <- (h0-val)/se
  
  # Following the paper of Retout et al., 2007 for the Wald-test:
  #powPred = round(100*(1 - stats::pnorm(norm.val-(100/rse)) + stats::pnorm(-norm.val-(100/rse))), digits=1)
  power_pred = 1 - stats::pnorm(wald_stat+norm.val) + stats::pnorm(wald_stat-norm.val)
  
  #needRSE = 100/(norm.val-stats::qnorm(1-power/100))
  need_se = abs(h0-val)/(norm.val-stats::qnorm(1-power))
  need_rse <- need_se/val

  #out$power = data.frame(Value=val, RSE=rse, predPower=powPred, wantPower=power, needRSE=needRSE)
  out$power = data.frame(Value=val, RSE=rse*100, power_pred=power_pred*100, power_want=power*100, need_rse=need_rse*100)
  
  # find the smallest n to achieve the wanted power.
  #if(find_min_n & nrow(poped.db$settings$prior_fim)==0){
  if(find_min_n){
    res <- optimize_n_rse(poped.db,bpop_idx=bpop_idx,need_rse=need_rse,use_percent = FALSE)
    out$power$min_N_tot=res["n"]
  }
  
  return(out)
}