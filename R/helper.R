#---------- functions
dots <- function(...) {
  eval(substitute(alist(...)))
}

unbound_par <- function(
  x, lower=-Inf, upper=Inf,tol = 1e-8,
  logit=T)
{
  
  if (upper < lower) stop("'upper' must be greater than 'lower.'")
  if (all(is.finite(c(lower,upper))) & (tol > (upper - lower)/2)) {
    stop("'tol' must be less than half the distance between the upper and lower bounds.")
  }
  x <- pmax(x, lower + tol)
  x <- pmin(x, upper - tol)
  
  if(all(is.infinite(c(lower,upper)))) {
    return(x)
  } else if(is.infinite(lower)) {
    return(log(upper - x))
  } else if(is.infinite(upper)) {
    return(log(x - lower))
  } else if(logit){
    return(stats::qlogis((x - lower)/(upper-lower)))
  } else{ # probit
    return(stats::qnorm((x - lower)/(upper-lower)))
  } 
}

bound_par <- function(
  x, lower=-Inf, upper=Inf,
  logit=T)
{
  
  if (upper < lower) stop("'upper' must be greater than 'lower.'")
  
  if(all(is.infinite(c(lower,upper)))) {
    return(x)
  } else if(is.infinite(lower)) {
    return(upper - exp(x))
  } else if(is.infinite(upper)) {
    return(exp(x) + lower)
  } else if(logit){
    return(stats::plogis(x)*(upper - lower) + lower)
  } else{ # probit
    return(stats::pnorm(x)*(upper - lower) + lower)
  } 
}


# transform_back <- function(par,lower=-Inf,upper=Inf){
#   # FastImputation::BoundNormalizedVariable(
#   #   par,
#   #   constraints = 
#   #     list(lower=lower,
#   #          upper=upper))
#   bound_par(par,lower=lower,upper=upper)
# }

transform_back_par <- function(
  ps_tbl,
  par=ps_tbl$par,
  ps_transformed=ps_tbl$transformed,
  ps_lower_orig=ps_tbl$lower_orig,
  ps_upper_orig=ps_tbl$upper_orig)
{
  if(any(ps_transformed))
    par[ps_transformed] <- 
      mapply(bound_par,par[ps_transformed],
             ps_lower_orig[ps_transformed],
             ps_upper_orig[ps_transformed])
  return(par)
}


##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings (with value) and errors
##' @param expr an \R expression to evaluate
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler, The R Core Team
##' @keywords internal
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


get_fim_size <- function(poped.db) {
  numnotfixed_bpop = sum(poped.db$parameters$notfixed_bpop)
  numnotfixed_d    = sum(poped.db$parameters$notfixed_d)
  numnotfixed_covd = sum(poped.db$parameters$notfixed_covd)
  numnotfixed_docc  = sum(poped.db$parameters$notfixed_docc)
  numnotfixed_covdocc  = sum(poped.db$parameters$notfixed_covdocc)
  numnotfixed_sigma  = sum(poped.db$parameters$notfixed_sigma)
  numnotfixed_covsigma  = sum(poped.db$parameters$notfixed_covsigma)
  
  n_fixed_eff <- numnotfixed_bpop
  n_rand_eff <- numnotfixed_d+numnotfixed_covd+numnotfixed_docc+
    numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma
  fim_size <- n_fixed_eff+n_rand_eff
  return(fim_size)
}