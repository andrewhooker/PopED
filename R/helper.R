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