#---------- functions
dots <- function(...) {
  eval(substitute(alist(...)))
}

transform_back <- function(par,lower=-Inf,upper=Inf){
  FastImputation::BoundNormalizedVariable(
    par,
    constraints = 
      list(lower=lower,
           upper=upper))
}

transform_back_par <- function(
  ps_tbl,
  par=ps_tbl$par,
  ps_transformed=ps_tbl$transformed,
  ps_lower_orig=ps_tbl$lower_orig,
  ps_upper_orig=ps_tbl$upper_orig)
{
  if(any(ps_transformed))
    par[ps_transformed] <- 
      mapply(transform_back,par[ps_transformed],
             ps_lower_orig[ps_transformed],
             ps_upper_orig[ps_transformed])
  return(par)
}