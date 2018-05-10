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