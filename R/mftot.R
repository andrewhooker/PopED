#' Evaluate the Fisher Information Matrix (FIM)
#' 
#' Compute the FIM given specific model(s), parameters, design and methods. 
#' 
#' @param poped.db A PopED database.
#' @param bpop The fixed effects parameter values.  Supplied as a vector.
#' @param d A between subject variability matrix (OMEGA in NONMEM).
#' @param docc A between occasion variability matrix.
#' @param sigma A residual unexplained variability matrix (SIGMA in NONMEM).
#' @param model_switch A matrix that is the same size as xt, specifying which model each sample belongs to.
#' @param ni A vector of the number of samples in each group.
#' @param xt A matrix of sample times.  Each row is a vector of sample times for a group.
#' @param x A matrix for the discrete design variables.  Each row is a group.
#' @param a A matrix of covariates.  Each row is a group.
#' @param groupsize A vector of the number of individuals in each group.
#' 
#' @return As a list:
#' \item{ret}{The FIM}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso For an easier function to use, please see \code{\link{evaluate.fim}}.  
#' @family FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_mftot.R
#' @export
#' @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

mftot <- function(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db,...){
  m=size(ni,1)
  s=0
  for(i in 1:m){
    if((ni[i]!=0 && groupsize[i]!=0)){
      if((!isempty(x))){
        x_i = t(x[i,,drop=F])      
      } else {
        x_i =  zeros(0,1)
      }
      if((!isempty(a))){
        a_i = t(a[i,,drop=F])
      } else {
        a_i =  zeros(0,1)
      }
      # mf_all <- function(model_switch,xt,x,a,bpop,d,sigma,docc,poped.db){
      extra_args <- list(...)
      if(is.null(extra_args$loq) & is.null(extra_args$uloq)){ # no loq
        returnArgs <- mf_all(t(model_switch[i,1:ni[i,drop=F],drop=F]),
                             t(xt[i,1:ni[i,drop=F],drop=F]),
                             x_i,a_i,bpop,d,sigma,docc,poped.db) 
      } else { # handle LOQ 
        returnArgs <- mf_all_loq(t(model_switch[i,1:ni[i,drop=F],drop=F]),
                             t(xt[i,1:ni[i,drop=F],drop=F]),
                             x_i,a_i,bpop,d,sigma,docc,poped.db,...) 
      }
      if(is.null(returnArgs)) stop(sprintf('Unknown FIM-calculation type'))
      mf_tmp <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      s=s+groupsize[i]*mf_tmp
    }
  }

  return(list( ret= s,poped.db =poped.db )) 
}
