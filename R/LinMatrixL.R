#' The linearized matrix L
#' 
#' Function computes the derivative of the model with respect to the between subject variability terms in the model (b's and bocc's) evaluated at
#' a defined point 
#' (b_ind and bocc_ind).
#' 
#' @inheritParams mf
#' @param bpop The fixed effects parameter values.  Supplied as a vector.
#' @param b_ind The point at which to evaluate the derivative
#' @param bocc_ind The point at which to evaluate the derivative
#' @param globalStructure A PopED database.
#' 
#' @return As a list:
#' \item{y}{A matrix of size (samples per individual x number of random effects)}
#' \item{globalStructure}{A PopED database}
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixL <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure){
  
  if((globalStructure$NumRanEff==0)){
    y=0
  } else {
    returnArgs <- gradff(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure) 
    grad_ff_tmp <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    y=grad_ff_tmp%*%gradfg(x,a,bpop,b_ind,bocc_ind,globalStructure)
  }
  return(list( y= y,globalStructure=globalStructure)) 
}


