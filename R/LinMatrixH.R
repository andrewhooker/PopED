#' Model linearization with respect to epsilon.
#' 
#' The function performs a linearization of the model with respect to the residual variability.
#' Derivative of model w.r.t. eps evaluated at eps=0
#' 
#' @inheritParams mftot
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param xt_ind A vector of the individual/group sample times
#' @param b_ind vector of individual realization of the BSV terms b
#' @param bocc_ind Vector of individual realizations of the BOV terms bocc
#' 
#' @return A matrix of size (samples per individual x number of epsilons) 
#'  
#' @family FIM
#'     
# @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
# @example tests/testthat/examples_fcn_doc/examples_LinMatrixH.R
# @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixH <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db){
  #----------Model linearization with respect to epsilon.
  #
  # size of return is (samples per individual x number of epsilons) 
  #
  # derivative of model w$r.t. eps eval at e=0
  #
  NumEPS = size(poped.db$parameters$sigma,1)
  if((NumEPS==0)){
    y=0
  } else { 
    returnArgs <- gradf_eps(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,NumEPS,poped.db) 
    y <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
  }
  return(list( y= y,poped.db=poped.db)) 
}

#' Model linearization with respect to epsilon.
#' 
#' The function performs a linearization of the model with respect to the residual variability.
#' Derivative of model w.r.t. eps evaluated at eps=0 and b=b_ind.
#' 
#' @inheritParams mftot
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams LinMatrixH
#' @param num_eps The number of \code{eps()} in the model.
#' 
#' @return A matrix of size (samples per individual x number of epsilons) 
#'  
#' @family FIM
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_gradf_eps.R
#' @export
#' @keywords internal
#' 
gradf_eps <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,num_eps,poped.db){
  #----------Model linearization with respect to epsilon.
  #
  # size of return is (samples per individual x number of epsilons)
  #
  # derivative of model w$r.t. eps eval at e=0 and b=b_ind
  #
  #
  
  if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==1) ){#No interaction
    fg0=feval(poped.db$model$fg_pointer,x,a,bpop,zeros(size(b_ind)[1],size(b_ind)[2]),zeros(size(bocc_ind)[1],size(bocc_ind)[2]))
    
  } else {
    fg0=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind) #Interaction
  }
  
  e0=zeros(1,num_eps)
  dfeps_de0 = grad_all(poped.db$model$ferror_pointer,4,size(xt_ind,1),model_switch,xt_ind,fg0,e0,poped.db)
  
  return(list( dfeps_de0= dfeps_de0,poped.db=poped.db)) 
}
