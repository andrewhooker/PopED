#' Evaluate the Fisher Information Matrix (FIM)
#' 
#' Compute the FIM given specific model(s), parameters, design and methods. 
#' 
#' @param globalStructure A PopED database.
#' @param bpop The fixed effects parameter values.  Supplied as a vector.
#' @param d A between subject variability matrix (OMEGA in NONMEM).
#' @param docc A between occasion variability matrix.
#' @param sigma A residual unexplained variability matrix (SIGMA in NONMEM).
#' @param model_switch A matrix that is the same size as xt, specifying which model each sample belongs to.
#' @param ni A vector of the number of samples in each group.
#' @param xt A matrix of sample times.  Each row is a vector of sample times for a group.
#' @param x A matrix for the discrete design variables.  Each row is a group.
#' @param a A matrix of covariates.  Each row is a group.
#' @param groupsize A vector of the numer of individuals in each group.
#' 
#' @return As a list:
#' \item{ret}{The FIM}
#' \item{globalStructure}{A PopED database}
#' 
#' @seealso For an easier function to use, please see \code{\link{evaluate.fim}}.  
#' @family FIM
#' 


## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

mftot <- function(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure){
    
    returnArgs <- switch(globalStructure$iFIMCalculationType+1,
                         mftot0(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure), 
                         mftot1(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure), 
                         mftot2(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure),
                         mftot3(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure), 
                         mftot4(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure),
                         mftot5(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure),
                         mftot6(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure),
                         mftot7(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure)) 
    
    if(is.null(returnArgs)) stop(sprintf('Unknown FIM-calculation type'))
    ret <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    return(list( ret= ret,globalStructure =globalStructure )) 
}
