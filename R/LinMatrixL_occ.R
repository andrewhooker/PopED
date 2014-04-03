#' Model linearization with respect to occasion variablity parameters.
#' 
#' The function performs a linearization of the model with respect to the occation  variability parameter..
#' Derivative of model w.r.t. eta_occ, evaluated bocc_ind.
#' 
#' @inheritParams mftot
#' @inheritParams LinMatrixH
#' @param iCurrentOcc The current occasion.
#' 
#' @return A matrix of size (samples per individual x number of iovs)
#'  
#' @family FIM
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixL_occ <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,iCurrentOcc,globalStructure){
#
# size: (samples per individual x number of iovs)
#
if((globalStructure$NumOcc==0)){
	y=0
} else {
     returnArgs <- gradff(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure) 
grad_ff_tmp <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    y=grad_ff_tmp%*%gradfg_occ(x,a,bpop,b_ind,bocc_ind,iCurrentOcc,globalStructure)
}
return(list( y= y,globalStructure=globalStructure)) 
}


