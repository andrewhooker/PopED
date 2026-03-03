gradfg_occ <- function(x,a,bpop,b_ind,bocc_ind,currentOcc,poped.db){
#
#
# size: (number of g's x NumOccVariables)
#
# deriv of g's w$r.t. bocc's eval at b_ind, bocc_ind at occ currentOcc
#
  dfg_db0 = grad_all(poped.db$model$fg_pointer,5,poped.db$parameters$ng,x,a,bpop,b_ind,bocc_ind,poped.db,currentOcc = currentOcc,noPopED = TRUE)
  return(dfg_db0) 
}

