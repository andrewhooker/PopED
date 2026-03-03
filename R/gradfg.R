## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradfg <- function(x,a,bpop,b_ind,bocc_ind,poped.db){
#
#
# size: (number of g's x number of random effects)
#
# deriv of g's w.r.t. b's and bocc's eval at b_ind
#
  dfg_db0 = grad_all(poped.db$model$fg_pointer,4,poped.db$parameters$ng,x,a,bpop,b_ind,bocc_ind,poped.db,noPopED = TRUE)
  return(dfg_db0) 
}

