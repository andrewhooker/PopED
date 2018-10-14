## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradff <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db){
  #----------Model linearization with respect to random var.
  #
  # size of return is (samples per individual x number of g's)
  #
  # derivative of model w.r.t. g eval at b=b_ind
  #
  #
  fg0 = feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind)
  epsi0 = zeros(1,length(poped.db$parameters$notfixed_sigma))

  dff_dg0 = grad_all(poped.db$model$ferror_pointer,3,size(xt_ind,1),model_switch,xt_ind,fg0,epsi0,poped.db)

  return(list( dff_dg0= dff_dg0,poped.db=poped.db)) 
}
