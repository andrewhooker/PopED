## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m2 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db){
# M2 derivative of the vectorized variance w$r.t. bpops
# 
# the output is a matrix with dimensions (ind_samps^2 X nbpop)
# create a (n^2 x nbpop) matrix
  dv_dbeta = grad_bpop(helper_v_EBE,5,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db)
  return(list(dv_dbeta = dv_dbeta, poped.db = poped.db))
}
