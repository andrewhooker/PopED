## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m1 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,poped.db){
  #
  # function computes the derivative of the
  # linerarized model function w$r.t. bpop
  # for an individual
  #
  # the output is a matrix with dimensions (ind_samps X nbpop)
  df_dbeta = grad_bpop(helper_LinMatrix,5,size(xt_ind,1),model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc=NULL,poped.db)

  return(list( df_dbeta= df_dbeta,poped.db=poped.db)) 
}

