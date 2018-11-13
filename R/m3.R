## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m3 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,bUseVarSigmaDerivative,poped.db){
  #
  # size: (samps per subject^2 x (number of random effects + number of occasion variances + number of sigmas))

  dv_dd = NULL
  dv_covd = NULL
  dv_ddocc = NULL
  dv_covdocc = NULL
  dv_dsig = NULL
  dv_dcovsig = NULL

  if (sum(poped.db$parameters$notfixed_d)>0) {
    dv_dd = grad_bpop(helper_v_EBE,8,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db, subset=poped.db$parameters$notfixed_d)
  }

  if (sum(poped.db$parameters$notfixed_covd)>0) {
    dv_covd = grad_bpop(helper_v_EBE,8,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db, subset=poped.db$parameters$notfixed_covd, offdiag = TRUE)
  }

  if (sum(poped.db$parameters$notfixed_docc)>0) {
    dv_ddocc = grad_bpop(helper_v_EBE,10,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db, subset=poped.db$parameters$notfixed_docc)
  }

  if (sum(poped.db$parameters$notfixed_covdocc)>0) {
    dv_covdocc = grad_bpop(helper_v_EBE,10,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db, subset=poped.db$parameters$notfixed_covdocc, offdiag = TRUE)
  }

  if (sum(poped.db$parameters$notfixed_sigma)>0) {
    dv_dsig = grad_bpop(helper_v_EBE,9,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db, subset=poped.db$parameters$notfixed_sigma)
    if (bUseVarSigmaDerivative == FALSE) {
      dv_dsig = t(2*sqrt(diag(sigma))[poped.db$parameters$notfixed_sigma==1] * t(dv_dsig))
    }
  }

  if (sum(poped.db$parameters$notfixed_covsigma)>0) {
    dv_dcovsig = grad_bpop(helper_v_EBE,9,size(xt_ind,1)^2,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db, subset=poped.db$parameters$notfixed_covsigma, offdiag = TRUE)
  }

  dv_db = do.call(cbind, list(dv_dd, dv_covd, dv_ddocc, dv_covdocc, dv_dsig, dv_dcovsig))

  return(list(dv_db = dv_db, poped.db = poped.db)) 
}
