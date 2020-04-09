#' The Fisher Information Matrix (FIM) for one individual
#' 
#' Compute the FIM for one individual given specific model(s), parameters, design and methods. 
#' 
#' @param xt A vector of sample times.  
#' @param model_switch A vector that is the same size as xt, specifying which model each sample belongs to.
#' @param x A vector for the discrete design variables.
#' @param a A vector of covariates.  
#' @inheritParams mftot
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped.db}{A PopED database}
#' 
#' @family FIM
#' 
# @example tests/testthat/examples_fcn_doc/warfarin_basic.R
# @example tests/testthat/examples_fcn_doc/examples_mf3.R
# @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf3 <- function(model_switch,xt,x,a,bpop,d,sigma,docc,poped.db){

  numnotfixed_bpop = sum(poped.db$parameters$notfixed_bpop)
  numnotfixed_d    = sum(poped.db$parameters$notfixed_d)
  numnotfixed_covd = sum(poped.db$parameters$notfixed_covd)
  numnotfixed_docc  = sum(poped.db$parameters$notfixed_docc)
  numnotfixed_covdocc  = sum(poped.db$parameters$notfixed_covdocc)
  numnotfixed_sigma  = sum(poped.db$parameters$notfixed_sigma)
  numnotfixed_covsigma  = sum(poped.db$parameters$notfixed_covsigma)
  
  n_fixed_eff <- numnotfixed_bpop
  n_rand_eff <- numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma

  n=size(xt,1)
  ret = 0
  
  for(i in 1:poped.db$settings$iFOCENumInd){
    b_ind = poped.db$parameters$b_global[,i,drop=F]
    bocc_ind = poped.db$parameters$bocc_global[[i]]
    
    if((poped.db$settings$bCalculateEBE) ){#Calculate an EBE
      epsi0 = zeros(1,length(poped.db$parameters$notfixed_sigma))
      g=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind)
      returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt,g,epsi0,poped.db) 
      mean_data <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      start_bind = t(b_ind)
      b_ind = ind_estimates(mean_data,bpop,d,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt,x,a,b_ind,bocc_ind,poped.db)
      #        b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(b_ind),(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
      
      #b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(zeros(size(b_ind)[1],size(b_ind)[2])),!(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
      poped.db$mean_data = mean_data
    }
    
    if(n_fixed_eff!=0){
      returnArgs <- m1(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,poped.db)
      m1_tmp <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
    } else {
      m1_tmp <- 0
    }
    if(n_rand_eff!=0){
      bUseVarSigmaDerivative = poped.db$settings$iFIMCalculationType != 4
      returnArgs <- m3(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,bUseVarSigmaDerivative,poped.db)
      m3_tmp <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
    } else {
      m3_tmp <- 0
    }
    if(n_fixed_eff!=0 & n_rand_eff!=0 & poped.db$settings$iFIMCalculationType %in% c(0,5,6)){
      returnArgs <- m2(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
      m2_tmp <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
    } else {
      m2_tmp <- 0
    }

    returnArgs <-  v(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
    v_tmp <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]

    if (!poped.db$settings$iFIMCalculationType %in% c(5,7)) {
      f1=zeros(n+n*n,n_fixed_eff+n_rand_eff)
      f1[1:n,1:n_fixed_eff] <- m1_tmp
      f1[(n+1):(n+n*n),1:n_fixed_eff] <- m2_tmp
      if(n_rand_eff!=0) f1[(n+1):(n+n*n),(n_fixed_eff+1):(n_fixed_eff+n_rand_eff)] <- m3_tmp
      
      if(any(v_tmp!=0)){#If there are some non-zero elements to v_tmp
        f2=zeros(n+n*n,n+n*n)
        
        v_tmp_inv = inv(v_tmp,pseudo_on_fail = T)
        f2[1:n,1:n] = v_tmp_inv
        
        tmp_m4_inv <- 1/2*kronecker(v_tmp_inv,v_tmp_inv)
        f2[(n + 1):(n + n*n), (n + 1):(n + n*n)] = tmp_m4_inv
        ret = ret + t(f1) %*% f2 %*% f1
      } else {
        ret = ret + t(f1) %*% f1
      }
    } else {
      v_tmp_inv = inv(v_tmp)
      if(n_rand_eff!=0) dim(m3_tmp) = c(n,n,n_rand_eff)
      
      tmp_fim = zeros(n_fixed_eff + n_rand_eff, n_fixed_eff + n_rand_eff)
      tmp_fim[1:n_fixed_eff,1:n_fixed_eff] = 2*t(m1_tmp) %*% v_tmp_inv %*% m1_tmp

      if (is.matrix(m2_tmp)) {
        dim(m2_tmp) = c(n,n,n_fixed_eff)
        for (m in 1:n_fixed_eff) {
          for (k in 1:n_fixed_eff) {
            tmp_fim[m,k] = tmp_fim[m,k] + trace_matrix(m2_tmp[,,m]%*%v_tmp_inv%*%m2_tmp[,,k]%*%v_tmp_inv)
          }
        }
        if(n_rand_eff!=0){
          for(m in 1:n_rand_eff){
            for(k in 1:n_fixed_eff){
              num = trace_matrix(m3_tmp[,,m]%*%v_tmp_inv%*%m2_tmp[,,k]%*%v_tmp_inv)
              tmp_fim[n_fixed_eff + m, k]=num
              tmp_fim[k, n_fixed_eff + m]=num
            }
          }
        }
      }
      if(n_rand_eff!=0){
        for (m in 1:n_rand_eff) {
          for (k in 1:n_rand_eff) {
            tmp_fim[n_fixed_eff + m, n_fixed_eff + k] = trace_matrix(m3_tmp[,,m] %*% v_tmp_inv %*% m3_tmp[,,k] %*% v_tmp_inv)
          }
        }
      }
      ret = ret + 1/2*tmp_fim
    }
  }
  ret = ret/poped.db$settings$iFOCENumInd
  
  return(list(ret = ret,poped.db = poped.db)) 
}
