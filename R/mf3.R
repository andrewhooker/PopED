#' The reduced  Fisher Information Matrix (FIM) for one individual
#' 
#' Compute the reduced  FIM for one individual given specific model(s), parameters, design and methods. 
#' This computation assumes that there is no correlation in the FIM between the fixed and random effects, 
#' and set these elements in the FIM to zero.
#' 
#' @param xt A vector of sample times.  
#' @inheritParams mf
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot1}}.  
#' @family FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_mf3.R
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf3 <- function(model_switch,xt,x,a,bpop,d,sigma,docc,poped.db){
  #Calculate the reduced FIM
  
  numnotfixed_bpop = sum(poped.db$notfixed_bpop)
  numnotfixed_d    = sum(poped.db$notfixed_d)
  numnotfixed_covd = sum(poped.db$notfixed_covd)
  numnotfixed_docc  = sum(poped.db$notfixed_docc)
  numnotfixed_covdocc  = sum(poped.db$notfixed_covdocc)
  numnotfixed_sigma  = sum(poped.db$notfixed_sigma)
  numnotfixed_covsigma  = sum(poped.db$notfixed_covsigma)
  
  n=size(xt,1)
  ret = 0
  
  for(i in 1:poped.db$iFOCENumInd){
    b_ind = poped.db$b_global[,i,drop=F]
    bocc_ind = poped.db$bocc_global[[i]]
    
    if((poped.db$bCalculateEBE) ){#Calculate an EBE
      epsi0 = zeros(1,length(poped.db$notfixed_sigma))
      g=feval(poped.db$fg_pointer,x,a,bpop,b_ind,bocc_ind)
      returnArgs <- feval(poped.db$ferror_pointer,model_switch,xt,g,epsi0,poped.db) 
      mean_data <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      start_bind = t(b_ind)
      b_ind = ind_estimates(mean_data,bpop,d,sigma,start_bind,(poped.db$iApproximationMethod==2),FALSE,model_switch,xt,x,a,b_ind,bocc_ind,poped.db)
      #        b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(b_ind),(poped.db$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
      
      #b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(zeros(size(b_ind)[1],size(b_ind)[2])),!(poped.db$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
      poped.db$mean_data = mean_data
    }
    
    
    f1=zeros(n+n*n,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
    returnArgs <- m1(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,poped.db) 
    f1[1:n,1:numnotfixed_bpop] <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    returnArgs <- m3(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,poped.db) 
    f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)] <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    f2=zeros(n+n*n,n+n*n)
    returnArgs <-  v(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
    v_tmp <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    if((matrix_any(v_tmp)!=0) ){#If the inverse is not empty
      f2[1:n,1:n]=inv(v_tmp)
      tmp_m4=m4(v_tmp,n)
      f2[(n+1):(n+n*n),(n+1):(n+n*n)]=inv(tmp_m4)
    }
    if((all(f2==0))){
      ret = ret+t(f1)*f1
    } else {
      ret = ret+t(f1)%*%f2%*%f1
    }
  }
  ret = ret/poped.db$iFOCENumInd
  
  return(list( ret= ret,poped.db=poped.db)) 
}
