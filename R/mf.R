#' The full  Fisher Information Matrix (FIM) for one individual
#' 
#' Compute the full  FIM for one individual given specific model(s), parameters, design and methods. 
#' This computation makes no assumption that fixed and random effects are uncorrelated.  
#' 
#' @param xt_ind A vector of sample times.  
#' @param model_switch A vector that is the same size as xt, specifying which model each sample belongs to.
#' @param x A vector for the discrete design variables.
#' @param a A vector of covariates.  
#' 
#' @inheritParams mftot
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{globalStructure}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot0}}.  
#'  
#' @family FIM

## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf <- function(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure){

numnotfixed_bpop = sum(globalStructure$notfixed_bpop)
numnotfixed_d    = sum(globalStructure$notfixed_d)
numnotfixed_covd = sum(globalStructure$notfixed_covd)
numnotfixed_docc  = sum(globalStructure$notfixed_docc)
numnotfixed_covdocc  = sum(globalStructure$notfixed_covdocc)
numnotfixed_sigma  = sum(globalStructure$notfixed_sigma)
numnotfixed_covsigma  = sum(globalStructure$notfixed_covsigma)

n=size(xt_ind,1)
ret = 0

for(i in 1:globalStructure$iFOCENumInd){
    b_ind = globalStructure$b_global[,i,drop=F]
    bocc_ind = globalStructure$bocc_global[[i]]

    if((globalStructure$bCalculateEBE) ){#Calculate an EBE
        epsi0 = zeros(1,length(globalStructure$notfixed_sigma))
        g=feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_ind)
         returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt_ind,g,epsi0,globalStructure) 
mean_data <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
        start_bind = t(b_ind)
        b_ind = ind_estimates(mean_data,bpop,d,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
        #        b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(b_ind),(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)

        #b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(zeros(size(b_ind)[1],size(b_ind)[2])),!(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
        globalStructure$mean_data = mean_data
    }
    
    f1=zeros(n+n*n,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
     returnArgs <- m1(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,globalStructure) 
m1_tmp <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    f1[1:n,1:numnotfixed_bpop]=m1_tmp
     returnArgs <- m2(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,globalStructure) 
f1[(n+1):(n+n*n),1:numnotfixed_bpop] <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    if((numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma>0)){
         returnArgs <- m3(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,globalStructure) 
f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)] <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    }
    f2=zeros(n+n*n,n+n*n)
     returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,globalStructure) 
v_tmp <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    if((matrix_any(v_tmp)!=0)){
        #f2[1:n,1:n]=v_tmp\diag_matlab(1,size(v_tmp))
        #v_tmp\diag_matlab(1,size(v_tmp))
        f2[1:n,1:n]=inv(v_tmp)
        #m4_tmp = m4(v_tmp,n)
        #f2[(n+1):(n+n*n),(n+1):(n+n*n)]=m4_tmp\diag_matlab(1,size(m4_tmp))
        f2[(n+1):(n+n*n),(n+1):(n+n*n)]=inv(m4(v_tmp,n))
    }
    
    if((sum(sum(f2!=0))==0) ){#Only zeros in f2, uses FIM = m1'*m1
        ret = ret+t(m1_tmp)*m1_tmp
    } else {
        ret=ret+t(f1)%*%f2%*%f1
    }

}
#ret
ret = ret/globalStructure$iFOCENumInd



return(list( ret= ret,globalStructure=globalStructure)) 
}

