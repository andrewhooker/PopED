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
#' \item{globalStructure}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot1}}.  
#' @family FIM
#' 
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf3 <- function(model_switch,xt,x,a,bpop,d,sigma,docc,globalStructure){
#Calculate the reduced FIM

numnotfixed_bpop = sum(globalStructure$notfixed_bpop)
numnotfixed_d    = sum(globalStructure$notfixed_d)
numnotfixed_covd = sum(globalStructure$notfixed_covd)
numnotfixed_docc  = sum(globalStructure$notfixed_docc)
numnotfixed_covdocc  = sum(globalStructure$notfixed_covdocc)
numnotfixed_sigma  = sum(globalStructure$notfixed_sigma)
numnotfixed_covsigma  = sum(globalStructure$notfixed_covsigma)

n=size(xt,1)
ret = 0

for(i in 1:globalStructure$iFOCENumInd){
    b_ind = globalStructure$b_global[,i,drop=F]
    bocc_ind = globalStructure$bocc_global[[i]]
    
    if((globalStructure$bCalculateEBE) ){#Calculate an EBE
        epsi0 = zeros(1,length(globalStructure$notfixed_sigma))
        g=feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_ind)
         returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt,g,epsi0,globalStructure) 
mean_data <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
        start_bind = t(b_ind)
        b_ind = ind_estimates(mean_data,bpop,d,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt,x,a,b_ind,bocc_ind,globalStructure)
#        b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(b_ind),(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
        
        #b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(zeros(size(b_ind)[1],size(b_ind)[2])),!(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
        globalStructure$mean_data = mean_data
    }
    
    
    f1=zeros(n+n*n,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
     returnArgs <- m1(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,globalStructure) 
f1[1:n,1:numnotfixed_bpop] <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
     returnArgs <- m3(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,globalStructure) 
f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)] <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]

    f2=zeros(n+n*n,n+n*n)
     returnArgs <-  v(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,globalStructure) 
v_tmp <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
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
ret = ret/globalStructure$iFOCENumInd

return(list( ret= ret,globalStructure=globalStructure)) 
}
