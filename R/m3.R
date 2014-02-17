## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m3 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,bUseVarSigmaDerivative,globalStructure){
  
  #
  # size: (samps per subject^2 x (number of random effects + number of occasion variances + number of sigmas))
  
  NumSigma = size(sigma,1)
  NumDocc = size(docc,1)
  
  if((isempty(sigma))){
    NumSigma=0
  }
  
  ns=size(xt_ind,1)^2
  
  dv_db_new=zeros(ns,sum(globalStructure$notfixed_d)+sum(globalStructure$notfixed_covd)+sum(globalStructure$notfixed_docc)+sum(globalStructure$notfixed_covdocc)+sum(globalStructure$notfixed_sigma)+sum(globalStructure$notfixed_covsigma))
  
  if((globalStructure$m2_switch[1] == 30) ){#Automatic differentiation of M3 dosen't work with ud variance term
    returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure) 
    l <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    returnArgs <- LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure) 
    h <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    returnArgs <- LinMatrixLH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,size(sigma,1),globalStructure) 
    lh <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    locc=cell(1,globalStructure$NumOcc)
    for(i in 1:globalStructure$NumOcc){
      if(globalStructure$NumOcc==0) next
      returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,i,globalStructure) 
      locc[[i]] <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
    }
    j=1
    #Differentiate the variance w$r.t iiv
    for(i in 1:globalStructure$NumRanEff){
      if((globalStructure$notfixed_d[i]==1)){
        dv_db_new[,j]=reshape_matlab(l[,i,drop=F]*t(l[,i,drop=F])+diag_matlab(diag_matlab(lh[,(i-1)*NumSigma+1:i*NumSigma,drop=F]*sigma*t(lh[,(i-1)*NumSigma+1:i*NumSigma,drop=F]))),ns,1) #Last term is interaction
        j=j+1
      }
    }
    
    if((sum(globalStructure$notfixed_covd)!=0)){
      for(i in 1:length(globalStructure$notfixed_covd) ){
        if((globalStructure$notfixed_covd[i]==1)){
          returnArgs <- get_cov_matrix_index(d,i)   
          m <- returnArgs[[1]]
          n <- returnArgs[[2]]
          if((m==-1 || n==-1)){
            stop(sprintf('Wrong index in get_covariance_matrix_index d, PopED is }ing!'))
          }
          lh1 = diag_matlab(diag_matlab(lh[,(m-1)*NumSigma+1:m*NumSigma,drop=F]*sigma*t(lh[,(n-1)*NumSigma+1:n*NumSigma,drop=F]))) #Interaction term
          lh2 = diag_matlab(diag_matlab(lh[,(n-1)*NumSigma+1:n*NumSigma,drop=F]*sigma*t(lh[,(m-1)*NumSigma+1:m*NumSigma,drop=F]))) #Interaction term
          dv_db_new[,j]=reshape_matlab(l[,m,drop=F]*t(l[,n,drop=F])+l[,n,drop=F]*t(l[,m,drop=F])+lh1+lh2,ns,1)
          j=j+1
        }
      }
    }
    
    #Differentiate the variance w$r.t occasion_variability
    for(i in 1:NumDocc){
      if((globalStructure$notfixed_docc[i]==1)){
        tmp = zeros(size(xt_ind,1),size(xt_ind,1))
        for(k in 1:globalStructure$NumOcc){
          tmp = tmp+locc[[k]][,i,drop=F]*t(locc[[k]][,i,drop=F])
        }
        dv_db_new[,j]=reshape_matlab(tmp,ns,1)
        j=j+1
      }
    }
    
    #Differentiate the variance w$r.t. occasion covariances
    if((sum(globalStructure$notfixed_covdocc)!=0)      ){
      for(i in 1:length(globalStructure$notfixed_covdocc) ){
        if((globalStructure$notfixed_covdocc[i]==1)){
          returnArgs <- get_cov_matrix_index(docc,i)   
          m <- returnArgs[[1]]
          n <- returnArgs[[2]]
          if((m==-1 || n==-1)){
            stop(sprintf('Wrong index in get_covariance_matrix_index docc, PopED is }ing!'))
          }
          tmp = zeros(size(xt_ind,1),size(xt_ind,1))
          for(k in 1:globalStructure$NumOcc){
            tmp = tmp+locc[[k]][,m,drop=F]*t(locc[[k]][,n,drop=F])+locc[[k]][,n,drop=F]*t(locc[[k]][,m,drop=F])
          }
          dv_db_new[,j]=reshape_matlab(tmp,ns,1)
          j=j+1
        }
      }
    }
    
    #Differentiate the variance w$r.t. sigma
    for(i in 1:NumSigma){
      if((globalStructure$notfixed_sigma[i]==1)){
        tmp_lh = zeros(size(xt_ind,1),globalStructure$NumRanEff)
        for(k in 1:globalStructure$NumRanEff ){#Only use the Random Eff interacting with sigma_i
          tmp_lh[,k] = lh[,i+(k-1)*NumSigma,drop=F]
        }
        if((bUseVarSigmaDerivative) ){#Derivative w$r.t. sigma as variance
          dv_db_new[,j]= reshape_matlab(diag_matlab(diag_matlab(h(,i)*t(h(,i))))+diag_matlab(diag_matlab(tmp_lh*d*t(tmp_lh))),ns,1)
        } else { #Derivarite w$r.t. sigma as stdev
          dv_db_new[,j]= reshape_matlab(2*sqrt(sigma[i,i])*diag_matlab(diag_matlab(h(,i)*t(h(,i))))+2*sqrt(sigma[i,i])*diag_matlab(diag_matlab(tmp_lh*d*t(tmp_lh))),ns,1)
        }
        j=j+1
      }
    }
    
    if((sum(globalStructure$notfixed_covsigma)!=0)){
      for(i in 1:length(globalStructure$notfixed_covsigma) ){
        if((globalStructure$notfixed_covsigma[i]==1)){
          returnArgs <- get_cov_matrix_index(sigma,i)   
          m <- returnArgs[[1]]
          n <- returnArgs[[2]]
          if((m==-1 || n==-1)){
            stop(sprintf('Wrong index in get_covariance_matrix_index sigma, PopED is }ing!'))
          }
          tmp_lh_m = zeros(size(xt_ind,1),globalStructure$NumRanEff)
          tmp_lh_n = zeros(size(xt_ind,1),globalStructure$NumRanEff)
          
          for(k in 1:globalStructure$NumRanEff ){#Only use the Random Eff interacting with sigma_m or sigma_n
            tmp_lh_m[,k] = lh[,m+(k-1)*NumSigma,drop=F]
            tmp_lh_n[,k] = lh[,n+(k-1)*NumSigma,drop=F]
          }
          dv_db_new[,j]=reshape_matlab(diag_matlab(diag_matlab(h(,m)*t(h(,n))))+diag_matlab(diag_matlab(h(,n)*t(h(,m))))+diag_matlab(diag_matlab(tmp_lh_m*d*t(tmp_lh_n)))+diag_matlab(diag_matlab(tmp_lh_n*d*t(tmp_lh_m))),ns,1)
          j=j+1
        }
      }
    }
    
    ret=dv_db_new
    
  } else { #If complex or central differentiation
    
    k=1
    
    for(i in 1:globalStructure$NumRanEff){
      if((globalStructure$notfixed_d[i]==1)){
        d_plus=d
        
        # Central approximation
        d_plus[i,i]=d_plus[i,i]+globalStructure$hm2
        d_minus=d
        d_minus[i,i]=d_minus[i,i]-globalStructure$hm2
        
        if((globalStructure$bCalculateEBE)){
          start_bind = t(b_ind)%*%zeros(size(t(b_ind)))%*%t(b_ind)
          b_ind_plus = ind_estimates(globalStructure$mean_data,bpop,d_plus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          b_ind_minus = ind_estimates(globalStructure$mean_data,bpop,d_minus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          # b_ind_plus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d+globalStructure$hm2))
          # b_ind_minus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d-globalStructure$hm2))
          #b_ind_plus = b_ind
          #b_ind_minus = b_ind
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d_plus,sigma,docc,globalStructure) 
        v_plus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d_minus,sigma,docc,globalStructure) 
        v_minus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        dv=v_plus-v_minus
        if((!isempty(dv))){
          ir=dv/(2*globalStructure$hm2)
          ir=reshape_matlab(ir,ns,1)
          dv_db_new[,k]=ir
        }
        k=k+1
      }
    }
    
    for(i in 1:length(globalStructure$notfixed_covd)){
      if((globalStructure$notfixed_covd[i]==1)){
        
        d_plus=update_offdiag(d,i,globalStructure$hm2)
        d_minus=update_offdiag(d,i,-globalStructure$hm2)
        
        if((globalStructure$bCalculateEBE)){
          start_bind = t(b_ind)%*%zeros(size(t(b_ind)))%*% t(b_ind)
          b_ind_plus = ind_estimates(globalStructure$mean_data,bpop,d_plus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          b_ind_minus = ind_estimates(globalStructure$mean_data,bpop,d_minus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          # b_ind_plus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d+globalStructure$hm2))
          # b_ind_minus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d-globalStructure$hm2))
          #b_ind_plus = b_ind
          #b_ind_minus = b_ind
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d_plus,sigma,docc,globalStructure) 
        v_plus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d_minus,sigma,docc,globalStructure) 
        v_minus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        dv=v_plus-v_minus
        if((!isempty(dv))){
          ir=dv/(2*globalStructure$hm2)
          ir=reshape_matlab(ir,ns,1)
          dv_db_new[,k]=ir
        }
        k=k+1
      }
    }
    
    if((!isempty(docc))){
      for(i in 1:NumDocc){
        if((globalStructure$notfixed_docc[i]==1)){
          docc_plus=docc
          
          # Central approximation
          docc_plus[i,i]=docc_plus[i,i]+globalStructure$hm2
          docc_minus=docc
          docc_minus[i,i]=docc_minus[i,i]-globalStructure$hm2
          
          if((globalStructure$bCalculateEBE)){
            start_bind = t(b_ind)
            warning('EBE calculation with occasions is not available in the current version!')
            b_ind_plus = b_ind#ind_estimates(globalStructure$mean_data,bpop,d_plus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
            b_ind_minus = b_ind#ind_estimates(globalStructure$mean_data,bpop,d_minus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          } else {
            b_ind_plus = b_ind
            b_ind_minus = b_ind
          }
          
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma,docc_plus,globalStructure) 
          v_plus <- returnArgs[[1]]
          globalStructure <- returnArgs[[2]]
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma,docc_minus,globalStructure) 
          v_minus <- returnArgs[[1]]
          globalStructure <- returnArgs[[2]]
          dv=v_plus-v_minus
          if((!isempty(dv))){
            ir=dv/(2*globalStructure$hm2)
            ir=reshape_matlab(ir,ns,1)
            dv_db_new[,k]=ir
          }
          k=k+1
        }
      }
      
      for(i in 1:length(globalStructure$notfixed_covdocc)){
        if((globalStructure$notfixed_covdocc[i]==1)){
          
          docc_plus=update_offdiag(docc,i,globalStructure$hm2)
          docc_minus=update_offdiag(docc,i,-globalStructure$hm2)
          
          if((globalStructure$bCalculateEBE)){
            start_bind = t(b_ind)
            warning('EBE calculation with covariance of occasions is not available in the current version!')
            b_ind_plus = b_ind#ind_estimates(globalStructure$mean_data,bpop,d_plus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
            b_ind_minus = b_ind#ind_estimates(globalStructure$mean_data,bpop,d_minus,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          } else {
            b_ind_plus = b_ind
            b_ind_minus = b_ind
          }
          
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma,docc_plus,globalStructure) 
          v_plus <- returnArgs[[1]]
          globalStructure <- returnArgs[[2]]
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma,docc_minus,globalStructure) 
          v_minus <- returnArgs[[1]]
          globalStructure <- returnArgs[[2]]
          dv=v_plus-v_minus
          if((!isempty(dv))){
            ir=dv/(2*globalStructure$hm2)
            ir=reshape_matlab(ir,ns,1)
            dv_db_new[,k]=ir
          }
          k=k+1
        }
      }
      
      
    }
    for(i in 1:NumSigma){
      #     [S,R]=cov2corr(sigma) If off-diag covariances should be updated when
      #     differentiating a variance term
      if((globalStructure$notfixed_sigma[i]==1)){
        sigma_plus=sigma
        
        # Central approximation
        sigma_plus[i,i]=sigma_plus[i,i]+globalStructure$hm2
        #        sigma_plus = corr2cov(sqrt(diag_matlab(sigma_plus)),R)
        sigma_minus=sigma
        sigma_minus[i,i]=sigma_minus[i,i]-globalStructure$hm2
        #       sigma_minus = corr2cov(sqrt(diag_matlab(sigma_minus)),R)
        
        if((globalStructure$bCalculateEBE)){
          start_bind = t(b_ind)
          b_ind_plus = ind_estimates(globalStructure$mean_data,bpop,d,sigma_plus,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          b_ind_minus = ind_estimates(globalStructure$mean_data,bpop,d,sigma_minus,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma_plus,docc,globalStructure) 
        v_plus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma_minus,docc,globalStructure) 
        v_minus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        dv=v_plus-v_minus
        
        if((!isempty(dv))){
          if((bUseVarSigmaDerivative) ){#Derivative w$r.t. sigma as variance
            ir=dv/(2*globalStructure$hm2)
          } else {
            ir=2*sqrt(sigma[i,i])*dv/(2*globalStructure$hm2) #Derivative w$r.t. sigma as stdev
          }
          ir=reshape_matlab(ir,ns,1)
          dv_db_new[,k]=ir
        }
        k=k+1
      }
    }
    
    for(i in 1:length(globalStructure$notfixed_covsigma)){
      if(any(size(globalStructure$notfixed_covsigma)==0)) next
      if((globalStructure$notfixed_covsigma[i]==1)){
        sigma_plus=update_offdiag(sigma,i,globalStructure$hm2)
        sigma_minus=update_offdiag(sigma,i,-globalStructure$hm2)
        
        if((globalStructure$bCalculateEBE)){
          start_bind = t(b_ind)
          b_ind_plus = ind_estimates(globalStructure$mean_data,bpop,d,sigma_plus,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
          b_ind_minus = ind_estimates(globalStructure$mean_data,bpop,d,sigma_minus,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma_plus,docc,globalStructure) 
        v_plus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma_minus,docc,globalStructure) 
        v_minus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        dv=v_plus-v_minus
        
        if((!isempty(dv))){
          #if (bUseVarSigmaDerivative) #Derivative w$r.t. sigma as variance
          ir=dv/(2*globalStructure$hm2)
          #else
          #    ir=2*sqrt(sigma[i,i])*dv/(2*globalStructure$hm2) #Derivative w$r.t. sigma as stdev
          #end
          ir=reshape_matlab(ir,ns,1)
          dv_db_new[,k]=ir
        }
        k=k+1
      }
    }
    
    ret = dv_db_new
  }
  return(list( ret= ret,globalStructure=globalStructure)) 
}
