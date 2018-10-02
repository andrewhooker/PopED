## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m3 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,bUseVarSigmaDerivative,poped.db){
  
  #
  # size: (samps per subject^2 x (number of random effects + number of occasion variances + number of sigmas))
  
  NumSigma = size(sigma,1)
  NumDocc = size(docc,1)
  
  if((isempty(sigma))){
    NumSigma=0
  }
  
  ns=size(xt_ind,1)^2
  
  dv_db_new=zeros(ns,sum(poped.db$parameters$notfixed_d)+
                    sum(poped.db$parameters$notfixed_covd)+
                    sum(poped.db$parameters$notfixed_docc)+
                    sum(poped.db$parameters$notfixed_covdocc)+
                    sum(poped.db$parameters$notfixed_sigma)+
                    sum(poped.db$parameters$notfixed_covsigma))
  
  #If complex or central differentiation
    k=1
    if(!isempty(d)){
      for(i in 1:poped.db$parameters$NumRanEff){
        if((poped.db$parameters$notfixed_d[i]==1)){
          d_plus=d
          
          # Central approximation
          d_plus[i,i]=d_plus[i,i]+poped.db$settings$hm2
          d_minus=d
          d_minus[i,i]=d_minus[i,i]-poped.db$settings$hm2
          
          if((poped.db$settings$bCalculateEBE)){
            start_bind = t(b_ind)%*%zeros(size(t(b_ind)))%*%t(b_ind)
            b_ind_plus = ind_estimates(poped.db$mean_data,bpop,d_plus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
            b_ind_minus = ind_estimates(poped.db$mean_data,bpop,d_minus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
            # b_ind_plus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d+poped.db$settings$hm2))
            # b_ind_minus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d-poped.db$settings$hm2))
            #b_ind_plus = b_ind
            #b_ind_minus = b_ind
          } else {
            b_ind_plus = b_ind
            b_ind_minus = b_ind
          }
          
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d_plus,sigma,docc,poped.db) 
          v_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d_minus,sigma,docc,poped.db) 
          v_minus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          dv=v_plus-v_minus
          if((!isempty(dv))){
            ir=dv/(2*poped.db$settings$hm2)
            ir=as.vector(ir)
            dv_db_new[,k]=ir
          }
          k=k+1
        }
      }  
      
      for(i in 1:length(poped.db$parameters$notfixed_covd)){
        if((poped.db$parameters$notfixed_covd[i]==1)){
          
          d_plus=update_offdiag(d,i,poped.db$settings$hm2)
          d_minus=update_offdiag(d,i,-poped.db$settings$hm2)
          
          if((poped.db$settings$bCalculateEBE)){
            start_bind = t(b_ind)%*%zeros(size(t(b_ind)))%*% t(b_ind)
            b_ind_plus = ind_estimates(poped.db$mean_data,bpop,d_plus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
            b_ind_minus = ind_estimates(poped.db$mean_data,bpop,d_minus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
            # b_ind_plus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d+poped.db$settings$hm2))
            # b_ind_minus = b_ind/sqrt(diag_matlab(d))*sqrt(diag_matlab(d-poped.db$settings$hm2))
            #b_ind_plus = b_ind
            #b_ind_minus = b_ind
          } else {
            b_ind_plus = b_ind
            b_ind_minus = b_ind
          }
          
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d_plus,sigma,docc,poped.db) 
          v_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d_minus,sigma,docc,poped.db) 
          v_minus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          dv=v_plus-v_minus
          if((!isempty(dv))){
            ir=dv/(2*poped.db$settings$hm2)
            ir=as.vector(ir)
            dv_db_new[,k]=ir
          }
          k=k+1
        }
      }
    }
    
    if((!isempty(docc))){
      for(i in 1:NumDocc){
        if((poped.db$parameters$notfixed_docc[i]==1)){
          docc_plus=docc
          
          # Central approximation
          docc_plus[i,i]=docc_plus[i,i]+poped.db$settings$hm2
          docc_minus=docc
          docc_minus[i,i]=docc_minus[i,i]-poped.db$settings$hm2
          
          if((poped.db$settings$bCalculateEBE)){
            start_bind = t(b_ind)
            warning('EBE calculation with occasions is not available in the current version!')
            b_ind_plus = b_ind#ind_estimates(poped.db$mean_data,bpop,d_plus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
            b_ind_minus = b_ind#ind_estimates(poped.db$mean_data,bpop,d_minus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
          } else {
            b_ind_plus = b_ind
            b_ind_minus = b_ind
          }
          
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma,docc_plus,poped.db) 
          v_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma,docc_minus,poped.db) 
          v_minus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          dv=v_plus-v_minus
          if((!isempty(dv))){
            ir=dv/(2*poped.db$settings$hm2)
            ir=as.vector(ir)
            dv_db_new[,k]=ir
          }
          k=k+1
        }
      }
      
      if(!isempty(poped.db$parameters$notfixed_covdocc)){
        for(i in 1:length(poped.db$parameters$notfixed_covdocc)){
          if((poped.db$parameters$notfixed_covdocc[i]==1)){
            
            docc_plus=update_offdiag(docc,i,poped.db$settings$hm2)
            docc_minus=update_offdiag(docc,i,-poped.db$settings$hm2)
            
            if((poped.db$settings$bCalculateEBE)){
              start_bind = t(b_ind)
              warning('EBE calculation with covariance of occasions is not available in the current version!')
              b_ind_plus = b_ind#ind_estimates(poped.db$mean_data,bpop,d_plus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
              b_ind_minus = b_ind#ind_estimates(poped.db$mean_data,bpop,d_minus,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
            } else {
              b_ind_plus = b_ind
              b_ind_minus = b_ind
            }
            
            returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma,docc_plus,poped.db) 
            v_plus <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            returnArgs <- v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma,docc_minus,poped.db) 
            v_minus <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            dv=v_plus-v_minus
            if((!isempty(dv))){
              ir=dv/(2*poped.db$settings$hm2)
              ir=as.vector(ir)
              dv_db_new[,k]=ir
            }
            k=k+1
          }
        }
      }
      
    }
    for(i in 1:NumSigma){
      #     [S,R]=cov2corr(sigma) If off-diag covariances should be updated when
      #     differentiating a variance term
      if((poped.db$parameters$notfixed_sigma[i]==1)){
        sigma_plus=sigma
        
        # Central approximation
        sigma_plus[i,i]=sigma_plus[i,i]+poped.db$settings$hm2
        #        sigma_plus = corr2cov(sqrt(diag_matlab(sigma_plus)),R)
        sigma_minus=sigma
        sigma_minus[i,i]=sigma_minus[i,i]-poped.db$settings$hm2
        #       sigma_minus = corr2cov(sqrt(diag_matlab(sigma_minus)),R)
        
        if((poped.db$settings$bCalculateEBE)){
          start_bind = t(b_ind)
          b_ind_plus = ind_estimates(poped.db$mean_data,bpop,d,sigma_plus,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
          b_ind_minus = ind_estimates(poped.db$mean_data,bpop,d,sigma_minus,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma_plus,docc,poped.db) 
        v_plus <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma_minus,docc,poped.db) 
        v_minus <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        dv=v_plus-v_minus
        
        if((!isempty(dv))){
          if((bUseVarSigmaDerivative) ){#Derivative w$r.t. sigma as variance
            ir=dv/(2*poped.db$settings$hm2)
          } else {
            ir=2*sqrt(sigma[i,i])*dv/(2*poped.db$settings$hm2) #Derivative w$r.t. sigma as stdev
          }
          ir=as.vector(ir)
          dv_db_new[,k]=ir
        }
        k=k+1
      }
    }
    
    for(i in 1:length(poped.db$parameters$notfixed_covsigma)){
      if(any(size(poped.db$parameters$notfixed_covsigma)==0)) next
      if((poped.db$parameters$notfixed_covsigma[i]==1)){
        sigma_plus=update_offdiag(sigma,i,poped.db$settings$hm2)
        sigma_minus=update_offdiag(sigma,i,-poped.db$settings$hm2)
        
        if((poped.db$settings$bCalculateEBE)){
          start_bind = t(b_ind)
          b_ind_plus = ind_estimates(poped.db$mean_data,bpop,d,sigma_plus,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
          b_ind_minus = ind_estimates(poped.db$mean_data,bpop,d,sigma_minus,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,d,sigma_plus,docc,poped.db) 
        v_plus <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,d,sigma_minus,docc,poped.db) 
        v_minus <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        dv=v_plus-v_minus
        
        if((!isempty(dv))){
          #if (bUseVarSigmaDerivative) #Derivative w$r.t. sigma as variance
          ir=dv/(2*poped.db$settings$hm2)
          #else
          #    ir=2*sqrt(sigma[i,i])*dv/(2*poped.db$settings$hm2) #Derivative w$r.t. sigma as stdev
          #end
          ir=as.vector(ir)
          dv_db_new[,k]=ir
        }
        k=k+1
      }
    }
    
    ret = dv_db_new
  return(list( ret= ret,poped.db=poped.db)) 
}
