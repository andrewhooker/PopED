## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

grad_bpop <- function(func,select_par,nout,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db,subset=poped.db$parameters$notfixed_bpop, offdiag = FALSE){
  #----------Model linearization with respect to pop parameters
  #
  # use helper function to check for/include EBEs
  #
  dx_dbpop = grad_all(func,select_par,nout,model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db,subset=subset,noPopED = FALSE,offdiag=offdiag)
  
  return(dx_dbpop=dx_dbpop) 
}

# helper for m2
helper_v_EBE <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) {
  
  if((poped.db$settings$bCalculateEBE)){
    #zeros(size(b_ind)[1],size(b_ind)[2])
    b_ind_x = ind_estimates(poped.db$mean_data,bpop,d,sigma,t(b_ind),(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
  } else {
    b_ind_x = b_ind
  }
  
  vx <- v(model_switch,xt_ind,x,a,bpop,b_ind_x,bocc_ind,d,sigma,docc,poped.db)[[1]]
  dim(vx) = c(prod(dim(vx)),1)
  return(list(vx=vx, poped.db=poped.db))
}

# helper for m1
helper_LinMatrix <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) {
  
  epsi0 = zeros(1,length(poped.db$parameters$notfixed_sigma))
  
  # create linearized model
  if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==3) ){#FO, FOI
    b_ind=zeros(poped.db$parameters$NumRanEff,1)
  }
  
  if((poped.db$settings$bCalculateEBE)){
    b_ind = ind_estimates(poped.db$mean_data,bpop,d,sigma,t(b_ind),(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
  }
  
  g_p=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind)
  
  returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_p,epsi0,poped.db)
  ferror <- returnArgs[[1]]
  
  if ( (poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==3 || (isempty(b_ind) && isempty(bocc_ind))) ) {
    #FO, FOI
    if((poped.db$settings$bUseSecondOrder)) {
      hess_eta = zeros(length(xt_ind),1)
      for (o in 1:length(xt_ind)) {
        hessian_eta = hessian_eta_complex(model_switch[o],xt_ind[o],x,a,bpop,b_ind,bocc_ind,poped.db)
        hess_eta[o] = 1/2*trace_matrix(hessian_eta*d)
      }
      ferror = ferror + hess_eta
    }
    
  } else {
    #FOCE, FOCEI
    returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db)
    l_plus <- returnArgs[[1]]

    if((isempty(b_ind)) ){#No IIV present
      l_plus = zeros(size(xt_ind,1), 1)
    } else {
      l_plus = l_plus%*%b_ind
    }

    occ_add_plus = zeros(size(xt_ind,1), 1)
    if (poped.db$parameters$NumOcc!=0) {
      for (m in 1:poped.db$parameters$NumOcc) {
        returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,m,poped.db)
        l_plus_occ <- returnArgs[[1]]
        occ_add_plus = occ_add_plus + l_plus_occ*(bocc_ind[,m])
      }
    }
    ferror = ferror-(l_plus+occ_add_plus)
  }
  return(list(ferror = ferror, poped.db = poped.db))
}
