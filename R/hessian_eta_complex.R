#Hessian over eta, evaluated at eta_hat => laplace approximation possible
hessian_eta_complex <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure,return_gradient=F){
  
  bAutomatic = FALSE
  epsi0 = zeros(1,length(globalStructure$notfixed_sigma))
  n=length(b_ind)
  hess=zeros(n)   # Memory for the Hessian matrix
  g=zeros(n,1)    # Memory for the gradient vector
  
  if((bAutomatic)){
    stop("Automatic differentiation not yet implemented in R version of PopED")
    #     if((globalStructure$Engine$Type==2) ){#FreeMat
    #         stop(sprintf('Automatic differentiation is not available in PopED with FreeMat'))
    #     }
    #     b_init = hessianinit(b_ind)
    #     fg_init=feval(globalStructure$fg_pointer,x,a,bpop,b_init,bocc_ind)
    #      returnArgs <-  feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg_init,epsi0,globalStructure) 
    # val <- returnArgs[[1]]
    # globalStructure <- returnArgs[[2]]
    #     hess = val$hx
  } else {
    h=1E-04
    h2=h*h
    for(k in 1:n){
      eta_plus = b_ind
      eta_plus[k] = eta_plus[k]+h*1i
      
      if((return_gradient)){
        g_plus = feval(globalStructure$fg_pointer,x,a,bpop,eta_plus,bocc_ind)
        ff_plus = feval(globalStructure$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,globalStructure)
        g[k]=Im(ff_plus)/h             # the kth gradient
      }
      for(l in k:n                       ){# Hessian (off-diagonal)
        eta_plus2 = eta_plus
        eta_plus2[l] = eta_plus2[l]+h
        g_plus = feval(globalStructure$fg_pointer,x,a,bpop,eta_plus2,bocc_ind)
        ff_plus = feval(globalStructure$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,globalStructure)
        
        eta_plus2[l]=eta_plus[l]-h
        g_plus = feval(globalStructure$fg_pointer,x,a,bpop,eta_plus2,bocc_ind)
        ff_minus = feval(globalStructure$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,globalStructure)
        
        hess[k,l]=sum(Im(ff_plus-ff_minus)/h2/2)    # Hessian (central + complex step)
        hess[l,k]=hess[k,l]                           #Make hessian symmetric
      }
    }
  }
  return(list( hess= hess,g=g)) 
}
