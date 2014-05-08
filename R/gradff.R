## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradff <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure){
  #----------Model linearization with respect to random var.
  #
  # size of return is (samples per individual x number of g's)
  #
  # derivative of model w$r.t. g eval at b=b_ind
  #
  #
  dff_dg0=zeros(size(xt_ind,1),globalStructure$ng)
  fg0=feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_ind)
  
  epsi0 = zeros(1,length(globalStructure$notfixed_sigma))
  
  #Central approximation
  if((globalStructure$gradff_switch[1] == 1)){
    for(i in 1:globalStructure$ng){
      g_plus=fg0
      g_minus=fg0
      g_plus[i]=g_plus[i]+globalStructure$hlf
      g_minus[i]=g_minus[i]-globalStructure$hlf
      returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,globalStructure) 
      ferror_plus <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
      returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt_ind,g_minus,epsi0,globalStructure) 
      ferror_minus <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
      dff_dg0[,i]=(ferror_plus-ferror_minus)/(2.0*globalStructure$hlf)
    }
  } else {
    #Complex approximation
    if((globalStructure$gradff_switch[1]==0)){
      for(i in 1:globalStructure$ng){
        g_plus=fg0
        g_plus[i] = complex(real=g_plus[i],imaginary=globalStructure$hlf)
        returnArgs <-  feval(globalStructure$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,globalStructure) 
        ferror_pointer_plus <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
        dff_dg0[,i]=Im(ferror_pointer_plus)/globalStructure$hlf
      }
    } else {
      if((globalStructure$gradff_switch[1] == 20) ){#Analytic derivative
        for(k in 1:size(xt_ind,1)){
          dff_dg0[k,] = eval(sprintf('analytic_dff_dg%d(model_switch,xt_ind[k],fg0)',model_switch[k]))
        }
      } else {
        if((globalStructure$gradff_switch[1] == 30) ){#Calculate using automatic differentiation (INTLab)
          stop("Automatic differentiation not currently implemented in PopED for R")
          #                 fg_init = gradientinit(fg0)
          #                  returnArgs <-  feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg_init,epsi0,globalStructure) 
          # val <- returnArgs[[1]]
          # globalStructure <- returnArgs[[2]]
          #                 dff_dg0 = val$dx
        } else {
          stop(sprintf('Unknown derivative option for gradff'))
        }
      }
    }
  }
  return(list( dff_dg0= dff_dg0,globalStructure=globalStructure)) 
}
