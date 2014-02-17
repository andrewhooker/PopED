## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixH <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure){
#----------Model linearization with respect to epsilon.
#
# size of return is (samples per individual x number of epsilons) 
#
# derivative of model w$r.t. eps eval at e=0
#
  NumEPS = size(globalStructure$sigma,1)
  if((NumEPS==0)){
	y=0
  } else { 
     returnArgs <- gradf_eps(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,NumEPS,globalStructure) 
y <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
  }
return(list( y= y,globalStructure=globalStructure)) 
}


gradf_eps <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,num_eps,globalStructure){
#----------Model linearization with respect to epsilon.
#
# size of return is (samples per individual x number of epsilons)
#
# derivative of model w$r.t. eps eval at e=0 and b=b_ind
#
#

dfeps_de0=zeros(size(xt_ind,1),num_eps)

if((globalStructure$iApproximationMethod==0 || globalStructure$iApproximationMethod==1) ){#No interaction
    fg0=feval(globalStructure$fg_pointer,x,a,bpop,zeros(size(b_ind)[1],size(b_ind)[2]),zeros(size(bocc_ind)[1],size(bocc_ind)[2]))

} else {
    fg0=feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_ind) #Interaction
}

e0=zeros(1,num_eps)

#Central approximation
if((globalStructure$hle_switch==1)){
    for(i in 1:num_eps){
        e_plus=e0
        e_minus=e0
        e_plus[i] = e_plus[i]+globalStructure$hle
        e_minus[i]= e_minus[i]-globalStructure$hle
         returnArgs <-  feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg0,e_plus,globalStructure) 
ferror_plus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
         returnArgs <-  feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg0,e_minus,globalStructure) 
ferror_minus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
        dfeps_de0[,i]=(ferror_plus-ferror_minus)/(2*globalStructure$hle)
    }
} else {
    #Complex approximation
    if((globalStructure$hle_switch==0)){
        for(i in 1:num_eps){
            e_plus=e0
e_plus[i] = complex(real=e_plus[i],imaginary=globalStructure$hle)
             returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg0,e_plus,globalStructure) 
ferror_plus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
            dfeps_de0[,i]=Im(ferror_plus)/globalStructure$hle
        }
    } else {
        if((globalStructure$hle_switch==30) ){#Automatic differentiation (INTLab)
          stop("Automatic differentiation not yet implemented in PopED for R")
#             e_init = gradientinit(e0)
#              returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg0,e_init,globalStructure) 
# ferror_val <- returnArgs[[1]]
# globalStructure <- returnArgs[[2]]
#             dfeps_de0 = ferror_val$dx
        } else {
            if((globalStructure$hle_switch!=20)){
                stop(sprintf('Unknown derivative option for gradf_eps'))
            }
        }
    }
}
return(list( dfeps_de0= dfeps_de0,globalStructure=globalStructure)) 
}
