## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

LinMatrixLH <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,NumEPS,globalStructure){
#----------Model linearization with respect to epsilon.
#
# size of return is (samples per individual x (number of sigma x number of omega)) 
#
# derivative of model w$r.t. sigma then eta, eval at e=0 and eta
#


y = zeros(size(xt_ind,1),globalStructure$NumRanEff*NumEPS)
if((globalStructure$iApproximationMethod==0 || globalStructure$iApproximationMethod==1) ){#No interaction
    #return
    return(list( y= y,globalStructure=globalStructure)) 
}
if((globalStructure$NumRanEff==0)){
   #return
   return(list( y= y,globalStructure=globalStructure)) 
}

if((globalStructure$hle_switch==30) ){#Automatic differentiation (INTLab)
  stop("Automatic differentiation not yet implemented in R version of PopED")
#     e0=zeros(1,NumEPS)
#     init_vec = matrix(c(t(b_ind), e0),nrow=1,byrow=T)
#     deriv_vec = hessianinit(init_vec)
#      returnArgs <-  new_ferror_file(model_switch,deriv_vec,xt_ind,x,a,bpop,bocc_ind,globalStructure) 
# deriv <- returnArgs[[1]]
# globalStructure <- returnArgs[[2]]
#     cellDeriv = cell(1,globalStructure$NumRanEff)
#      for(i in 1:globalStructure$NumRanEff){
#         tmp=zeros(size(xt_ind,1),NumEPS)
#         for(j in 1:size(xt_ind,1) ){#for each sample
#             val = deriv(j)$hx
#             tmp(j,)=val(globalStructure$NumRanEff+1:globalStructure$NumRanEff+NumEPS,i)
#         }
#         cellDeriv[[i]] = tmp
#      }
#     
#      for(i in 1:globalStructure$NumRanEff){
#         tmp = cellDeriv[[i]]
#         y[,(i-1)*NumEPS+1:i*NumEPS]=tmp(,1:NumEPS)
#      }
#     #return
#     return(list( y= y,globalStructure=globalStructure)) 
}

if((globalStructure$hle_switch==20)){
    stop(sprintf('Analytic derivative with interaction is not yet available!'))
}

for(i in 1:globalStructure$NumRanEff){
    b_ind_plus=b_ind
    b_ind_minus=b_ind
    b_ind_plus[i] = b_ind_plus[i]+globalStructure$hle
    b_ind_minus[i]= b_ind_minus[i]-globalStructure$hle
     returnArgs <-  LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind_plus,bocc_ind,globalStructure) 
lin_plus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
     returnArgs <-  LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind_minus,bocc_ind,globalStructure) 
lin_minus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    temp = (lin_plus-lin_minus)/(2*globalStructure$hle)
    y[,(i-1)*NumEPS+1:i*NumEPS]=temp[,1:NumEPS,drop=F]
}
return(list( y= y,globalStructure=globalStructure)) 
}

#Helper function to get the hessian for the AD derivative
new_ferror_file <- function(model_switch,deriv_vec,xt_ind,x,a,bpop,bocc_ind,globalStructure){
    fg0=feval(globalStructure$fg_pointer,x,a,bpop,deriv_vec(1:globalStructure$NumRanEff),bocc_ind) #Interaction
     returnArgs <- feval(globalStructure$ferror_pointer,model_switch,xt_ind,fg0,deriv_vec(globalStructure$NumRanEff+1:length(deriv_vec)),globalStructure) 
f_error <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
return(list( f_error= f_error,globalStructure=globalStructure)) 
}
