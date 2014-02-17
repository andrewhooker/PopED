## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m2 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,globalStructure){
# M2 derivative of the vectorized variance w$r.t. bpops
# 
# the output is a matrix with dimensions (ind_samps^2 X nbpop)
# create a (n^2 x nbpop) matrix
ns=size(xt_ind,1)^2
dv_dbeta=zeros(ns,sum(globalStructure$notfixed_bpop))

k=1
if((globalStructure$m2_switch[1]==0 || globalStructure$m2_switch[1] == 1)){
#     if (globalStructure$m2_switch[1] == 0 && (globalStructure$gradff_switch[1]==0 || globalStructure$gradfg_switch[1] == 0 || globalStructure$hle_switch==0))
#         stop(sprintf('Complex differentiation to derive M2 cant be used with complex differentiaton of gradff, gradfg or hle\n'))
#     }
    
    for(i in 1:globalStructure$nbpop){
        if((globalStructure$notfixed_bpop[i]==1)){
            bpop_plus=bpop
            
            #if (globalStructure$m2_switch[1] == 1) #Central approximation or complex always uses central
                bpop_plus[i]=bpop_plus[i]+globalStructure$hm2
                bpop_minus=bpop
                bpop_minus[i]=bpop_minus[i]-globalStructure$hm2
                
                if((globalStructure$bCalculateEBE)){
                    #zeros(size(b_ind)[1],size(b_ind)[2])
                    start_bind = t(b_ind)
                    b_ind_plus = ind_estimates(globalStructure$mean_data,bpop_plus,d,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
                    b_ind_minus = ind_estimates(globalStructure$mean_data,bpop_minus,d,sigma,start_bind,(globalStructure$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,globalStructure)
                } else {
                    b_ind_plus = b_ind
                    b_ind_minus = b_ind
                }
                
                 returnArgs <-  v(model_switch,xt_ind,x,a,bpop_plus,b_ind_plus,bocc_ind,d,sigma,docc,globalStructure) 
v_plus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
                 returnArgs <-  v(model_switch,xt_ind,x,a,bpop_minus,b_ind_minus,bocc_ind,d,sigma,docc,globalStructure) 
v_minus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
                dv=v_plus-v_minus
                if((!isempty(dv))){
                    #Central
                    ir=dv/(2*globalStructure$hm2)
                    ir=reshape_matlab(ir,ns,1)
                    dv_dbeta[,k]=ir
                }
            #else #Complex differentiation
            #    stop(sprintf('Complex derivative cant be used to derive M2\n'))
#                 bpop_plus[i]=bpop_plus[i]+j*globalStructure$hm2
#                 [v_plus,globalStructure] = v(model_switch,xt_ind,x,a,bpop_plus,b_ind,bocc_ind,d,sigma,docc,globalStructure)
#                 dv=v_plus
#                 if (!isempty(dv))
#                     #Complex
#                     ir=Im(dv)/(globalStructure$hm2)
#                     ir=reshape_matlab(ir,ns,1)
#                     dv_dbeta[,k]=ir
#                 }
           # }
            k=k+1
        }
    }
} else {
    if((globalStructure$m2_switch[1]==20) ){#Analytic derivative
      stop("Analytic derivatives not implemented in PopED for R")
      #         dv_dbeta=analytic_dvary_dbpop[model_switch,xt_ind,x,a,bpop,d]
#         for(i in globalStructure$nbpop:-1:1){
#             if((globalStructure$notfixed_bpop[i]==0)){
#                 dv_dbeta[,i]=matrix(0,0,0)
#             }
#         }
    } else {
        if((globalStructure$m2_switch[1]==30) ){#Automatic differentiation (INTLab), only works with "normal" variance
          stop("Automatic differentiation not implemented in PopED for R")
#              returnArgs <-  m2_ad_wrapper(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,globalStructure) 
# dv_dbeta <- returnArgs[[1]]
# globalStructure <- returnArgs[[2]]
        } else {
         stop(sprintf('Unknown derivative option for m2'))
        }
    }
}
return(list( dv_dbeta= dv_dbeta,globalStructure=globalStructure)) 
}

  

