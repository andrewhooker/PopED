gradfg_occ <- function(x,a,bpop,b_ind,bocc_ind,currentOcc,globalStructure){
#
#
# size: (number of g's x NumOccVariables)
#
# deriv of g's w$r.t. bocc's eval at b_ind, bocc_ind at occ currentOcc
#

dfg_db0=zeros(globalStructure$ng,size(bocc_ind,1))

#Central approximation
if((globalStructure$gradfg_switch[1] == 1)){
    for(k in 1:size(bocc_ind,1)){
        bocc_plus = bocc_ind
        bocc_minus = bocc_ind
        bocc_plus[k,currentOcc] = bocc_plus[k,currentOcc] + globalStructure$hlg
        bocc_minus[k,currentOcc] =bocc_minus[k,currentOcc] - globalStructure$hlg
        dfg_db0[,k]=(feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_plus)-feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_minus))/(2.0*globalStructure$hlg)
    }
} else {
    #Complex approximation
    if((globalStructure$gradfg_switch[1] == 0)){
        for(k in 1:size(bocc_ind,1)){
            bocc_plus = bocc_ind
            bocc_plus[k,currentOcc] = complex(real=bocc_plus[k,currentOcc], imaginary=globalStructure$hlg)
            dfg_db0[,k]=Im(feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_plus))/globalStructure$hlg
        }
    } else {
        #Calculate the analytic solution at b=b_ind
        if((globalStructure$gradfg_switch[1] == 20)){
            stop("Automatic calculation of analytic derivatives not currently implemented in PopED for R")
            #dfg_db0 = analytic_dfg_db1(x,a,bpop,b_ind)
        } else {
            if((globalStructure$gradfg_switch[1] == 30) ){#Automatic differentiation (INTLab)
              stop("Automatic differentiation not currently implemented in PopED for R")
#                 bocc_init = gradientinit(bocc_ind)
#                 val = globalStructure$fg_pointer(x,a,bpop,b_ind,bocc_init)
#                 dfg_db0 = val$dx(,(currentOcc-1)*globalStructure$NumDocc+1:currentOcc*globalStructure$NumDocc)
            } else {
                stop(sprintf('Unknown derivative option for gradfg_occ'))
            }
        }
    }
}
return( dfg_db0) 
}

