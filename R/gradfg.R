## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradfg <- function(x,a,bpop,b_ind,bocc_ind,globalStructure){
#
#
# size: (number of g's x number of random effects)
#
# deriv of g's w$r.t. b's and bocc's eval at b_ind
#

dfg_db0=zeros(globalStructure$ng,globalStructure$NumRanEff)

#Central approximation
if((globalStructure$gradfg_switch[1] == 1)){
    for(i in 1:globalStructure$NumRanEff){
        b_plus = b_ind
        b_minus = b_ind
        b_plus[i]=b_plus[i] + globalStructure$hlg
        b_minus[i]=b_minus[i]-globalStructure$hlg
        dfg_db0[,i]=(feval(globalStructure$fg_pointer,x,a,bpop,b_plus,bocc_ind)-feval(globalStructure$fg_pointer,x,a,bpop,b_minus,bocc_ind))/(2.0*globalStructure$hlg)
    }
} else {
    #Complex approximation
    if((globalStructure$gradfg_switch[1] == 0)){
        for(i in 1:globalStructure$NumRanEff){
            b_plus = b_ind
b_plus[i] = complex(real=b_plus[i],imaginary=globalStructure$hlg)
            dfg_db0[,i]=Im(feval(globalStructure$fg_pointer,x,a,bpop,b_plus,bocc_ind))/globalStructure$hlg
        }
    } else {
        #Calculate the analytic solution at b=b_ind
        if((globalStructure$gradfg_switch[1] == 20)){
          stop("Automatic computation of analytic derivatives not currently implemented in PopED for R")
          #dfg_db0 = analytic_dfg_db1(x,a,bpop,b_ind)
        } else {
            if((globalStructure$gradfg_switch[1] == 30) ){#Calculate using automatic differentiation (INTLab)
              stop("Automatic differentiation not currently implemented in PopED for R")
              #                 b_init = gradientinit(b_ind)
#                 val = feval(globalStructure$fg_pointer,x,a,bpop,b_init,bocc_ind)
#                 dfg_db0 = val$dx
            } else {
                stop(sprintf('Unknown derivative option for gradfg'))
            }
        }
    }
}
return( dfg_db0) 
}

