## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

mf_all <- function(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure){

    returnArgs <- switch(globalStructure$iFIMCalculationType+1,
                         mf(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure), #Default (with no assumption that bpop and b are uncorrelated)
                         mf3(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure), #Reduced FIM
                         stop("Not yet implemented"), #Weighted models
                         stop("Not yet implemented"), #Loc models
                         mf5(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure), #Reduced FIM with derivative of SD sigma
                         mf6(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure), #FULL FIM parameterized with A,B,C matrices & derivative of variance
                         mf7(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure), #Calculate one model switch at a time, good for large matrices
                         mf8(model_switch,xt_ind,x,a,bpop,d,sigma,docc,globalStructure)) #Reduced FIM parameterized with A,B,C matrices & derivative of variance
    
    if(is.null(returnArgs)) stop(sprintf('Unknown FIM-calculation type'))
    ret <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    return(list( ret= ret,globalStructure =globalStructure )) 
}


