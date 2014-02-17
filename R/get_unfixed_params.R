## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

get_unfixed_params <- function(globalStructure,params=NULL){

#Return all the unfixed parameters with the specified order
#(bpop,d,covd,docc,covdocc,sigma,covsigma) all = vector of all unfixed
#params
#var derivative is a vector of 1 and 0, 1 means derivative of parameter is
#taken w$r.t. variance otherwise w$r.t. sd
# If params is supplied then the parameter is taken from this vector 
# instead of globalStructure

if(is.null(params)){
    bpop = globalStructure$gbpop[,2,drop=F]
    d = globalStructure$gd[,2,drop=F]
    covd = globalStructure$covd
    docc = globalStructure$docc[,2,drop=F]
    covdocc = globalStructure$covdocc
    sigma = diag_matlab(globalStructure$sigma)
    covsigma = zeros(1,length(sigma)*(length(sigma)-1)/2)
    k=1
    for(i in 1:size(globalStructure$sigma,1)){
        for(j in 1:size(globalStructure$sigma,2)){
            if((i<j)){
                covsigma[k] = globalStructure$sigma[i,j]
                k=k+1
            }
        }
    }
} else {
    nbpop = length(globalStructure$notfixed_bpop)
    nd = length(globalStructure$notfixed_d)
    ncovd = length(globalStructure$notfixed_covd)
    ndocc = length(globalStructure$notfixed_docc)
    ncovdocc = length(globalStructure$notfixed_covdocc)
    nsigma = length(globalStructure$notfixed_sigma)
    ncovsigma = length(globalStructure$notfixed_covsigma)

    bpop = params[1:nbpop]
    d=params[(1+nbpop):(nbpop+nd)]
    covd=params[(1+nbpop+nd):(nbpop+nd+ncovd)]
    docc=params[(1+nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc)]
    covdocc=params[(1+nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc)]
    sigma=params[(1+nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma)]
    covsigma=params[(1+nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma)]
}

bpop = bpop[globalStructure$notfixed_bpop==1]
d = d[globalStructure$notfixed_d==1]
covd = covd[globalStructure$notfixed_covd==1]
docc = docc[globalStructure$notfixed_docc==1]
covdocc = covdocc[globalStructure$notfixed_covdocc==1]
sigma = sigma[globalStructure$notfixed_sigma==1]
covsigma = covsigma[globalStructure$notfixed_covsigma==1]

all = matrix(c(bpop, d, covd, docc, covdocc, sigma, covsigma),ncol=1,byrow=T)

if((globalStructure$iFIMCalculationType!=4)){
    var_derivative = matrix(1,size(all))
} else {
    var_derivative = matrix(c(rep(1,length(bpop)), rep(1,length(d)), rep(1,length(covd)), rep(1,length(docc)), rep(1,length(covdocc)), rep(0,length(sigma)), rep(1,length(covsigma))),ncol=1,byrow=T)
}

return(list( bpop= bpop,d=d,covd=covd,docc=docc,covdocc=covdocc,sigma=sigma,covsigma=covsigma,all=all,var_derivative =var_derivative )) 
}

