## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

get_fim_size <- function(globalStructure){
#Returns the size of FIM, i$e. col or row size
numnotfixed_bpop = sum(globalStructure$notfixed_bpop)
numnotfixed_d    = sum(globalStructure$notfixed_d)
numnotfixed_covd = sum(globalStructure$notfixed_covd)
numnotfixed_docc  = sum(globalStructure$notfixed_docc)
numnotfixed_covdocc  = sum(globalStructure$notfixed_covdocc)
numnotfixed_sigma  = sum(globalStructure$notfixed_sigma)
numnotfixed_covsigma  = sum(globalStructure$notfixed_covsigma)

n = numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma

return( n ) 
}
