#' Downsize a general design to a specific design
#' 
#' Function takes a design with potentially empty design 
#' variables and resuces the design so that a FIM can be calculated.
#' 
#' @param globalStructure A PopED database 
#' @return A list containing:
#' \item{ni}{A vector of the number of samples in each group.}
#' \item{xt}{A matrix of sample times.  Each row is a vector of sample times for a group.}
#' \item{model_switch}{A matrix that is the same size as xt, specifying which model each sample belongs to.}
#' \item{x}{A matrix for the discrete design variables.  Each row is a group.}
#' \item{a}{A matrix of covariates.  Each row is a group.}
#' \item{bpop}{A matrix of fixed effect parameter values.}
#' 
#' @family poped_input


## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker



downsizing_general_design <- function(globalStructure){
# ------------- downsizing of general design

ni=globalStructure$gni[1:globalStructure$m,,drop=F]
xt=globalStructure$gxt[1:globalStructure$m,1:globalStructure$maxni,drop=F]
model_switch = globalStructure$global_model_switch[1:globalStructure$m,1:globalStructure$maxni,drop=F]

if((globalStructure$nx!=0)){
    x=globalStructure$gx[1:globalStructure$m,1:globalStructure$nx,drop=F]
} else {
    x = zeros(globalStructure$m,0)
}
if((globalStructure$na!=0)){
    a=globalStructure$ga[1:globalStructure$m,1:globalStructure$na,drop=F]
    maxa=globalStructure$gmaxa[1:globalStructure$m,1:globalStructure$na,drop=F]
    mina=globalStructure$gmina[1:globalStructure$m,1:globalStructure$na,drop=F]
} else {
    a = zeros(globalStructure$m,0)
    maxa = matrix(0,0,0)
    mina = matrix(0,0,0)
}
bpop=globalStructure$gbpop[1:globalStructure$nbpop,1:3,drop=F]
n=t(ni)%*%matrix(1,globalStructure$m,1)

d=globalStructure$gd[1:globalStructure$NumRanEff,1:3,drop=F]
maxxt=globalStructure$gmaxxt[1:globalStructure$m,1:globalStructure$maxni,drop=F]
minxt=globalStructure$gminxt[1:globalStructure$m,1:globalStructure$maxni,drop=F]

return(list( ni= ni, xt= xt, model_switch= model_switch, x= x, a= a, bpop= bpop, 
          n= 
          n, d= d, maxxt= maxxt, minxt= minxt,maxa=maxa,mina =mina )) 
}
