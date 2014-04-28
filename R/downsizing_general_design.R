#' Downsize a general design to a specific design
#' 
#' Function takes a design with potentially empty design 
#' variables and resuces the design so that a FIM can be calculated using \code{\link{mftot}}.
#' 
#' @param poped.db A PopED database 
#' @return A list containing:
#' \item{ni}{A vector of the number of samples in each group.}
#' \item{xt}{A matrix of sample times.  Each row is a vector of sample times for a group.}
#' \item{model_switch}{A matrix that is the same size as xt, specifying which model each sample belongs to.}
#' \item{x}{A matrix for the discrete design variables.  Each row is a group.}
#' \item{a}{A matrix of covariates.  Each row is a group.}
#' \item{bpop}{A matrix of fixed effect parameter values.}
#' 
#' @family poped_input
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_downsizing_general_design.R

## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker



downsizing_general_design <- function(poped.db){
# ------------- downsizing of general design

ni=poped.db$gni[1:poped.db$m,,drop=F]
xt=poped.db$gxt[1:poped.db$m,1:poped.db$maxni,drop=F]
model_switch = poped.db$global_model_switch[1:poped.db$m,1:poped.db$maxni,drop=F]

if((poped.db$nx!=0)){
    x=poped.db$gx[1:poped.db$m,1:poped.db$nx,drop=F]
} else {
    x = zeros(poped.db$m,0)
}
if((poped.db$na!=0)){
    a=poped.db$ga[1:poped.db$m,1:poped.db$na,drop=F]
    maxa=poped.db$gmaxa[1:poped.db$m,1:poped.db$na,drop=F]
    mina=poped.db$gmina[1:poped.db$m,1:poped.db$na,drop=F]
} else {
    a = zeros(poped.db$m,0)
    maxa = matrix(0,0,0)
    mina = matrix(0,0,0)
}
bpop=poped.db$gbpop[1:poped.db$nbpop,1:3,drop=F]
n=t(ni)%*%matrix(1,poped.db$m,1)

d=poped.db$gd[1:poped.db$NumRanEff,1:3,drop=F]
maxxt=poped.db$gmaxxt[1:poped.db$m,1:poped.db$maxni,drop=F]
minxt=poped.db$gminxt[1:poped.db$m,1:poped.db$maxni,drop=F]

return(list( ni= ni, xt= xt, model_switch= model_switch, x= x, a= a, bpop= bpop, 
          n= 
          n, d= d, maxxt= maxxt, minxt= minxt,maxa=maxa,mina =mina )) 
}
