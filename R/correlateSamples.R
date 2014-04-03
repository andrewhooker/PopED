#' Correlate samples
#' 
#' @param R The R matrix
#' @param Omega The Omega matrix
correlateSamples <- function(R,Omega){

if(((size(Omega,1)==1 && size(Omega,2)==1) || !any(any(Omega-diag_matlab(diag_matlab(Omega)))))){
    ret = t(R)
    return()
}

std=sqrt(diag_matlab(Omega))
C=Omega/(std%*%t(std))
R=t(R)
#T

T=cor(R)              #correlation matrix of R
P=t(chol(C))     #Cholesky factorization of omega
Q=t(chol(T))          #Cholesky factorization of T
S=P/Q                 #calc S s$t. R*S=Omega
Rs=R*t(S)

ret=zeros(size(R))                  #allocate space for return values
for(j in 1:length(Omega)){
    r=rank(Rs[,j])       #determine rank of Rs entries
    RS=sort(R[,j])     #sort original entries of R
    ret[,j] = RS[r]       #permute entries of R according to Rs
}
return( ret) 
}

# ranking <- function(x){                                                                                                                                              
#  returnArgs <- sort(x)                                                              
# s <- returnArgs[[1]]
# i <- returnArgs[[2]]
# r[i,1]=t(1:length(x))
# return( r) 
# }

## %This function fixes a bug in FreeMat 4.0
## function ret=diag(a)
##     if (~isempty(a) && size(a,1)==1 && size(a,2)==1)
##         ret=builtin('diag',[a]);
##     else
##         ret=builtin('diag',a);
##     end
## end


