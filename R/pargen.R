#' Parameter simulation
#' 
#' Function generates random samples for a list of parameters
#' 
#' @param par A matrix describing the parameters. Each row is a parameter and 
#'   the matrix has three columns: 
#'   \enumerate{ 
#'   \item First column - Type of
#'   distribution (0-fixed, 1-normal, 2-uniform, 3-user specified, 4-lognormal,
#'   5-Truncated normal). 
#'   \item Second column - Mean of distribution. 
#'   \item Third
#'   column - Variance or range of distribution. 
#'   }
#' @param user_dist_pointer A text string of the name of a function that
#'   generates random samples from a user defined distribution.
#' @param sample_size The number of random samples per parameter to generate
#' @param bLHS Logical, indicating if Latin Hypercube Sampling should be used.
#' @param sample_number The sample number to extract.
#' @param globalStructure A PopED database.
#'   
#' @return A matrix of random samples of size (sample_size x
#'   number_of_parameters)
#'   

## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

pargen <- function(par,user_dist_pointer,sample_size,bLHS,sample_number,globalStructure){

nvar=size(par,1)
ret=zeros(sample_size,nvar)

if((bLHS==0) ){#Random Sampling
    for(k in 1:sample_size){
        np=size(par,1)
        if(np!=0){
            n=randn(np,1) # normal mean=0 and sd=1
            u=rand(np,1)*2-1 # uniform from -1 to 1
            t=par[,1] # type of distribution
            c2=par[,3] # variance or range of distribution

            bUserSpecifiedDistribution = (sum(t==3)>=1) #If at least one user specified distribution
            ret[k,] = (t==0)*par[,2,drop=F] + (t==2)*(par[,2,drop=F]+u*c2/2) + (t==1)*(par[,2,drop=F]+n*c2^(1/2))+(t==4)*par[,2,drop=F]*exp(n*c2^(1/2))
            if((sum(t==5)>0) ){#Truncated normal
                for(i in 1:size(par,1)){
                    if((t(i)==5)){
                        ret[k,i]=getTruncatedNormal(par[i,2],c2[i])
                    }
                }
            }

            if((bUserSpecifiedDistribution)){
                if((isempty(sample_number))){
                    ret[k,,drop=F] = feval(user_dist_pointer,ret[k,,drop=F],t,k,globalStructure)
                } else {
                    ret[k,,drop=F] = feval(user_dist_pointer,ret[k,,drop=F],t,sample_number,globalStructure)
                }
            }
        }
    }
} else if(nvar!=0){ #LHS
    ran=rand(sample_size,nvar)
    # method of Stein
    for(j in 1:nvar){
        idx=randperm(sample_size)
        P=(idx-ran[,j])/sample_size       # probability of the cdf
        returnArgs <- switch(par[j,1]+1,
                             par[j,2],  #point
                             par[j,2]+qnorm(P)*sqrt(par[j,3]), # normal
                             par[j,2]-par[j,3]/2 + P*par[j,3], #uniform
                             ret[,j], #Do nothing
                             par[j,2]*exp(qnorm(P)*sqrt(par[j,3])) #log-normal
                             )
        if(is.null(returnArgs)) stop(sprintf('Unknown distribution for the inverse probability function used in Latin Hypercube Sampling'))
        ret[,j]=returnArgs
    }
    
    bUserSpecifiedDistribution = (sum(par[,1,drop=F]==3)>=1) #If at least one user specified distribution
    
    if((bUserSpecifiedDistribution)){
        for(k in 1:sample_size){
            if((isempty(sample_number))){
                ret[k,,drop=F] = feval(user_dist_pointer,ret[k,,drop=F],par[,1,drop=F],k,globalStructure)
            } else {
                ret[k,,drop=F] = feval(user_dist_pointer,ret[k,,drop=F],par[,1,drop=F],sample_number,globalStructure)
            }
        }
    }
}
return( ret) 
}

#'  Generate a random sample from a truncated normal distribution.
#'  
#'  @param mean the mean of the normal distribution
#'  @param variance The variance of the normal distribution
#'    
#'  @return A random sample from the specified truncated normal distribution
#'    
getTruncatedNormal <- function(mean,variance){
while(TRUE){
    n = mean+randn(1,1)*sqrt(variance)
    if((sign(n)==sign(mean))){
        break
    }
}
return( n) 
}

#sign.2 <- function(x){
#  s = x/abs(x)*(x!=0)
#return( s) 
#}

