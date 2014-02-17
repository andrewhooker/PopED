## Function written to match MATLAB function
## Author: Andrew Hooker

inv<- function(mat){
    #return(solve(mat))
    return(chol2inv(chol(mat)))
}
