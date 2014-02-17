## Function written to match MATLAB function
## Author: Andrew Hooker

ones <- function(dim1,dim2=NULL){
    if(is.null(dim2)) dim2 <- dim1
    mat <- matrix(1,dim1,dim2)
    return(mat)
}
