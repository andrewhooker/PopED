#' Create a matrix of zeros.
#' 
#' Create a matrix of zeros of size (dim1 x dim2). 
#' 
#' @param dim1 The dimension of the matrix (if square) or the number of rows.   
#' @param dim2 The number of columns 
#' @return A matrix of zeros.
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_zeros.R
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

zeros <- function(dim1,dim2=NULL){
  mat <- ones(dim1=dim1,dim2=dim2)*0
  return(mat)
}
