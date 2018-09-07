## Function written to match MATLAB function
## Author: Andrew Hooker

#' Compute the inverse of a matrix
#'
#' Function computes the inverse of a matrix.
#'
#' @param mat A matrix
#' @param method Which method to use. 1 is Cholesky \code{chol2inv(chol(mat)}, 
#'   2 is using \code{solve(mat)} and 3 is the Moore-Penrose generalized inverse (pseudoinverse).
#' @param tol The tolerance at which we should identify a singular value as zero (used in pseudoinverse calculation).
#' @param pseudo_on_fail If another method fails should the Moore-Penrose generalized inverse (pseudoinverse) be used?
#' @param ... Not used.
#'
#' @return The inverse matrix
#' @export
#' @keywords internal

inv<- function(mat, method=1, tol = .Machine$double.eps, pseudo_on_fail = TRUE,...){
  
  invar <- NULL
  
  # Cholesky Decomposition costs approx. 50% more than calculating the determinant
  # So instead of testing the determinant we can immediately test the decomposition
  if(method==1){
    invvar <- try(chol2inv(chol(mat)), silent = TRUE)

    # If the result is numeric, it was successful
    if (is.numeric(invvar)) return(invvar)
  }
  if(method==2){
    invvar <- try(solve(mat), silent = TRUE)
    
    # If the result is numeric, it was successful
    if (is.numeric(invvar)) return(invvar)
  } 
  
  # Otherwise calculate the Moore-Penrose generalized inverse (pseudoinverse)
  if(method==3 || (pseudo_on_fail && !is.numeric(invvar))){
    message("Problems inverting the matrix. Results could be misleading.")
    Xsvd <- svd(mat)
    #if (is.complex(mat)) Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * max(dim(mat))*max(Xsvd$d), 0)
    if (all(Positive)) return(Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u)))
    else if (!any(Positive)) return(array(0, dim(mat)[2L:1L]))
    else return(Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE])))
  }  
  
  stop(invvar)
  
}
