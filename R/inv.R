## Function written to match MATLAB function
## Author: Andrew Hooker

#' Compute the inverse of a matrix
#'
#' Function computes the inverse of a matrix.
#'
#' @param mat A matrix
#' @param method Which method to use. 1 is cholesky, 2 is using solve and 3 is the Moore-Penrose generalized inverse.
#' @param tol The tolerance at wich we should identify a sigular value as zero.
#'
#' @return The inverse matrix
#' @export
#' @keywords internal

inv<- function(mat,method=1,tol = .Machine$double.eps){
  if(method==1) return(chol2inv(chol(mat)))
  if(method==2) return(solve(mat))
  if(method==3){ # Moore-Penrose generalized inverse
    Xsvd <- svd(mat)
    #if (is.complex(mat)) Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * max(dim(mat))*max(Xsvd$d), 0)
    if (all(Positive)) return(Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u)))
    else if (!any(Positive)) return(array(0, dim(mat)[2L:1L]))
    else return(Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE])))
  }
}
