# Compute the inverse of a matrix

Function computes the inverse of a matrix.

## Usage

``` r
inv(mat, method = 1, tol = .Machine$double.eps, pseudo_on_fail = TRUE, ...)
```

## Arguments

- mat:

  A matrix

- method:

  Which method to use. 1 is Cholesky `chol2inv(chol(mat)`, 2 is using
  `solve(mat)` and 3 is the Moore-Penrose generalized inverse
  (pseudoinverse).

- tol:

  The tolerance at which we should identify a singular value as zero
  (used in pseudoinverse calculation).

- pseudo_on_fail:

  If another method fails should the Moore-Penrose generalized inverse
  (pseudoinverse) be used?

- ...:

  Not used.

## Value

The inverse matrix
