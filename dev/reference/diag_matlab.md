# Function written to match MATLAB's diag function

There are some differences between tha MATLAB and the R version of diag.
Specifically, if a 1xN or a Nx1 matrix is supplied to the R
[`diag`](https://rdrr.io/r/base/diag.html) function then just the first
element of this vector is returned. This function tries to match the
MATLAB version in handling vectors (matricies with one dimension equal
to one), and will return a diagonal matrix in these situations.

## Usage

``` r
diag_matlab(mat)
```

## Arguments

- mat:

  Either a vector to make into a diagonal matrix or a matrix you want to
  extract the diagonal from

## Value

Either a diagonal matrix or the diagonal of a matrix.

## See also

Other MATLAB:
[`cell()`](https://andrewhooker.github.io/PopED/dev/reference/cell.md),
[`feval()`](https://andrewhooker.github.io/PopED/dev/reference/feval.md),
[`fileparts()`](https://andrewhooker.github.io/PopED/dev/reference/fileparts.md),
[`isempty()`](https://andrewhooker.github.io/PopED/dev/reference/isempty.md),
[`ones()`](https://andrewhooker.github.io/PopED/dev/reference/ones.md),
[`rand()`](https://andrewhooker.github.io/PopED/dev/reference/rand.md),
[`randn()`](https://andrewhooker.github.io/PopED/dev/reference/randn.md),
[`size()`](https://andrewhooker.github.io/PopED/dev/reference/size.md),
[`tic()`](https://andrewhooker.github.io/PopED/dev/reference/tic.md),
[`toc()`](https://andrewhooker.github.io/PopED/dev/reference/toc.md),
[`zeros()`](https://andrewhooker.github.io/PopED/dev/reference/zeros.md)

## Examples

``` r
diag_matlab(3)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
diag_matlab(c(1,2,3))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    2    0
#> [3,]    0    0    3
diag_matlab(cbind(1,2,3))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    2    0
#> [3,]    0    0    3
diag_matlab(rbind(1,2,3))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    2    0
#> [3,]    0    0    3

diag_matlab(matrix(c(1, 2, 3),6,6))
#> [1] 1 2 3 1 2 3

# here is where the R default does something different
diag(cbind(1,2,3))
#> [1] 1
diag(rbind(1,2,3))
#> [1] 1
```
