# Create a matrix of zeros.

Create a matrix of zeros of size (dim1 x dim2).

## Usage

``` r
zeros(dim1, dim2 = NULL)
```

## Arguments

- dim1:

  The dimension of the matrix (if square) or the number of rows.

- dim2:

  The number of columns

## Value

A matrix of zeros.

## See also

Other MATLAB:
[`cell()`](https://andrewhooker.github.io/PopED/dev/reference/cell.md),
[`diag_matlab()`](https://andrewhooker.github.io/PopED/dev/reference/diag_matlab.md),
[`feval()`](https://andrewhooker.github.io/PopED/dev/reference/feval.md),
[`fileparts()`](https://andrewhooker.github.io/PopED/dev/reference/fileparts.md),
[`isempty()`](https://andrewhooker.github.io/PopED/dev/reference/isempty.md),
[`ones()`](https://andrewhooker.github.io/PopED/dev/reference/ones.md),
[`rand()`](https://andrewhooker.github.io/PopED/dev/reference/rand.md),
[`randn()`](https://andrewhooker.github.io/PopED/dev/reference/randn.md),
[`size()`](https://andrewhooker.github.io/PopED/dev/reference/size.md),
[`tic()`](https://andrewhooker.github.io/PopED/dev/reference/tic.md),
[`toc()`](https://andrewhooker.github.io/PopED/dev/reference/toc.md)

## Examples

``` r
zeros(3)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
zeros(0,3)
#>      [,1] [,2] [,3]
zeros(4,7)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    0    0    0    0    0    0    0
#> [2,]    0    0    0    0    0    0    0
#> [3,]    0    0    0    0    0    0    0
#> [4,]    0    0    0    0    0    0    0
zeros(1,4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
```
