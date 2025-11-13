# Create a matrix of ones

Create a matrix of ones of size (dim1 x dim2).

## Usage

``` r
ones(dim1, dim2 = NULL)
```

## Arguments

- dim1:

  The dimension of the matrix (if square) or the number of rows.

- dim2:

  The number of columns

## Value

A matrix of ones

## See also

Other MATLAB:
[`cell()`](https://andrewhooker.github.io/PopED/dev/reference/cell.md),
[`diag_matlab()`](https://andrewhooker.github.io/PopED/dev/reference/diag_matlab.md),
[`feval()`](https://andrewhooker.github.io/PopED/dev/reference/feval.md),
[`fileparts()`](https://andrewhooker.github.io/PopED/dev/reference/fileparts.md),
[`isempty()`](https://andrewhooker.github.io/PopED/dev/reference/isempty.md),
[`rand()`](https://andrewhooker.github.io/PopED/dev/reference/rand.md),
[`randn()`](https://andrewhooker.github.io/PopED/dev/reference/randn.md),
[`size()`](https://andrewhooker.github.io/PopED/dev/reference/size.md),
[`tic()`](https://andrewhooker.github.io/PopED/dev/reference/tic.md),
[`toc()`](https://andrewhooker.github.io/PopED/dev/reference/toc.md),
[`zeros()`](https://andrewhooker.github.io/PopED/dev/reference/zeros.md)

## Examples

``` r
ones(4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    1    1    1
#> [2,]    1    1    1    1
#> [3,]    1    1    1    1
#> [4,]    1    1    1    1
ones(3,4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    1    1    1
#> [2,]    1    1    1    1
#> [3,]    1    1    1    1
```
