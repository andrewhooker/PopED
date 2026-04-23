# Function written to match MATLAB's rand function

Generate random samples from a uniform distribution \[0,1\] and return
in matrix form

## Usage

``` r
rand(dim1, dim2 = NULL)
```

## Arguments

- dim1:

  The dimension of the matrix (if square), otherwise the number of rows.

- dim2:

  The number of columns, if different from the number of rows.

## Value

Matrix of random generated samples.

## See also

Other MATLAB:
[`cell()`](https://andrewhooker.github.io/PopED/dev/reference/cell.md),
[`diag_matlab()`](https://andrewhooker.github.io/PopED/dev/reference/diag_matlab.md),
[`feval()`](https://andrewhooker.github.io/PopED/dev/reference/feval.md),
[`fileparts()`](https://andrewhooker.github.io/PopED/dev/reference/fileparts.md),
[`isempty()`](https://andrewhooker.github.io/PopED/dev/reference/isempty.md),
[`ones()`](https://andrewhooker.github.io/PopED/dev/reference/ones.md),
[`randn()`](https://andrewhooker.github.io/PopED/dev/reference/randn.md),
[`size()`](https://andrewhooker.github.io/PopED/dev/reference/size.md),
[`tic()`](https://andrewhooker.github.io/PopED/dev/reference/tic.md),
[`toc()`](https://andrewhooker.github.io/PopED/dev/reference/toc.md),
[`zeros()`](https://andrewhooker.github.io/PopED/dev/reference/zeros.md)

## Examples

``` r
rand(2,3)
#>            [,1]      [,2]      [,3]
#> [1,] 0.78562035 0.3140058 0.8876423
#> [2,] 0.01472684 0.9888032 0.8638064

rand(5)
#>            [,1]      [,2]      [,3]       [,4]      [,5]
#> [1,] 0.69633302 0.6753241 0.6039208 0.57017268 0.9722578
#> [2,] 0.55576633 0.6834282 0.6109661 0.71155293 0.8258911
#> [3,] 0.50001558 0.9068502 0.4766922 0.06932691 0.5327581
#> [4,] 0.02953794 0.1530525 0.9027106 0.10538600 0.9491165
#> [5,] 0.17232315 0.6742126 0.6682546 0.28168070 0.9006530
```
