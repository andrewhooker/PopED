# Function written to match MATLAB's randn function

Generate random samples from a standardized normal distribution and
return in matrix form.

## Usage

``` r
randn(dim1, dim2 = NULL)
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
[`rand()`](https://andrewhooker.github.io/PopED/dev/reference/rand.md),
[`size()`](https://andrewhooker.github.io/PopED/dev/reference/size.md),
[`tic()`](https://andrewhooker.github.io/PopED/dev/reference/tic.md),
[`toc()`](https://andrewhooker.github.io/PopED/dev/reference/toc.md),
[`zeros()`](https://andrewhooker.github.io/PopED/dev/reference/zeros.md)

## Examples

``` r
randn(2,3)
#>            [,1]      [,2]       [,3]
#> [1,] -0.1056558  1.613284  0.6253141
#> [2,]  1.4898955 -1.623811 -0.7983959

randn(5)
#>             [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] -2.00454302 -0.8669026  2.2927026 -0.2468592 -0.8973077
#> [2,]  1.41783276  1.8538769 -1.9795730 -1.0357461 -0.4940153
#> [3,]  1.21293195  0.9119217  0.9319057 -1.1906234  1.3887291
#> [4,]  0.99323055 -0.7318803  0.3306204 -1.3221883  1.7442997
#> [5,] -0.04224099 -0.2358161  0.9562781  0.3576470  0.2759333
```
