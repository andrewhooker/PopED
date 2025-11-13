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
#>           [,1]        [,2]       [,3]
#> [1,] 0.6969792  1.00608416 -0.1798247
#> [2,] 0.3245009 -0.08181554  1.9912799

randn(5)
#>            [,1]       [,2]        [,3]       [,4]        [,5]
#> [1,]  0.4886010  0.7845366 -0.29945086 -0.6210085 -1.22803296
#> [2,] -0.2734432 -0.6285187 -1.13219142 -1.4092699 -0.44928916
#> [3,] -0.4182730 -0.1832261 -0.43740550  0.8429470 -0.71521453
#> [4,] -0.6120256  1.0192527 -1.43049739 -0.3497332  0.16371979
#> [5,]  1.3561098  0.2856346 -0.04297386  1.0899990 -0.02747351
```
