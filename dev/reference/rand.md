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
#>           [,1]      [,2]      [,3]
#> [1,] 0.9403509 0.6145529 0.2367191
#> [2,] 0.1659312 0.7401846 0.7116930

rand(5)
#>             [,1]      [,2]       [,3]      [,4]      [,5]
#> [1,] 0.912880327 0.7414747 0.05472421 0.2370465 0.3529122
#> [2,] 0.009274727 0.3445153 0.70319566 0.1136393 0.7833523
#> [3,] 0.531906917 0.1167471 0.69419626 0.3805123 0.7474697
#> [4,] 0.376585201 0.6221801 0.52934335 0.8246649 0.2475224
#> [5,] 0.926520199 0.5967847 0.15631651 0.6082446 0.5962495
```
