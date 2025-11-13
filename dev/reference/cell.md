# Create a cell array (a matrix of lists)

Create a cell array as in MATLAB.

## Usage

``` r
cell(...)
```

## Arguments

- ...:

  Dimensions for the cell array.

## Value

A list of empty lists.

## Note

This is a modified version of the same function in the matlab R-package.

## See also

Other MATLAB:
[`diag_matlab()`](https://andrewhooker.github.io/PopED/dev/reference/diag_matlab.md),
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
cell(3)
#>      [,1]      [,2]      [,3]     
#> [1,] numeric,0 numeric,0 numeric,0
#> [2,] numeric,0 numeric,0 numeric,0
#> [3,] numeric,0 numeric,0 numeric,0
cell(2,3)
#>      [,1]      [,2]      [,3]     
#> [1,] numeric,0 numeric,0 numeric,0
#> [2,] numeric,0 numeric,0 numeric,0

## define possible values of 2 categorical design variable
x.space <- cell(1,2)
x.space[1,1] <- list(seq(10,100,10))
x.space[1,2] <- list(seq(10,300,10))
x.space
#>      [,1]       [,2]      
#> [1,] numeric,10 numeric,30
x.space[1,1]
#> [[1]]
#>  [1]  10  20  30  40  50  60  70  80  90 100
#> 
x.space[1,2]
#> [[1]]
#>  [1]  10  20  30  40  50  60  70  80  90 100 110 120 130 140 150 160 170 180 190
#> [20] 200 210 220 230 240 250 260 270 280 290 300
#> 
```
