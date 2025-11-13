# Function written to match MATLAB's size function

Function written to match MATLAB's size function

## Usage

``` r
size(obj, dimension.index = NULL)
```

## Arguments

- obj:

  An object you want to know the various dimensions of. Typically a
  matrix.

- dimension.index:

  Which dimension you are interested in.

## Value

The dimensions of the object or specific dimension you are interested
in.

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
[`tic()`](https://andrewhooker.github.io/PopED/dev/reference/tic.md),
[`toc()`](https://andrewhooker.github.io/PopED/dev/reference/toc.md),
[`zeros()`](https://andrewhooker.github.io/PopED/dev/reference/zeros.md)

## Examples

``` r
size(c(2,3,4,5,6))
#> [1] 1 5

size(10)
#> [1] 1 1

size(zeros(4,7))
#> [1] 4 7
```
