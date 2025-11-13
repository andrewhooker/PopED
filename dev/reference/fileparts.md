# MATLAB fileparts function

Get the various parts of a file with path string.

## Usage

``` r
fileparts(filename.with.path)
```

## Arguments

- filename.with.path:

  A string of a filename with a path

## Value

A list with the following components:

- pathname:

  The path name

- filename:

  The file name

- fileext:

  The file extension

## Note

This is a modified version of the same function in the matlab R-package.

## See also

Other MATLAB:
[`cell()`](https://andrewhooker.github.io/PopED/dev/reference/cell.md),
[`diag_matlab()`](https://andrewhooker.github.io/PopED/dev/reference/diag_matlab.md),
[`feval()`](https://andrewhooker.github.io/PopED/dev/reference/feval.md),
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
fileparts("ggg/ttt/lll.R")
#> $pathname
#> [1] "ggg/ttt"
#> 
#> $filename
#> [1] "lll"
#> 
#> $fileext
#> [1] ".R"
#> 
```
