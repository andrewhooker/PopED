# Timer function (as in MATLAB)

Function to stop a timer. Start with tic().

## Usage

``` r
toc(echo = TRUE, name = ".poped_savedTime")
```

## Arguments

- echo:

  Print time to screen?

- name:

  The saved name of the time object.

## Note

This is a modified version of the same function in the matlab R-package.

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
[`zeros()`](https://andrewhooker.github.io/PopED/dev/reference/zeros.md)

## Examples

``` r
tic()
toc()
#> Elapsed time: 0.001 seconds.

tic(name="foo")
toc()
#> Elapsed time: 0.002 seconds.
tic()
toc()
#> Elapsed time: 0 seconds.
toc()
#> Elapsed time: 0.001 seconds.
tic()
toc(name="foo")
#> Elapsed time: 0.002 seconds.
```
