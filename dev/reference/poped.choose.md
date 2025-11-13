# Choose between `arg1` and `arg2`

Function chooses `arg1` unless it is `NULL` in which case `arg2` is
chosen.

## Usage

``` r
poped.choose(arg1, arg2)
```

## Arguments

- arg1:

  The first argument

- arg2:

  The second argument

## See also

Other poped_input:
[`convert_variables()`](https://andrewhooker.github.io/PopED/dev/reference/convert_variables.md),
[`create.poped.database()`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md),
[`create_design()`](https://andrewhooker.github.io/PopED/dev/reference/create_design.md),
[`create_design_space()`](https://andrewhooker.github.io/PopED/dev/reference/create_design_space.md),
[`downsizing_general_design()`](https://andrewhooker.github.io/PopED/dev/reference/downsizing_general_design.md)

## Examples

``` r
poped.choose(2,5)
#> [1] 2

poped.choose("foo",66)
#> [1] "foo"

poped.choose(NULL,"hello")
#> [1] "hello"
```
