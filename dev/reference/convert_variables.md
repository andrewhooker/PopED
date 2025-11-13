# Create global variables in the PopED database

Function takes design variables from input files and converts them to
the global variables needed in PopED. Typically not used by the user.
Instead use the function
[`create.poped.database`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md).

## Usage

``` r
convert_variables(poped.db)
```

## Arguments

- poped.db:

  A PopED database

## Value

A PopED database

## See also

Other poped_input:
[`create.poped.database()`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md),
[`create_design()`](https://andrewhooker.github.io/PopED/dev/reference/create_design.md),
[`create_design_space()`](https://andrewhooker.github.io/PopED/dev/reference/create_design_space.md),
[`downsizing_general_design()`](https://andrewhooker.github.io/PopED/dev/reference/downsizing_general_design.md),
[`poped.choose()`](https://andrewhooker.github.io/PopED/dev/reference/poped.choose.md)
