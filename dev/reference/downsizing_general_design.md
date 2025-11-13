# Downsize a general design to a specific design

Function takes a design with potentially empty design variables and
rescues the design so that a FIM can be calculated using
[`mftot`](https://andrewhooker.github.io/PopED/dev/reference/mftot.md).

## Usage

``` r
downsizing_general_design(poped.db)
```

## Arguments

- poped.db:

  A PopED database

## Value

A list containing:

- ni:

  A vector of the number of samples in each group.

- xt:

  A matrix of sample times. Each row is a vector of sample times for a
  group.

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

- x:

  A matrix for the discrete design variables. Each row is a group.

- a:

  A matrix of covariates. Each row is a group.

- bpop:

  A matrix of fixed effect parameter values.

## See also

Other poped_input:
[`convert_variables()`](https://andrewhooker.github.io/PopED/dev/reference/convert_variables.md),
[`create.poped.database()`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md),
[`create_design()`](https://andrewhooker.github.io/PopED/dev/reference/create_design.md),
[`create_design_space()`](https://andrewhooker.github.io/PopED/dev/reference/create_design_space.md),
[`poped.choose()`](https://andrewhooker.github.io/PopED/dev/reference/poped.choose.md)
