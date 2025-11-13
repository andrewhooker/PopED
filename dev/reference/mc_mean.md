# Compute the monte-carlo mean of a function

Function computes the monte-carlo mean of a function by varying the
parameter inputs to the function

## Usage

``` r
mc_mean(
  ofv_fcn,
  poped.db,
  bpopdescr = poped.db$parameters$bpop,
  ddescr = poped.db$parameters$d,
  doccdescr = poped.db$parameters$d,
  user_distribution_pointer = poped.db$model$user_distribution_pointer,
  ED_samp_size = poped.db$settings$ED_samp_size,
  bLHS = poped.db$settings$bLHS,
  ...
)
```

## Arguments

- ofv_fcn:

  A function with poped.db as the first input

- poped.db:

  A PopED database.

- bpopdescr:

  Matrix defining the fixed effects, per row (row number =
  parameter_number) we should have:

  - column 1 the type of the distribution for E-family designs (0 =
    Fixed, 1 = Normal, 2 = Uniform, 3 = User Defined Distribution, 4 =
    lognormal and 5 = truncated normal)

  - column 2 defines the mean.

  - column 3 defines the variance of the distribution (or length of
    uniform distribution).

- ddescr:

  Matrix defining the diagonals of the IIV (same logic as for the
  `bpopdescr`).

- doccdescr:

  Matrix defining the IOV. per row (row number = parameter_number) we
  should have:

  - column 1 the type of the distribution for E-family designs (0 =
    Fixed, 1 = Normal, 2 = Uniform, 3 = User Defined Distribution, 4 =
    lognormal and 5 = truncated normal)

  - column 2 defines the mean of the variance.

  - column 3 defines the variance of the distribution (or length of
    uniform distribution).

- user_distribution_pointer:

  Function name for user defined distributions for E-family designs

- ED_samp_size:

  Sample size for E-family sampling

- bLHS:

  How to sample from distributions in E-family calculations. 0=Random
  Sampling, 1=LatinHyperCube –

- ...:

  Other arguments passed to the function.

## Value

The mean of the function evaluated at different parameter values.
