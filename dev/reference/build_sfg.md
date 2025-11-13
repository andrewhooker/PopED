# Build PopED parameter function from a model function

Build PopED parameter function from a model function

## Usage

``` r
build_sfg(
  model = "ff.PK.1.comp.oral.sd.CL",
  covariates = c("dose", "tau"),
  par_names = NULL,
  etas = "exp",
  no_etas = c("F", "Favail"),
  env = parent.frame()
)
```

## Arguments

- model:

  A string of text describing the model function name

- covariates:

  A list of covariate names to be filtered out of the model

- par_names:

  A list of parameter names in the model file. If not supplied then all
  undefined variables in the model file are extracted and the covariate
  names are filtered out of that list.

- etas:

  Can be "exp", "prop", "add" or "none". Either one value for all
  parameters or a list defining the model per parameter.

- no_etas:

  Parameters that should not have etas associated with them.

- env:

  The environment to create the function in.

## Value

A parameter model function to be used as input to PopED calculations.

## Examples

``` r
build_sfg(model="ff.PK.1.comp.oral.md.CL")
#> function (x, a, bpop, b, bocc) 
#> parameters = c(CL = bpop[1] * exp(b[1]), Favail = bpop[2], KA = bpop[3] * 
#>     exp(b[2]), V = bpop[4] * exp(b[3]), DOSE = a[1], TAU = a[2])
#> <environment: 0x56095b45c170>

etas <- c(Favail="exp",KA="exp",V="add",CL="exp")
build_sfg(model="ff.PK.1.comp.oral.md.CL",etas = etas)
#> function (x, a, bpop, b, bocc) 
#> parameters = c(CL = bpop[1] * exp(b[1]), Favail = bpop[2] * exp(b[2]), 
#>     KA = bpop[3] * exp(b[3]), V = bpop[4] + b[4], DOSE = a[1], 
#>     TAU = a[2])
#> <environment: 0x56095b45c170>
```
