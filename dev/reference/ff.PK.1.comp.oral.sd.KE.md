# Structural model: one-compartment, oral absorption, single bolus dose, parameterized using KE.

This is a structural model function that encodes a model that is
one-compartment, oral absorption, single bolus dose, parameterized using
KE. The function is suitable for input to the
[`create.poped.database`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md)
function using the `ff_fun` or `ff_file` argument.

## Usage

``` r
ff.PK.1.comp.oral.sd.KE(model_switch, xt, parameters, poped.db)
```

## Arguments

- model_switch:

  a vector of values, the same size as `xt`, identifying which model
  response should be computed for the corresponding xt value. Used for
  multiple response models.

- xt:

  a vector of independent variable values (often time).

- parameters:

  A named list of parameter values.

- poped.db:

  a poped database. This can be used to extract information that may be
  needed in the model file.

## Value

A list consisting of:

1.  y the values of the model at the specified points.

2.  poped.db A (potentially modified) poped database.

## See also

Other models:
[`feps.add()`](https://andrewhooker.github.io/PopED/dev/reference/feps.add.md),
[`feps.add.prop()`](https://andrewhooker.github.io/PopED/dev/reference/feps.add.prop.md),
[`feps.prop()`](https://andrewhooker.github.io/PopED/dev/reference/feps.prop.md),
[`ff.PK.1.comp.oral.md.CL()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.md.CL.md),
[`ff.PK.1.comp.oral.md.KE()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.md.KE.md),
[`ff.PK.1.comp.oral.sd.CL()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.sd.CL.md),
[`ff.PKPD.1.comp.oral.md.CL.imax()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PKPD.1.comp.oral.md.CL.imax.md),
[`ff.PKPD.1.comp.sd.CL.emax()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PKPD.1.comp.sd.CL.emax.md)

Other structural_models:
[`ff.PK.1.comp.oral.md.CL()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.md.CL.md),
[`ff.PK.1.comp.oral.md.KE()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.md.KE.md),
[`ff.PK.1.comp.oral.sd.CL()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.sd.CL.md),
[`ff.PKPD.1.comp.oral.md.CL.imax()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PKPD.1.comp.oral.md.CL.imax.md),
[`ff.PKPD.1.comp.sd.CL.emax()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PKPD.1.comp.sd.CL.emax.md)

## Examples

``` r
library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.sd.KE
#> function (model_switch, xt, parameters, poped.db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         y = (DOSE * Favail * KA/(V * (KA - KE))) * (exp(-KE * 
#>             xt) - exp(-KA * xt))
#>         return(list(y = y, poped.db = poped.db))
#>     })
#> }
#> <bytecode: 0x564615ac9be0>
#> <environment: namespace:PopED>

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(KE=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.KE,
                                  fg_fun=sfg,
                                  fError_fun=feps.prop,
                                  bpop=c(KE=0.15/8, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(KE=0.07, V=0.02, KA=0.6), 
                                  sigma=0.01,
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70)

##  create plot of model without variability 
plot_model_prediction(poped.db)


## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
#>              [,1]       [,2]       [,3]         [,4]         [,5]       [,6]
#> [1,] 1248651.3457 339.370439 144.758375    0.0000000     0.000000  0.0000000
#> [2,]     339.3704  20.724252  -1.777254    0.0000000     0.000000  0.0000000
#> [3,]     144.7584  -1.777254  51.742071    0.0000000     0.000000  0.0000000
#> [4,]       0.0000   0.000000   0.000000 3010.9773834    40.490260  0.1151092
#> [5,]       0.0000   0.000000   0.000000   40.4902598 27487.656006  3.1586320
#> [6,]       0.0000   0.000000   0.000000    0.1151092     3.158632 41.8319046
#> [7,]       0.0000   0.000000   0.000000  784.2206814 10869.344969 70.0662365
#>              [,7]
#> [1,]      0.00000
#> [2,]      0.00000
#> [3,]      0.00000
#> [4,]    784.22068
#> [5,]  10869.34497
#> [6,]     70.06624
#> [7,] 807063.71904
det(FIM)
#> [1] 3.690649e+24
get_rse(FIM,poped.db)
#>         KE          V         KA       d_KE        d_V       d_KA SIGMA[1,1] 
#>   4.784638   2.756206  13.925829  26.037875  30.238758  25.770773  11.163214 
```
