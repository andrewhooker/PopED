# RUV model: Proportional.

This is a residual unexplained variability (RUV) model function that
encodes the model described above. The function is suitable for input to
the
[`create.poped.database`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md)
function using the `fError_file` argument.

## Usage

``` r
feps.prop(model_switch, xt, parameters, epsi, poped.db)
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

- epsi:

  A matrix with the same number of rows as the `xt` vector, columns
  match the numbers defined in this function.

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
[`ff.PK.1.comp.oral.md.CL()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.md.CL.md),
[`ff.PK.1.comp.oral.md.KE()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.md.KE.md),
[`ff.PK.1.comp.oral.sd.CL()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.sd.CL.md),
[`ff.PK.1.comp.oral.sd.KE()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PK.1.comp.oral.sd.KE.md),
[`ff.PKPD.1.comp.oral.md.CL.imax()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PKPD.1.comp.oral.md.CL.imax.md),
[`ff.PKPD.1.comp.sd.CL.emax()`](https://andrewhooker.github.io/PopED/dev/reference/ff.PKPD.1.comp.sd.CL.emax.md)

Other RUV_models:
[`feps.add()`](https://andrewhooker.github.io/PopED/dev/reference/feps.add.md),
[`feps.add.prop()`](https://andrewhooker.github.io/PopED/dev/reference/feps.add.prop.md)

## Examples

``` r
library(PopED)

############# START #################
## Create PopED database
## (warfarin example)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.sd.CL
#> function (model_switch, xt, parameters, poped.db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         y = (DOSE * Favail * KA/(V * (KA - CL/V))) * (exp(-CL/V * 
#>             xt) - exp(-KA * xt))
#>         return(list(y = y, poped.db = poped.db))
#>     })
#> }
#> <bytecode: 0x55c6a83a2c18>
#> <environment: namespace:PopED>

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

## -- Define model, parameters, initial design
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(prop=0.01),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  a=c(DOSE=70))

############# END ###################
## Create PopED database
## (warfarin example)
#####################################


##  create plot of model without variability 
plot_model_prediction(poped.db)


## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
#>             [,1]       [,2]      [,3]         [,4]         [,5]        [,6]
#> [1,] 19821.28445 -21.836551 -8.622140 0.000000e+00     0.000000  0.00000000
#> [2,]   -21.83655  20.656071 -1.807099 0.000000e+00     0.000000  0.00000000
#> [3,]    -8.62214  -1.807099 51.729039 0.000000e+00     0.000000  0.00000000
#> [4,]     0.00000   0.000000  0.000000 3.107768e+03    10.728786  0.02613561
#> [5,]     0.00000   0.000000  0.000000 1.072879e+01 27307.089308  3.26560786
#> [6,]     0.00000   0.000000  0.000000 2.613561e-02     3.265608 41.81083599
#> [7,]     0.00000   0.000000  0.000000 5.215403e+02 11214.210707 71.08763902
#>              [,7]
#> [1,]      0.00000
#> [2,]      0.00000
#> [3,]      0.00000
#> [4,]    521.54030
#> [5,]  11214.21071
#> [6,]     71.08764
#> [7,] 806176.95068
det(FIM)
#> [1] 5.996147e+22
get_rse(FIM,poped.db)
#>        CL         V        KA      d_CL       d_V      d_KA  sig_prop 
#>  4.738266  2.756206 13.925829 25.627205 30.344316 25.777327 11.170784 
```
