# Evaluate a design

This function evaluates the design defined in a poped database.

## Usage

``` r
evaluate_design(poped.db, ...)
```

## Arguments

- poped.db:

  A poped database

- ...:

  Extra parameters passed to
  [`calc_ofv_and_fim`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_fim.md)
  and
  [`get_rse`](https://andrewhooker.github.io/PopED/dev/reference/get_rse.md)

## Value

A list of elements evaluating the current design.

## See also

Other evaluate_design:
[`evaluate.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md),
[`evaluate_power()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate_power.md),
[`get_rse()`](https://andrewhooker.github.io/PopED/dev/reference/get_rse.md),
[`model_prediction()`](https://andrewhooker.github.io/PopED/dev/reference/model_prediction.md),
[`plot_efficiency_of_windows()`](https://andrewhooker.github.io/PopED/dev/reference/plot_efficiency_of_windows.md),
[`plot_model_prediction()`](https://andrewhooker.github.io/PopED/dev/reference/plot_model_prediction.md)

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
#> <bytecode: 0x560952eadfa8>
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

evaluate_design(poped.db)
#> $ofv
#> [1] 52.44799
#> 
#> $fim
#>                   CL          V        KA         d_CL          d_V        d_KA
#> CL       19821.28445 -21.836551 -8.622140 0.000000e+00     0.000000  0.00000000
#> V          -21.83655  20.656071 -1.807099 0.000000e+00     0.000000  0.00000000
#> KA          -8.62214  -1.807099 51.729039 0.000000e+00     0.000000  0.00000000
#> d_CL         0.00000   0.000000  0.000000 3.107768e+03    10.728786  0.02613561
#> d_V          0.00000   0.000000  0.000000 1.072879e+01 27307.089308  3.26560786
#> d_KA         0.00000   0.000000  0.000000 2.613561e-02     3.265608 41.81083599
#> sig_prop     0.00000   0.000000  0.000000 5.215403e+02 11214.210707 71.08763902
#>              sig_prop
#> CL            0.00000
#> V             0.00000
#> KA            0.00000
#> d_CL        521.54030
#> d_V       11214.21071
#> d_KA         71.08764
#> sig_prop 806176.95068
#> 
#> $rse
#>        CL         V        KA      d_CL       d_V      d_KA  sig_prop 
#>  4.738266  2.756206 13.925829 25.627205 30.344316 25.777327 11.170784 
#> 

```
