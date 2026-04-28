# Model linearization with respect to epsilon and eta.

The function performs a linearization of the model with respect to the
residual variability and then the between subject variability.
Derivative of model w.r.t. eps then eta, evaluated at eps=0 and b=b_ind.

## Usage

``` r
LinMatrixLH(
  model_switch,
  xt_ind,
  x,
  a,
  bpop,
  b_ind,
  bocc_ind,
  NumEPS,
  poped.db
)
```

## Arguments

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

- xt_ind:

  A vector of the individual/group sample times

- x:

  A matrix for the discrete design variables. Each row is a group.

- a:

  A matrix of covariates. Each row is a group.

- bpop:

  The fixed effects parameter values. Supplied as a vector.

- b_ind:

  vector of individual realization of the BSV terms b

- bocc_ind:

  Vector of individual realizations of the BOV terms bocc

- NumEPS:

  The number of eps() terms in the model.

- poped.db:

  A PopED database.

## Value

A matrix of size (samples per individual x (number of sigma x number of
omega))

## See also

Other FIM:
[`LinMatrixH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixH.md),
[`LinMatrixL_occ()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixL_occ.md),
[`calc_ofv_and_fim()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_fim.md),
[`ed_laplace_ofv()`](https://andrewhooker.github.io/PopED/dev/reference/ed_laplace_ofv.md),
[`ed_mftot()`](https://andrewhooker.github.io/PopED/dev/reference/ed_mftot.md),
[`efficiency()`](https://andrewhooker.github.io/PopED/dev/reference/efficiency.md),
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md),
[`evaluate.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md),
[`gradf_eps()`](https://andrewhooker.github.io/PopED/dev/reference/gradf_eps.md),
[`mf3()`](https://andrewhooker.github.io/PopED/dev/reference/mf3.md),
[`mf7()`](https://andrewhooker.github.io/PopED/dev/reference/mf7.md),
[`mftot()`](https://andrewhooker.github.io/PopED/dev/reference/mftot.md),
[`ofv_criterion()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_criterion.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)

## Examples

``` r
library(PopED)

############# START #################
## Create PopED database
## (warfarin model for optimization)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error  
## to avoid sample times at very low concentrations (time 0 or very late samples).

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
#> <bytecode: 0x56460e6c7a50>
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

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.add.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(prop=0.01,add=0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0.01,
                                  maxxt=120,
                                  a=c(DOSE=70),
                                  mina=c(DOSE=0.01),
                                  maxa=c(DOSE=100))

############# END ###################
## Create PopED database
## (warfarin model for optimization)
#####################################


#for the FOI approximation
ind=1
poped.db$settings$iApproximationMethod=3 # FOI approximation method

LinMatrixLH(model_switch=t(poped.db$design$model_switch[ind,,drop=FALSE]),
          xt_ind=t(poped.db$design$xt[ind,,drop=FALSE]),
          x=zeros(0,1),
          a=t(poped.db$design$a[ind,,drop=FALSE]),
          bpop=poped.db$parameters$bpop[,2,drop=FALSE],
          b_ind=zeros(poped.db$parameters$NumRanEff,1),
          bocc_ind=zeros(poped.db$parameters$NumDocc,1),
          NumEPS=size(poped.db$parameters$sigma,1),
          poped.db)["y"]
#> $y
#>             [,1] [,2]       [,3] [,4]        [,5] [,6]
#> [1,] -0.01736389    0 -3.4080716    0  2.63882249    0
#> [2,] -0.05954792    0 -5.4115556    0  3.17591065    0
#> [3,] -0.18102853    0 -7.2011552    0  2.27256214    0
#> [4,] -0.74460438    0 -7.2016748    0 -0.01922018    0
#> [5,] -2.44998688    0 -3.2358671    0 -0.10864643    0
#> [6,] -2.97791125    0 -1.5623369    0 -0.08675727    0
#> [7,] -3.07661785    0  0.7649203    0 -0.04417244    0
#> [8,] -2.09673889    0  1.1568727    0 -0.01795897    0
#> 

  
```
