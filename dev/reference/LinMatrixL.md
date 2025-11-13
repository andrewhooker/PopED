# The linearized matrix L

Function computes the derivative of the model with respect to the
between subject variability terms in the model (b's and bocc's)
evaluated at a defined point (b_ind and bocc_ind).

## Usage

``` r
LinMatrixL(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, poped.db)
```

## Arguments

- model_switch:

  A vector that is the same size as xt, specifying which model each
  sample belongs to.

- x:

  A vector for the discrete design variables.

- a:

  A vector of covariates.

- bpop:

  The fixed effects parameter values. Supplied as a vector.

- b_ind:

  The point at which to evaluate the derivative

- bocc_ind:

  The point at which to evaluate the derivative

- poped.db:

  A PopED database.

## Value

As a list:

- y:

  A matrix of size (samples per individual x number of random effects)

- poped.db:

  A PopED database

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


#for the FO approximation
ind=1
LinMatrixL(model_switch=t(poped.db$design$model_switch[ind,,drop=FALSE]),
          xt_ind=t(poped.db$design$xt[ind,,drop=FALSE]),
          x=zeros(0,1),
          a=t(poped.db$design$a[ind,,drop=FALSE]),
          bpop=poped.db$parameters$bpop[,2,drop=FALSE],
          b_ind=zeros(poped.db$parameters$NumRanEff,1),
          bocc_ind=zeros(poped.db$parameters$NumDocc,1),
          poped.db)["y"]
#> $y
#>             [,1]       [,2]        [,3]
#> [1,] -0.01736446 -3.4080713  2.63882264
#> [2,] -0.05954832 -5.4115558  3.17591023
#> [3,] -0.18102648 -7.2011569  2.27256206
#> [4,] -0.74460344 -7.2016770 -0.01921862
#> [5,] -2.44998833 -3.2358678 -0.10864692
#> [6,] -2.97791129 -1.5623370 -0.08675634
#> [7,] -3.07661786  0.7649213 -0.04417255
#> [8,] -2.09673861  1.1568729 -0.01795922
#> 
```
