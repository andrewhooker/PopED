# Compute an objective function and gradient

Compute an objective function and gradient with respect to the
optimization parameters. This function can be passed to the Broyden
Fletcher Goldfarb Shanno (BFGS) method for nonlinear minimization with
box constraints implemented in
[`bfgsb_min`](https://andrewhooker.github.io/PopED/dev/reference/bfgsb_min.md).

## Usage

``` r
calc_ofv_and_grad(
  x,
  optxt,
  opta,
  model_switch,
  aa,
  axt,
  groupsize,
  ni,
  xtopto,
  xopto,
  aopto,
  bpop,
  d,
  sigma,
  docc_full,
  poped.db,
  only_fim = FALSE
)
```

## Arguments

- x:

  A matrix for the discrete design variables. Each row is a group.

- optxt:

  If sampling times are optimized

- opta:

  If continuous design variables are optimized

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

- aa:

  The aa value

- axt:

  the axt value

- groupsize:

  A vector of the number of individuals in each group.

- ni:

  A vector of the number of samples in each group.

- xtopto:

  the xtopto value

- xopto:

  the xopto value

- aopto:

  the aopto value

- bpop:

  Matrix defining the fixed effects, per row (row number =
  parameter_number) we should have:

  - column 1 the type of the distribution for E-family designs (0 =
    Fixed, 1 = Normal, 2 = Uniform, 3 = User Defined Distribution, 4 =
    lognormal and 5 = truncated normal)

  - column 2 defines the mean.

  - column 3 defines the variance of the distribution (or length of
    uniform distribution).

  Can also just supply the parameter values as a vector
  [`c()`](https://rdrr.io/r/base/c.html) if no uncertainty around the
  parameter value is to be used. The parameter order of 'bpop' is
  defined in the 'fg_fun' or 'fg_file'. If you use named arguments in
  'bpop' then the order of this vector can be rearranged to match the
  'fg_fun' or 'fg_file'. See \`reorder_parameter_vectors\`.

- d:

  Matrix defining the diagonals of the IIV (same logic as for the fixed
  effects matrix bpop to define uncertainty). One can also just supply
  the parameter values as a [`c()`](https://rdrr.io/r/base/c.html). The
  parameter order of 'd' is defined in the 'fg_fun' or 'fg_file'. If you
  use named arguments in 'd' then the order of this vector can be
  rearranged to match the 'fg_fun' or 'fg_file'. See
  \`reorder_parameter_vectors\`.

- sigma:

  Matrix defining the variances can covariances of the residual
  variability terms of the model. can also just supply the diagonal
  parameter values (variances) as a
  [`c()`](https://rdrr.io/r/base/c.html).

- docc_full:

  A between occasion variability matrix.

- poped.db:

  A PopED database.

- only_fim:

  Should the gradient be calculated?

## Value

A list containing:

- f:

  The objective function.

- g:

  The gradient.

## See also

Other Optimize:
[`Doptim()`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md),
[`LEDoptim()`](https://andrewhooker.github.io/PopED/dev/reference/LEDoptim.md),
[`RS_opt()`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md),
[`a_line_search()`](https://andrewhooker.github.io/PopED/dev/reference/a_line_search.md),
[`bfgsb_min()`](https://andrewhooker.github.io/PopED/dev/reference/bfgsb_min.md),
[`calc_autofocus()`](https://andrewhooker.github.io/PopED/dev/reference/calc_autofocus.md),
[`mfea()`](https://andrewhooker.github.io/PopED/dev/reference/mfea.md),
[`optim_ARS()`](https://andrewhooker.github.io/PopED/dev/reference/optim_ARS.md),
[`optim_LS()`](https://andrewhooker.github.io/PopED/dev/reference/optim_LS.md),
[`poped_optim()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim.md),
[`poped_optim_1()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_1.md),
[`poped_optim_2()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_2.md),
[`poped_optim_3()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_3.md),
[`poped_optimize()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optimize.md)

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


opta=TRUE
aa=opta*poped.db$settings$cfaa*matrix(1,poped.db$design$m,size(poped.db$design$a,2))
aa
#>       [,1]
#> [1,] 0.001

optxt=TRUE
axt=optxt*poped.db$settings$cfaxt*matrix(1,poped.db$design$m,max(poped.db$design_space$maxni))
axt
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]
#> [1,] 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001

calc_ofv_and_grad(x=c(poped.db$design$xt,poped.db$design$a),
                  optxt=optxt, opta=opta, 
                  model_switch=poped.db$design$model_switch,
                  aa=aa,
                  axt=axt,
                  groupsize=poped.db$design$groupsize,
                  ni=poped.db$design$ni,
                  xtopto=poped.db$design$xt,
                  xopto=poped.db$design$x,
                  aopto=poped.db$design$a,
                  bpop=poped.db$parameters$param.pt.val$bpop,
                  d=poped.db$parameters$param.pt.val$d,
                  sigma=poped.db$parameters$param.pt.val$sigma,
                  docc_full=poped.db$parameters$param.pt.val$docc,
                  poped.db,
                  only_fim=FALSE)
#> $f
#> [1] -55.39645
#> 
#> $g
#>               [,1]
#>  [1,]  0.079631185
#>  [2,] -0.025697327
#>  [3,] -0.312973785
#>  [4,]  0.009484133
#>  [5,]  0.019485307
#>  [6,] -0.001675713
#>  [7,] -0.009544061
#>  [8,] -0.001929802
#>  [9,] -0.035844547
#> 

if (FALSE) { # \dontrun{
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)
  
} # }



```
