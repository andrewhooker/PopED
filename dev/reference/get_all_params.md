# Extract all model parameters from the PopED database.

Extract all model parameters from the PopED database.

## Usage

``` r
get_all_params(poped.db)
```

## Arguments

- poped.db:

  A PopED database.

## Value

A list containing:

- bpop:

  A vector of fixed effect parameter values.

- d:

  A vector of between subject variability parameters

- covd:

  A vector of the covariances of the between subject variability
  parameters. Row major format of the lower triangular portion of the D
  (OMEGA) matrix

- docc:

  A vector of the between occasion variability (BOV) terms in the model

- covdocc:

  A vector of the covariances between the BOV terms. Row major of the
  lower triangular portion of the BOV matrix.

- sigma:

  A vector of the residual unexplained variances (RUV)

- covsigma:

  A vector of the covariances between the RUV terms

- all:

  A vector with all of the above, in the order of this list.

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


get_all_params(poped.db)
#> $bpop
#>        [,1]
#> CL     0.15
#> V      8.00
#> KA     1.00
#> Favail 1.00
#> 
#> $d
#>    [,1]
#> CL 0.07
#> V  0.02
#> KA 0.60
#> 
#> $covd
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> 
#> $docc
#>      [,1]
#> 
#> $covdocc
#>     
#> [1,]
#> 
#> $sigma
#> [1] 0.01
#> 
#> $covsigma
#>     
#> [1,]
#> 
#> $all
#>       [,1]
#>  [1,] 0.15
#>  [2,] 8.00
#>  [3,] 1.00
#>  [4,] 1.00
#>  [5,] 0.07
#>  [6,] 0.02
#>  [7,] 0.60
#>  [8,] 0.00
#>  [9,] 0.00
#> [10,] 0.00
#> [11,] 0.01
#> 

```
