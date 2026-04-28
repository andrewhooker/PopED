# Display a summary of output from poped_optim

Display a summary of output from poped_optim

## Usage

``` r
# S3 method for class 'poped_optim'
summary(object, ...)
```

## Arguments

- object:

  An object returned from
  [`poped_optim`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim.md)
  to summarize.

- ...:

  Additional arguments. Passed to
  [`blockfinal`](https://andrewhooker.github.io/PopED/dev/reference/blockfinal.md).

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
#> <bytecode: 0x55bed4d97ed8>
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

##############
# D-family Optimization
##############


# ARS+BFGS+LS optimization of dose
# optimization with just a few iterations
# only to check that things are working
out_1 <- poped_optim(poped.db,opt_a =TRUE,
                      control = list(ARS=list(iter=2),
                                     BFGS=list(maxit=2),
                                     LS=list(line_length=2)),
                      iter_max = 1)
#> ===============================================================================
#> Initial design evaluation
#> 
#> Initial OFV = 55.3964
#> 
#> Initial design
#> expected relative standard error
#> (%RSE, rounded to nearest integer)
#>    Parameter   Values   RSE_0
#>           CL     0.15       5
#>            V        8       3
#>           KA        1      14
#>         d_CL     0.07      30
#>          d_V     0.02      37
#>         d_KA      0.6      27
#>     sig_prop     0.01      32
#>      sig_add     0.25      26
#> 
#> ==============================================================================
#> Optimization of design parameters
#> 
#> * Optimize Covariates
#> 
#> ************* Iteration 1  for all optimization methods***********************
#> 
#> *******************************************
#> Running Adaptive Random Search Optimization
#> *******************************************
#> Initial OFV = 55.3964
#> 
#> Total iterations: 2 
#> Elapsed time: 0.022 seconds.
#> 
#> Final OFV =  56.01888 
#> Parameters: 98.84467 
#> 
#> *******************************************
#> Running BFGS Optimization
#> *******************************************
#> initial  value -56.018883 
#> final  value -56.019059 
#> stopped after 2 iterations
#> 
#> *******************************************
#> Running Line Search Optimization
#> *******************************************
#> 
#>    Initial parameters: 98.85978 
#>    Initial OFV: 56.01906 
#> 
#>    Searching parameter 1 
#>      Changed from 98.8598 to 100 ; OFV = 56.032 
#> 
#>    Elapsed time: 0.032 seconds.
#> 
#>    Final OFV =  56.03204 
#>    Parameters: 100 
#> 
#> *******************************************
#> Stopping criteria testing
#> (Compare between start of iteration and end of iteration)
#> *******************************************
#> Difference in OFV:  0.636
#> Relative difference in OFV:  1.15%
#> Efficiency: 
#>   ((exp(ofv_final) / exp(ofv_init))^(1/n_parameters)) = 1.0827
#> 
#>  Efficiency stopping criteria: 
#>   Is (1.0827 <= 1.001)?   No.
#>   Stopping criteria NOT achieved.
#> 
#> Stopping criteria NOT achieved.
#> 
#> ===============================================================================
#> FINAL RESULTS
#> 
#> Optimized Covariates:
#> Group 1: 100
#> 
#> OFV = 56.032
#> 
#> Efficiency: 
#>   ((exp(ofv_final) / exp(ofv_init))^(1/n_parameters)) = 1.0827
#> 
#> Expected relative standard error
#> (%RSE, rounded to nearest integer):
#>    Parameter   Values   RSE_0   RSE
#>           CL     0.15       5     5
#>            V        8       3     3
#>           KA        1      14    14
#>         d_CL     0.07      30    28
#>          d_V     0.02      37    34
#>         d_KA      0.6      27    26
#>     sig_prop     0.01      32    23
#>      sig_add     0.25      26    30
#> 
#> Total running time: 0.344 seconds


summary(out_1)
#> ===============================================================================
#> FINAL RESULTS
#> 
#> Optimized Covariates:
#> Group 1: 100
#> 
#> OFV = 56.032
#> 
#> Efficiency: 
#>   ((exp(ofv_final) / exp(ofv_init))^(1/n_parameters)) = 1.0827
#> 
#> Expected relative standard error
#> (%RSE, rounded to nearest integer):
#>    Parameter   Values   RSE_0   RSE
#>           CL     0.15       5     5
#>            V        8       3     3
#>           KA        1      14    14
#>         d_CL     0.07      30    28
#>          d_V     0.02      37    34
#>         d_KA      0.6      27    26
#>     sig_prop     0.01      32    23
#>      sig_add     0.25      26    30
#> 
#> Total running time: 0.344 seconds
```
