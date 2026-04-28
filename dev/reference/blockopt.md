# Summarize your optimization settings for optimization routines

Create some output to the screen and a text file that summarizes the
optimization settings you will use to optimize.

## Usage

``` r
blockopt(fn, poped.db, opt_method = "")
```

## Arguments

- fn:

  The file handle to write to.

- poped.db:

  A PopED database.

- opt_method:

  If "RS" (random search), "SG" (stochastic gradient) or "DO" (discrete
  optimization) then specific output is produced.

## See also

Other Helper:
[`blockexp()`](https://andrewhooker.github.io/PopED/dev/reference/blockexp.md),
[`blockfinal()`](https://andrewhooker.github.io/PopED/dev/reference/blockfinal.md),
[`blockheader()`](https://andrewhooker.github.io/PopED/dev/reference/blockheader.md)

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


blockopt(fn="",poped.db,opt_method="SG")
#> ==============================================================================
#> Optimization Settings
#> 
#> Stochastic Gradient :
#> Maximum number of cycles : 150
#> Epsilon for termination : 1e-08
#> First step factor for xt: 0.001
#> First step factor for a: 0.001
#> RS m0it: 50
#> 
#> NULL
blockopt(fn="",poped.db,opt_method="RS")
#> ==============================================================================
#> Optimization Settings
#> 
#> Random Search :
#> Number of cycles : 300
#> Locality factor for xt : 10
#> Locality factor for a  : 10
#> 
#> NULL
blockopt(fn="",poped.db,opt_method="DO")
#> ==============================================================================
#> Optimization Settings
#> 
#> Discrete Optimization  :
#> RS int it: 250
#> SG int it: 50
#> 
#> NULL
```
