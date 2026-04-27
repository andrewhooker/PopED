# Optimize using line search

The function performs a grid search sequentially along design variables.
The grid is defined by ls_step_size.

## Usage

``` r
a_line_search(
  poped.db,
  out_file = "",
  bED = FALSE,
  diff = 0,
  fmf_initial = 0,
  dmf_initial = 0,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5],
  ls_step_size = poped.db$settings$ls_step_size
)
```

## Arguments

- poped.db:

  A PopED database.

- out_file:

  The output file to write to.

- bED:

  If the algorithm should use E-family methods. Logical.

- diff:

  The OFV difference that is deemed significant for changing a design.
  If, by changing a design variable the difference between the new and
  old OFV is less than `diff` the change is not made.

- fmf_initial:

  The initial value of the FIM. If `0` then the FIM is calculated from
  poped.db.

- dmf_initial:

  The initial value of the objective function value (OFV). If `0` then
  the OFV is calculated from poped.db.

- opt_xt:

  Should the sample times be optimized?

- opt_a:

  Should the continuous design variables be optimized?

- opt_x:

  Should the discrete design variables be optimized?

- opt_samps:

  Are the number of sample times per group being optimized?

- opt_inds:

  Are the number of individuals per group being optimized?

- ls_step_size:

  Number of grid points in the line search.

## Value

A list containing:

- fmf:

  The FIM.

- dmf:

  The final value of the objective function value.

- best_changed:

  If the algorithm has found a better design than the starting design.

- xt:

  A matrix of sample times. Each row is a vector of sample times for a
  group.

- x:

  A matrix for the discrete design variables. Each row is a group.

- a:

  A matrix of covariates. Each row is a group.

- poped.db:

  A PopED database.

## See also

Other Optimize:
[`Doptim()`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md),
[`LEDoptim()`](https://andrewhooker.github.io/PopED/dev/reference/LEDoptim.md),
[`RS_opt()`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md),
[`bfgsb_min()`](https://andrewhooker.github.io/PopED/dev/reference/bfgsb_min.md),
[`calc_autofocus()`](https://andrewhooker.github.io/PopED/dev/reference/calc_autofocus.md),
[`calc_ofv_and_grad()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_grad.md),
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


# very sparse grid to evaluate (4 points for each design valiable)
output <- a_line_search(poped.db, opt_xt=TRUE, opt_a=TRUE, ls_step_size=4)
#> *****************************
#>             Line Search
#> 
#> Searching xt1 on group 1
#> group 1 -- xt[1] changed from  0.5 to  0.01
#>      OFV(MF) changed from 55.3964 to 55.7394 
#> group 1 -- xt[1] changed from  0.01 to  90.0025
#>      OFV(MF) changed from 55.7394 to 55.7436 
#> group 1 -- xt[1] changed from  90.0025 to  120
#>      OFV(MF) changed from 55.7436 to 55.8023 
#> Searching xt6 on group 1
#> group 1 -- xt[6] changed from  36 to  0.01
#>      OFV(MF) changed from 55.8023 to 55.8666 
#> group 1 -- xt[6] changed from  0.01 to  90.0025
#>      OFV(MF) changed from 55.8666 to 55.9043 
#> group 1 -- xt[6] changed from  90.0025 to  120
#>      OFV(MF) changed from 55.9043 to 55.9321 
#> Searching xt2 on group 1
#> Searching xt3 on group 1
#> Searching xt7 on group 1
#> group 1 -- xt[7] changed from  72 to  30.0075
#>      OFV(MF) changed from 55.9321 to 55.9484 
#> group 1 -- xt[7] changed from  30.0075 to  90.0025
#>      OFV(MF) changed from 55.9484 to 55.9661 
#> Searching xt8 on group 1
#> Searching xt5 on group 1
#> Searching xt4 on group 1
#>     OFV(MF): 55.9661
#> 
#> Best value for OFV(MF) = 55.9661
#> 
#> Best value for xt:
#> Group 1:      1      2      6     24     90    120    120    120
#> 
#> Searching a1 on individual/group 1
#> group 1 -- a[1] changed from  70 to  75.0025
#>      OFV(MF) changed from 55.9661 to 56.1502 
#> group 1 -- a[1] changed from  75.0025 to  100
#>      OFV(MF) changed from 56.1502 to 56.7067 
#>     OFV(MF): 56.7067
#> Best value for OFV(MF) = 56.7067
#> 
#> Best value for a: 
#> Group 1: 100 [0.01,100]
#> 
#> 
#> Line search run time: 0.302 seconds
#> ***************************
#> 

if (FALSE) { # \dontrun{  
  
  # longer run time
  output <- a_line_search(poped.db,opt_xt=TRUE)
  
  # output to a text file
  output <- a_line_search(poped.db,opt_xt=TRUE,out_file="tmp.txt")
  
} # }
```
