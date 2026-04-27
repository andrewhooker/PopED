# Modified Fedorov Exchange Algorithm

Optimize the objective function using a modified Fedorov exchange
algorithm. The function works for continuous and discrete optimization
variables. This function takes information from the PopED database
supplied as an argument. The PopED database supplies information about
the the model, parameters, design and methods to use. Some of the
arguments coming from the PopED database can be overwritten; if they are
supplied then they are used instead of the arguments from the PopED
database.

## Usage

``` r
mfea(
  poped.db,
  model_switch,
  ni,
  xt,
  x,
  a,
  bpopdescr,
  ddescr,
  maxxt,
  minxt,
  maxa,
  mina,
  fmf,
  dmf,
  EAStepSize = poped.db$settings$EAStepSize,
  ourzero = poped.db$settings$ourzero,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  trflag = T,
  ...
)
```

## Arguments

- poped.db:

  A PopED database.

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

- ni:

  A vector of the number of samples in each group.

- xt:

  A matrix of sample times. Each row is a vector of sample times for a
  group.

- x:

  A matrix for the discrete design variables. Each row is a group.

- a:

  A matrix of covariates. Each row is a group.

- bpopdescr:

  Matrix defining the fixed effects, per row (row number =
  parameter_number) we should have:

  - column 1 the type of the distribution for E-family designs (0 =
    Fixed, 1 = Normal, 2 = Uniform, 3 = User Defined Distribution, 4 =
    lognormal and 5 = truncated normal)

  - column 2 defines the mean.

  - column 3 defines the variance of the distribution (or length of
    uniform distribution).

- ddescr:

  Matrix defining the diagonals of the IIV (same logic as for the
  `bpopdescr`).

- maxxt:

  Matrix or single value defining the maximum value for each xt sample.
  If a single value is supplied then all xt values are given the same
  maximum value.

- minxt:

  Matrix or single value defining the minimum value for each xt sample.
  If a single value is supplied then all xt values are given the same
  minimum value

- maxa:

  Vector defining the max value for each covariate. If a single value is
  supplied then all a values are given the same max value

- mina:

  Vector defining the min value for each covariate. If a single value is
  supplied then all a values are given the same max value

- fmf:

  The initial value of the FIM. If set to zero then it is computed.

- dmf:

  The initial OFV. If set to zero then it is computed.

- EAStepSize:

  Exchange Algorithm StepSize

- ourzero:

  Value to interpret as zero in design

- opt_xt:

  Should the sample times be optimized?

- opt_a:

  Should the continuous design variables be optimized?

- opt_x:

  Should the discrete design variables be optimized?

- trflag:

  Should the optimization be output to the screen and to a file?

- ...:

  arguments passed to
  [`evaluate.fim`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md)
  and
  [`ofv_fim`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md).

## References

1.  J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and
    A.C. Hooker, "PopED: An extended, parallelized, nonlinear mixed
    effects models optimal design tool", Computer Methods and Programs
    in Biomedicine, 108, 2012.

## See also

Other Optimize:
[`Doptim()`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md),
[`LEDoptim()`](https://andrewhooker.github.io/PopED/dev/reference/LEDoptim.md),
[`RS_opt()`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md),
[`a_line_search()`](https://andrewhooker.github.io/PopED/dev/reference/a_line_search.md),
[`bfgsb_min()`](https://andrewhooker.github.io/PopED/dev/reference/bfgsb_min.md),
[`calc_autofocus()`](https://andrewhooker.github.io/PopED/dev/reference/calc_autofocus.md),
[`calc_ofv_and_grad()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_grad.md),
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

##############
# typically one will use poped_optimize 
# This then calls mfea 
##############

# optimization of covariate, with coarse grid
out_1 <- poped_optimize(poped.db,opt_a=1,
                              bUseExchangeAlgorithm=1,
                              EAStepSize=25,out_file = "")
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
#> MFEA - It. : 1
#> MFEA - It. : 1
#> Exchanged covariate 1 in group/ind 1 from 70 to 100
#> Exchanged covariate 1 in group/ind 1 from 70 to 100
#> Delta : 0.0114735   OFV. : 56.032
#> Delta : 0.0114735   OFV. : 56.032
#> MFEA - It. : 2
#> MFEA - It. : 2
#> Delta : 0   OFV. : 56.032
#> Delta : 0   OFV. : 56.032
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
#>         d_CL     0.07       0     0
#>          d_V     0.02      37    34
#>         d_KA      0.6       0     0
#>     sig_prop     0.01      32    23
#>      sig_add     0.25      26    30
#> 
#> Total running time: 0.028 seconds


if (FALSE) { # \dontrun{
  
  
  
  # MFEA optimization with only integer times allowed
  out_2 <- poped_optimize(poped.db,opt_xt=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(out_2$fmf,out_2$poped.db)
  plot_model_prediction(out_2$poped.db)
  
  
  ##############
  # If you really want to you can use mfea dirtectly
  ##############
  dsl <- downsizing_general_design(poped.db)
  
  output <- mfea(poped.db,
                 model_switch=dsl$model_switch,
                 ni=dsl$ni,
                 xt=dsl$xt,
                 x=dsl$x,
                 a=dsl$a,
                 bpopdescr=dsl$bpop,
                 ddescr=dsl$d,
                 maxxt=dsl$maxxt,
                 minxt=dsl$minxt,
                 maxa=dsl$maxa,
                 mina=dsl$mina,
                 fmf=0,dmf=0,
                 EAStepSize=1,
                 opt_xt=1)
  
  
} # }

```
