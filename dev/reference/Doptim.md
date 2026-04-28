# D-family optimization function

Optimize the objective function. There are 4 different optimization
algorithms used in this function

1.  Adaptive random search. See
    [`RS_opt`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md).

2.  Stochastic gradient.

3.  A Broyden Fletcher Goldfarb Shanno (BFGS) method for nonlinear
    minimization with box constraints.

4.  A line search. See
    [`a_line_search`](https://andrewhooker.github.io/PopED/dev/reference/a_line_search.md).

The optimization algorithms run in series, taking as input the output
from the previous method. The stopping rule used is to test if the line
search algorithm fids a better optimum then its initial value. If so,
then the chain of algorithms is run again. If line search is not used
then the argument `iter_tot` defines the number of times the chain of
algorithms is run. This function takes information from the PopED
database supplied as an argument. The PopED database supplies
information about the the model, parameters, design and methods to use.
Some of the arguments coming from the PopED database can be overwritten;
if they are supplied then they are used instead of the arguments from
the PopED database.

## Usage

``` r
Doptim(
  poped.db,
  ni,
  xt,
  model_switch,
  x,
  a,
  bpopdescr,
  ddescr,
  maxxt,
  minxt,
  maxa,
  mina,
  fmf = 0,
  dmf = 0,
  trflag = TRUE,
  bUseRandomSearch = poped.db$settings$bUseRandomSearch,
  bUseStochasticGradient = poped.db$settings$bUseStochasticGradient,
  bUseBFGSMinimizer = poped.db$settings$bUseBFGSMinimizer,
  bUseLineSearch = poped.db$settings$bUseLineSearch,
  sgit = poped.db$settings$sgit,
  ls_step_size = poped.db$settings$ls_step_size,
  BFGSConvergenceCriteriaMinStep = poped.db$settings$BFGSConvergenceCriteriaMinStep,
  BFGSProjectedGradientTol = poped.db$settings$BFGSProjectedGradientTol,
  BFGSTolerancef = poped.db$settings$BFGSTolerancef,
  BFGSToleranceg = poped.db$settings$BFGSToleranceg,
  BFGSTolerancex = poped.db$settings$BFGSTolerancex,
  iter_tot = poped.db$settings$iNumSearchIterationsIfNotLineSearch,
  iter_max = 10,
  ...
)
```

## Arguments

- poped.db:

  A PopED database.

- ni:

  A vector of the number of samples in each group.

- xt:

  A matrix of sample times. Each row is a vector of sample times for a
  group.

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

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

- trflag:

  Should the optimization be output to the screen and to a file?

- bUseRandomSearch:

  - **\*\*\*\*\*\*START OF Optimization algorithm SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  Use random search (1=TRUE, 0=FALSE)

- bUseStochasticGradient:

  Use Stochastic Gradient search (1=TRUE, 0=FALSE)

- bUseBFGSMinimizer:

  Use BFGS Minimizer (1=TRUE, 0=FALSE)

- bUseLineSearch:

  Use Line search (1=TRUE, 0=FALSE)

- sgit:

  Number of stochastic gradient iterations

- ls_step_size:

  Number of grid points in the line search.

- BFGSConvergenceCriteriaMinStep:

  BFGS Minimizer Convergence Criteria Minimum Step

- BFGSProjectedGradientTol:

  BFGS Minimizer Convergence Criteria Normalized Projected Gradient
  Tolerance

- BFGSTolerancef:

  BFGS Minimizer Line Search Tolerance f

- BFGSToleranceg:

  BFGS Minimizer Line Search Tolerance g

- BFGSTolerancex:

  BFGS Minimizer Line Search Tolerance x

- iter_tot:

  Number of iterations to use if line search is not used. Must be less
  than `iter_max` to be used.

- iter_max:

  If line search is used then the algorithm tests if line search (always
  run at the end of the optimization iteration) changes the design in
  any way. If not, the algorithm stops. If yes, then a new iteration is
  run unless `iter_max` iterations have already been run.

- ...:

  arguments passed to
  [`evaluate.fim`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md)
  and
  [`ofv_fim`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md).

## References

1.  M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a
    software for optimal experimental design in population kinetics",
    Computer Methods and Programs in Biomedicine, 74, 2004.

2.  J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and
    A.C. Hooker, "PopED: An extended, parallelized, nonlinear mixed
    effects models optimal design tool", Computer Methods and Programs
    in Biomedicine, 108, 2012.

## See also

Other Optimize:
[`LEDoptim()`](https://andrewhooker.github.io/PopED/dev/reference/LEDoptim.md),
[`RS_opt()`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md),
[`a_line_search()`](https://andrewhooker.github.io/PopED/dev/reference/a_line_search.md),
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


if (FALSE) { # \dontrun{
  
  ##############
  # typically one will use poped_optimize 
  # This then calls Doptim for continuous optimization problems
  ##############
  
  
  # RS+SG+LS optimization of sample times
  # optimization with just a few iterations
  # only to check that things are working
  output <- poped_optimize(poped.db,opt_xt=T,
                           rsit=5,sgit=5,ls_step_size=5)
  
  # RS+SG+LS optimization of sample times 
  # (longer run time than above but more likely to reach a maximum)
  output <- poped_optimize(poped.db,opt_xt=T)
  get_rse(output$fmf,output$poped.db)
  plot_model_prediction(output$poped.db)
  
  
  # Random search (just a few samples here)
  rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                              bUseRandomSearch= 1,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0)
  
  # line search, DOSE and sample time optimization
  ls.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 1,
                              ls_step_size=10)
  
  # Stochastic gradient search, DOSE and sample time optimization
  sg.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1, 
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 1,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0,
                              sgit=20)
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)

  ##############
  # If you really want to you can use Doptim dirtectly
  ##############
  dsl <- downsizing_general_design(poped.db)
  poped.db$settings$optsw[2] <- 1  # sample time optimization
  output <- Doptim(poped.db,dsl$ni, dsl$xt, dsl$model_switch, dsl$x, dsl$a, 
         dsl$bpop, dsl$d, dsl$maxxt, dsl$minxt,dsl$maxa,dsl$mina) 
  
} # }

```
