# Calculate the Fisher Information Matrix (FIM) and the OFV(FIM) for either point values or parameters or distributions.

This function computes the expectation of the FIM and OFV(FIM) for
either point values of parameter estimates or parameter distributions
given the model, parameters, distributions of parameter uncertainty,
design and methods defined in the PopED database.

## Usage

``` r
calc_ofv_and_fim(
  poped.db,
  ofv = 0,
  fim = 0,
  d_switch = poped.db$settings$d_switch,
  bpopdescr = poped.db$parameters$bpop,
  ddescr = poped.db$parameters$d,
  bpop = bpopdescr[, 2, drop = F],
  d = getfulld(ddescr[, 2, drop = F], poped.db$parameters$covd),
  docc_full = getfulld(poped.db$parameters$docc[, 2, drop = F],
    poped.db$parameters$covdocc),
  model_switch = poped.db$design$model_switch,
  ni = poped.db$design$ni,
  xt = poped.db$design$xt,
  x = poped.db$design$x,
  a = poped.db$design$a,
  fim.calc.type = poped.db$settings$iFIMCalculationType,
  use_laplace = poped.db$settings$iEDCalculationType,
  laplace.fim = FALSE,
  ofv_fun = poped.db$settings$ofv_fun,
  evaluate_fim = TRUE,
  ...
)
```

## Arguments

- poped.db:

  A PopED database.

- ofv:

  The current ofv. If other than zero then this value is simply returned
  unchanged.

- fim:

  The current FIM. If other than zero then this value is simply returned
  unchanged.

- d_switch:

  - **\*\*\*\*\*\*START OF CRITERION SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  D-family design (1) or ED-family design (0) (with or without parameter
  uncertainty)

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

- docc_full:

  A between occasion variability matrix.

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

- fim.calc.type:

  The method used for calculating the FIM. Potential values:

  - 0 = Full FIM. No assumption that fixed and random effects are
    uncorrelated.

  - 1 = Reduced FIM. Assume that there is no correlation in the FIM
    between the fixed and random effects, and set these elements in the
    FIM to zero.

  - 2 = weighted models (placeholder).

  - 3 = Not currently used.

  - 4 = Reduced FIM and computing all derivatives with respect to the
    standard deviation of the residual unexplained variation
    (sqrt(SIGMA) in NONMEM). This matches what is done in PFIM, and
    assumes that the standard deviation of the residual unexplained
    variation is the estimated parameter (NOTE: NONMEM estimates the
    variance of the residual unexplained variation by default).

  - 5 = Full FIM parameterized with A,B,C matrices & derivative of
    variance.

  - 6 = Calculate one model switch at a time, good for large matrices.

  - 7 = Reduced FIM parameterized with A,B,C matrices & derivative of
    variance.

- use_laplace:

  Should the Laplace method be used in calculating the expectation of
  the OFV?

- laplace.fim:

  Should an E(FIM) be calculated when computing the Laplace approximated
  E(OFV). Typically the FIM does not need to be computed and, if
  desired, this calculation is done using the standard MC integration
  technique, so can be slow.

- ofv_fun:

  User defined function used to compute the objective function. The
  function must have a poped database object as its first argument and
  have "..." in its argument list. Can be referenced as a function or as
  a file name where the function defined in the file has the same name
  as the file. e.g. "cost.txt" has a function named "cost" in it.

- evaluate_fim:

  Should the FIM be calculated?

- ...:

  Other arguments passed to the function.

## Value

A list containing the FIM and OFV(FIM) or the E(FIM) and E(OFV(FIM))
according to the function arguments.

## See also

Other FIM:
[`LinMatrixH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixH.md),
[`LinMatrixLH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixLH.md),
[`LinMatrixL_occ()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixL_occ.md),
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

Other E-family:
[`ed_laplace_ofv()`](https://andrewhooker.github.io/PopED/dev/reference/ed_laplace_ofv.md),
[`ed_mftot()`](https://andrewhooker.github.io/PopED/dev/reference/ed_mftot.md),
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md)

Other evaluate_FIM:
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md),
[`evaluate.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)

## Examples

``` r
library(PopED)

############# START #################
## Create PopED database
## (warfarin model for optimization
##  with parameter uncertainty)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error
## to avoid sample times at very low concentrations (time 0 or very late samoples).

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
#> <bytecode: 0x5641610772b0>
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

# Adding 10% log-normal Uncertainty to fixed effects (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                         bpop_vals,
                         ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
bpop_vals_ed_ln
#>          bpop_vals         
#> CL     4      0.15 0.000225
#> V      4      8.00 0.640000
#> KA     4      1.00 0.010000
#> Favail 0      1.00 0.000000

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.add.prop,
                                  bpop=bpop_vals_ed_ln, 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70,
                                  mina=0,
                                  maxa=100)

############# END ###################
## Create PopED database
## (warfarin model for optimization
##  with parameter uncertainty)
#####################################


calc_ofv_and_fim(poped.db)
#> $ofv
#> [1] 55.39645
#> 
#> $fim
#>             [,1]      [,2]      [,3]         [,4]         [,5]        [,6]
#> [1,] 17141.83891 20.838375 10.011000 0.000000e+00     0.000000  0.00000000
#> [2,]    20.83837 17.268051 -3.423641 0.000000e+00     0.000000  0.00000000
#> [3,]    10.01100 -3.423641 49.864697 0.000000e+00     0.000000  0.00000000
#> [4,]     0.00000  0.000000  0.000000 2.324341e+03     9.770352  0.03523364
#> [5,]     0.00000  0.000000  0.000000 9.770352e+00 19083.877564 11.72131703
#> [6,]     0.00000  0.000000  0.000000 3.523364e-02    11.721317 38.85137516
#> [7,]     0.00000  0.000000  0.000000 7.268410e+02  9656.158553 64.78095548
#> [8,]     0.00000  0.000000  0.000000 9.062739e+01   266.487127  2.94728469
#>              [,7]        [,8]
#> [1,]      0.00000    0.000000
#> [2,]      0.00000    0.000000
#> [3,]      0.00000    0.000000
#> [4,]    726.84097   90.627386
#> [5,]   9656.15855  266.487127
#> [6,]     64.78096    2.947285
#> [7,] 192840.20092 6659.569867
#> [8,]   6659.56987  475.500111
#> 

if (FALSE) { # \dontrun{
  
  calc_ofv_and_fim(poped.db,d_switch=0)
  calc_ofv_and_fim(poped.db,d_switch=0,use_laplace=TRUE)
  calc_ofv_and_fim(poped.db,d_switch=0,use_laplace=TRUE,laplace.fim=TRUE)

} # }
```
