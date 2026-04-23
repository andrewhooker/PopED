# Optimize the objective function using an adaptive random search algorithm for D-family and E-family designs.

Optimize the objective function using an adaptive random search
algorithm. Optimization can be performed for both D-family and E-family
designs. The function works for both discrete and continuous
optimization variables. This function takes information from the PopED
database supplied as an argument. The PopED database supplies
information about the the model, parameters, design and methods to use.
Some of the arguments coming from the PopED database can be overwritten;
by default these arguments are `NULL` in the function, if they are
supplied then they are used instead of the arguments from the PopED
database.

## Usage

``` r
RS_opt(
  poped.db,
  ni = NULL,
  xt = NULL,
  model_switch = NULL,
  x = NULL,
  a = NULL,
  bpopdescr = NULL,
  ddescr = NULL,
  maxxt = NULL,
  minxt = NULL,
  maxa = NULL,
  mina = NULL,
  fmf = 0,
  dmf = 0,
  trflag = TRUE,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  cfaxt = poped.db$settings$cfaxt,
  cfaa = poped.db$settings$cfaa,
  rsit = poped.db$settings$rsit,
  rsit_output = poped.db$settings$rsit_output,
  fim.calc.type = poped.db$settings$iFIMCalculationType,
  approx_type = poped.db$settings$iApproximationMethod,
  iter = NULL,
  d_switch = poped.db$settings$d_switch,
  use_laplace = poped.db$settings$iEDCalculationType,
  laplace.fim = FALSE,
  header_flag = TRUE,
  footer_flag = TRUE,
  out_file = NULL,
  compute_inv = TRUE,
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

- opt_xt:

  Should the sample times be optimized?

- opt_a:

  Should the continuous design variables be optimized?

- opt_x:

  Should the discrete design variables be optimized?

- cfaxt:

  First step factor for sample times

- cfaa:

  Stochastic Gradient search first step factor for covariates

- rsit:

  Number of Random search iterations

- rsit_output:

  Number of iterations in random search between screen output

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

- approx_type:

  Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI.

- iter:

  The number of iterations entered into the `blockheader_2` function.

- d_switch:

  - **\*\*\*\*\*\*START OF CRITERION SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  D-family design (1) or ED-family design (0) (with or without parameter
  uncertainty)

- use_laplace:

  Should the Laplace method be used in calculating the expectation of
  the OFV?

- laplace.fim:

  Should an E(FIM) be calculated when computing the Laplace approximated
  E(OFV). Typically the FIM does not need to be computed and, if
  desired, this calculation is done using the standard MC integration
  technique, so can be slow.

- header_flag:

  Should the header text be printed out?

- footer_flag:

  Should the footer text be printed out?

- out_file:

  Which file should the output be directed to? A string, a file handle
  using [`file`](https://rdrr.io/r/base/connections.html) or `""` will
  output to the screen.

- compute_inv:

  should the inverse of the FIM be used to compute expected RSE values?
  Often not needed except for diagnostic purposes.

- ...:

  arguments passed to
  [`evaluate.fim`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md)
  and
  [`ofv_fim`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md).

## References

1.  M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a
    software fir optimal experimental design in population kinetics",
    Computer Methods and Programs in Biomedicine, 74, 2004.

2.  J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and
    A.C. Hooker, "PopED: An extended, parallelized, nonlinear mixed
    effects models optimal design tool", Computer Methods and Programs
    in Biomedicine, 108, 2012.

## See also

Other Optimize:
[`Doptim()`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md),
[`LEDoptim()`](https://andrewhooker.github.io/PopED/dev/reference/LEDoptim.md),
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
#> <bytecode: 0x55b94b2139b0>
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


# Just a few iterations, optimize on DOSE and sample times using the full FIM
out_1 <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=3,fim.calc.type=0, out_file = "")
#> ===============================================================================
#> Initial design evaluation
#> 
#> Initial OFV = 57.0828
#> 
#> Initial design
#> expected relative standard error
#> (%RSE, rounded to nearest integer)
#>     Parameter   Values   RSE_0
#>            CL     0.15       5
#>             V        8       3
#>            KA        1       7
#>          d_CL     0.07      30
#>           d_V     0.02      37
#>          d_KA      0.6      28
#>    SIGMA[1,1]     0.01      33
#>    SIGMA[2,2]     0.25      26
#> 
#> ==============================================================================
#> Optimization of design parameters
#> 
#> * Optimize Sampling Schedule
#> * Optimize Covariates
#> 
#> *******************************
#> Initial Value
#>  OFV(mf) = 57.0828
#> *******************************
#> 
#> RS - It. : 3   OFV : 57.0828
#> 
#> *******************************
#> RS Results
#>  OFV(mf) = 57.0828
#> 
#> Optimized Sampling Schedule
#> Group 1:    0.5      1      2      6     24     36     72    120
#> 
#> Optimized Covariates:
#> Group 1: 70
#> 
#> *********************************
#> 
#> ===============================================================================
#> FINAL RESULTS
#> Optimized Sampling Schedule
#> Group 1:    0.5      1      2      6     24     36     72    120
#> 
#> Optimized Covariates:
#> Group 1: 70
#> 
#> OFV = 57.0828
#> 
#> Efficiency: 
#>   ((exp(ofv_final) / exp(ofv_init))^(1/n_parameters)) = 1
#> 
#> Expected relative standard error
#> (%RSE, rounded to nearest integer):
#>     Parameter   Values   RSE_0   RSE
#>            CL     0.15       5     5
#>             V        8       3     3
#>            KA        1       7     7
#>          d_CL     0.07       0     0
#>           d_V     0.02      37    37
#>          d_KA      0.6       0     0
#>    SIGMA[1,1]     0.01      33    33
#>    SIGMA[2,2]     0.25      26    26
#> 
#> Total running time: 0.034 seconds

if (FALSE) { # \dontrun{
  
  RS_opt(poped.db)
  
  RS_opt(poped.db,opt_xt=TRUE,rsit=100,compute_inv=F)
  RS_opt(poped.db,opt_xt=TRUE,rsit=20,d_switch=0)
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,d_switch=0,use_laplace=T)
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,d_switch=0,use_laplace=T,laplace.fim=T)
  
  ## Different headers and footers of output
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,out_file="foo.txt")
  output <- RS_opt(poped.db,opt_xt=TRUE,rsit=100,trflag=FALSE)
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,out_file="")
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,header_flag=FALSE)
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,footer_flag=FALSE)
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,header_flag=FALSE,footer_flag=FALSE)
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,header_flag=FALSE,footer_flag=FALSE,out_file="foo.txt")
  RS_opt(poped.db,opt_xt=TRUE,rsit=10,header_flag=FALSE,footer_flag=FALSE,out_file="") 

} # }
```
