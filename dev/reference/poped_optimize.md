# Retired optimization module for PopED

This function is an older version of
[`poped_optim`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim.md).
Please use
[`poped_optim`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim.md)
unless you have a specific reason to use this function instead.

## Usage

``` r
poped_optimize(
  poped.db,
  ni = NULL,
  xt = NULL,
  model_switch = NULL,
  x = NULL,
  a = NULL,
  bpop = NULL,
  d = NULL,
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
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5],
  cfaxt = poped.db$settings$cfaxt,
  cfaa = poped.db$settings$cfaa,
  rsit = poped.db$settings$rsit,
  rsit_output = poped.db$settings$rsit_output,
  fim.calc.type = poped.db$settings$iFIMCalculationType,
  ofv_calc_type = poped.db$settings$ofv_calc_type,
  approx_type = poped.db$settings$iApproximationMethod,
  bUseExchangeAlgorithm = poped.db$settings$bUseExchangeAlgorithm,
  iter = 1,
  d_switch = poped.db$settings$d_switch,
  ED_samp_size = poped.db$settings$ED_samp_size,
  bLHS = poped.db$settings$bLHS,
  use_laplace = poped.db$settings$iEDCalculationType,
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

- opt_samps:

  Are the number of sample times per group being optimized?

- opt_inds:

  Are the number of individuals per group being optimized?

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

- ofv_calc_type:

  OFV calculation type for FIM

  - 1 = "D-optimality". Determinant of the FIM: det(FIM)

  - 2 = "A-optimality". Inverse of the sum of the expected parameter
    variances: 1/trace_matrix(inv(FIM))

  - 4 = "lnD-optimality". Natural logarithm of the determinant of the
    FIM: log(det(FIM))

  - 6 = "Ds-optimality". Ratio of the Determinant of the FIM and the
    Determinant of the uninteresting rows and columns of the FIM:
    det(FIM)/det(FIM_u)

  - 7 = Inverse of the sum of the expected parameter RSE:
    1/sum(get_rse(FIM,poped.db,use_percent=FALSE))

- approx_type:

  Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI.

- bUseExchangeAlgorithm:

  Use Exchange algorithm (1=TRUE, 0=FALSE)

- iter:

  The number of iterations entered into the `blockheader_2` function.

- d_switch:

  - **\*\*\*\*\*\*START OF CRITERION SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  D-family design (1) or ED-family design (0) (with or without parameter
  uncertainty)

- ED_samp_size:

  Sample size for E-family sampling

- bLHS:

  How to sample from distributions in E-family calculations. 0=Random
  Sampling, 1=LatinHyperCube –

- use_laplace:

  Should the Laplace method be used in calculating the expectation of
  the OFV?

- ...:

  arguments passed to other functions. See
  [`Doptim`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md).

## Details

This function optimized the objective function. The function works for
both discrete and continuous optimization variables. This function takes
information from the PopED database supplied as an argument. The PopED
database supplies information about the the model, parameters, design
and methods to use. Some of the arguments coming from the PopED database
can be overwritten; if they are supplied then they are used instead of
the arguments from the PopED database.

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
[`poped_optim_3()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_3.md)

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

# below are a number of ways to optimize the problem

# RS+SG+LS optimization of DOSE and sample times
# optimization with just a few iterations
# only to check that things are working
out_1 <- poped_optimize(poped.db,opt_a=TRUE,opt_xt=TRUE,
                         rsit=2,sgit=2,ls_step_size=2, 
                         iter_max=1,out_file = "")
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
#> * Optimize Sampling Schedule
#> * Optimize Covariates
#> 
#> *******************************
#> Initial Value
#>  OFV(mf) = 55.3964
#> *******************************
#> 
#> RS - It. : 2   OFV : 55.3964
#> 
#> *******************************
#> RS Results
#>  OFV(mf) = 55.3964
#> 
#> Optimized Sampling Schedule
#> Group 1:    0.5      1      2      6     24     36     72    120
#> 
#> Optimized Covariates:
#> Group 1: 70
#> 
#> *********************************
#> 
#> Run time for random search: 0.014 seconds
#> 
#> SG - It. : 1  OFV : 55.45   Diff. :     1
#> SG - It. : 2  OFV : 55.61   Diff. : 0.002799
#> 
#> SG - Iteration 2 --------- FINAL -------------------------
#> Normalized gradient: Grad_xt(OFV)/OFV
#> -0.00262405729196299
#>  -0.000132025609212243
#>  0.0049562882696737
#>  -0.000109893316026729
#>  -0.000341253164434824
#>  2.72645074812749e-05
#>  0.000160461661950708
#>  3.06020481312946e-05
#> xt opt:
#> 0.210318514650852
#>  0.950298514650852
#>  2.28968148534915
#>  5.71031851465085
#>  23.7103185146509
#>  36.2896814853491
#>  72.2896814853491
#>  120
#> Normalized gradient: Grad_a(OFV)/OFV
#> 6.312896e-04
#> aopt:
#> 7.024140e+01
#> OFV(mf)    : 55.6095
#> diff       : 0.00279928
#> *************************************************************
#> Stochastic gradient run time: 0.358 seconds
#> 
#> *****************************
#>             Line Search
#> 
#> Searching xt3 on group 1
#> Searching xt8 on group 1
#> Searching xt1 on group 1
#> group 1 -- xt[1] changed from  0.210319 to  0.01
#>      OFV(MF) changed from 55.6095 to 55.8352 
#> group 1 -- xt[1] changed from  0.01 to  120
#>      OFV(MF) changed from 55.8352 to 55.8915 
#> Searching xt2 on group 1
#> Searching xt5 on group 1
#> Searching xt6 on group 1
#> group 1 -- xt[6] changed from  36.2897 to  0.01
#>      OFV(MF) changed from 55.8915 to 55.9642 
#> group 1 -- xt[6] changed from  0.01 to  120
#>      OFV(MF) changed from 55.9642 to 56.0271 
#> Searching xt7 on group 1
#> group 1 -- xt[7] changed from  72.2897 to  120
#>      OFV(MF) changed from 56.0271 to 56.0561 
#> Searching xt4 on group 1
#>     OFV(MF): 56.0561
#> 
#> Best value for OFV(MF) = 56.0561
#> 
#> Best value for xt:
#> Group 1: 0.9503   2.29   5.71  23.71    120    120    120    120
#> 
#> Searching a1 on individual/group 1
#> group 1 -- a[1] changed from  70.2414 to  100
#>      OFV(MF) changed from 56.0561 to 56.8149 
#>     OFV(MF): 56.8149
#> Best value for OFV(MF) = 56.8149
#> 
#> Best value for a: 
#> Group 1: 100 [0.01,100]
#> 
#> 
#> Line search run time: 0.202 seconds
#> ***************************
#> 
#> ===============================================================================
#> FINAL RESULTS
#> Optimized Sampling Schedule
#> Group 1: 0.9503   2.29   5.71  23.71    120    120    120    120
#> 
#> Optimized Covariates:
#> Group 1: 100
#> 
#> OFV = 56.8149
#> 
#> Efficiency: 
#>   ((exp(ofv_final) / exp(ofv_init))^(1/n_parameters)) = 1.194
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
#>     sig_prop     0.01      32    24
#>      sig_add     0.25      26    16
#> 
#> Total running time: 0.577 seconds

if (FALSE) { # \dontrun{
  
  # RS+SG+LS optimization of sample times 
  # (longer run time than above but more likely to reach a maximum)
  output <- poped_optimize(poped.db,opt_xt=T)
  get_rse(output$fmf,output$poped.db)
  plot_model_prediction(output$poped.db)
  
  # MFEA optimization with only integer times allowed
  mfea.output <- poped_optimize(poped.db,opt_xt=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(mfea.output$fmf,mfea.output$poped.db)
  plot_model_prediction(mfea.output$poped.db)
  
  # Examine efficiency of sampling windows
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5)
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1)
  
  # Random search (just a few samples here)
  rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                              bUseRandomSearch= 1,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0)
  get_rse(rs.output$fmf,rs.output$poped.db)
  
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
  # E-family Optimization
  ##############
  
  # Adding 10% log-normal Uncertainty to fixed effects (not Favail)
  bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
  bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                           bpop_vals,
                           ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
  bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
  bpop_vals_ed_ln
  
  ## -- Define initial design  and design space
  poped.db <- create.poped.database(
    ff_fun=ff.PK.1.comp.oral.sd.CL,
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
  
  # ED optimization using Random search (just a few samples here)
  output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=10,d_switch=0)
  get_rse(output$fmf,output$poped.db)
  
  # ED with laplace approximation, 
  # optimization using Random search (just a few samples here)
  output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=10,
                           d_switch=0,use_laplace=TRUE,laplace.fim=TRUE)
  get_rse(output$fmf,output$poped.db)
  
  
} # }
```
