# Result function for optimization routines

Create some output to the screen and a text file that summarizes the
problem you solved.

## Usage

``` r
blockfinal(
  fn,
  fmf,
  dmf,
  groupsize,
  ni,
  xt,
  x,
  a,
  model_switch,
  bpop,
  d,
  docc,
  sigma,
  poped.db,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  opt_inds = poped.db$settings$optsw[5],
  fmf_init = NULL,
  dmf_init = NULL,
  param_cvs_init = NULL,
  compute_inv = TRUE,
  out_file = NULL,
  trflag = TRUE,
  footer_flag = TRUE,
  run_time = NULL,
  ...
)
```

## Arguments

- fn:

  The file handle to write to.

- fmf:

  The initial value of the FIM. If set to zero then it is computed.

- dmf:

  The initial OFV. If set to zero then it is computed.

- groupsize:

  A vector of the number of individuals in each group.

- ni:

  A vector of the number of samples in each group.

- xt:

  A matrix of sample times. Each row is a vector of sample times for a
  group.

- x:

  A matrix for the discrete design variables. Each row is a group.

- a:

  A matrix of covariates. Each row is a group.

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

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

- docc:

  Matrix defining the IOV, the IOV variances and the IOV distribution as
  for d and bpop.

- sigma:

  Matrix defining the variances can covariances of the residual
  variability terms of the model. can also just supply the diagonal
  parameter values (variances) as a
  [`c()`](https://rdrr.io/r/base/c.html).

- poped.db:

  A PopED database.

- opt_xt:

  Should the sample times be optimized?

- opt_a:

  Should the continuous design variables be optimized?

- opt_x:

  Should the discrete design variables be optimized?

- opt_inds:

  Are the number of individuals per group being optimized?

- fmf_init:

  Initial FIM.

- dmf_init:

  Initial OFV.

- param_cvs_init:

  The initial design parameter RSE values in percent.

- compute_inv:

  should the inverse of the FIM be used to compute expected RSE values?
  Often not needed except for diagnostic purposes.

- out_file:

  Which file should the output be directed to? A string, a file handle
  using [`file`](https://rdrr.io/r/base/connections.html) or `""` will
  output to the screen.

- trflag:

  Should the optimization be output to the screen and to a file?

- footer_flag:

  Should the footer text be printed out?

- ...:

  arguments passed to
  [`evaluate.fim`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md)
  and
  [`ofv_fim`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md).

## See also

Other Helper:
[`blockexp()`](https://andrewhooker.github.io/PopED/dev/reference/blockexp.md),
[`blockheader()`](https://andrewhooker.github.io/PopED/dev/reference/blockheader.md),
[`blockopt()`](https://andrewhooker.github.io/PopED/dev/reference/blockopt.md)

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


FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)


blockfinal(fn="",fmf=FIM,
           dmf=dmf,
           groupsize=poped.db$design$groupsize,
           ni=poped.db$design$ni,
           xt=poped.db$design$xt,
           x=poped.db$design$x,a=poped.db$design$a,
           model_switch=poped.db$design$model_switch,
           poped.db$parameters$param.pt.val$bpop,
           poped.db$parameters$param.pt.val$d,
           poped.db$parameters$docc,
           poped.db$parameters$param.pt.val$sigma,
           poped.db,
           opt_xt=TRUE,
           fmf_init=FIM,
           dmf_init=dmf,
           param_cvs_init=get_rse(FIM,poped.db))
#> ===============================================================================
#> FINAL RESULTS
#> Optimized Sampling Schedule
#> Group 1:    0.5      1      2      6     24     36     72    120
#> 
#> OFV = 1.14386e+24
#> 
#> Efficiency: 
#>   ((exp(ofv_final) / exp(ofv_init))^(1/n_parameters)) = NaN
#> 
#> Expected relative standard error
#> (%RSE, rounded to nearest integer):
#>    Parameter   Values   RSE_0   RSE
#>           CL     0.15       5     5
#>            V        8       3     3
#>           KA        1      14    14
#>         d_CL     0.07      30    30
#>          d_V     0.02      37    37
#>         d_KA      0.6      27    27
#>     sig_prop     0.01      32    32
#>      sig_add     0.25      26    26
#> 
#> Total running time: 1.599 seconds

```
