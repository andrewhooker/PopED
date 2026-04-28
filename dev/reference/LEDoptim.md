# Optimization function for D-family, E-family and Laplace approximated ED designs

Optimize the objective function for D-family, E-family and Laplace
approximated ED designs. Right now there is only one optimization
algorithm used in this function

1.  Adaptive random search. See
    [`RS_opt`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md).

This function takes information from the PopED database supplied as an
argument. The PopED database supplies information about the the model,
parameters, design and methods to use. Some of the arguments coming from
the PopED database can be overwritten; if they are supplied then they
are used instead of the arguments from the PopED database.

## Usage

``` r
LEDoptim(
  poped.db,
  model_switch = NULL,
  ni = NULL,
  xt = NULL,
  x = NULL,
  a = NULL,
  bpopdescr = NULL,
  ddescr = NULL,
  maxxt = NULL,
  minxt = NULL,
  maxa = NULL,
  mina = NULL,
  ofv_init = 0,
  fim_init = 0,
  trflag = TRUE,
  header_flag = TRUE,
  footer_flag = TRUE,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  out_file = NULL,
  d_switch = FALSE,
  use_laplace = T,
  laplace.fim = FALSE,
  use_RS = poped.db$settings$bUseRandomSearch,
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

- ofv_init:

  The initial OFV. If set to zero then it is computed.

- fim_init:

  The initial value of the FIM. If set to zero then it is computed.

- trflag:

  Should the optimization be output to the screen and to a file?

- header_flag:

  Should the header text be printed out?

- footer_flag:

  Should the footer text be printed out?

- opt_xt:

  Should the sample times be optimized?

- opt_a:

  Should the continuous design variables be optimized?

- opt_x:

  Should the discrete design variables be optimized?

- out_file:

  Which file should the output be directed to? A string, a file handle
  using [`file`](https://rdrr.io/r/base/connections.html) or `""` will
  output to the screen.

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

- use_RS:

  should the function use a random search algorithm?

- ...:

  arguments passed to
  [`evaluate.fim`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md)
  and
  [`ofv_fim`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md).

## See also

Other Optimize:
[`Doptim()`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md),
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

# warfarin ed model

if (FALSE) { # \dontrun{
  
  LEDoptim(poped.db) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE) 

  LEDoptim(poped.db,opt_xt=T,rsit=10,laplace.fim=TRUE) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,use_laplace=FALSE) 
  
  ## testing header and footer
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           out_file="foobar.txt") 
  
  ff <- LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
                 trflag=FALSE) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           header_flag=FALSE) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           out_file="") 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           footer_flag=FALSE) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           footer_flag=FALSE, header_flag=FALSE) 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           footer_flag=FALSE, header_flag=FALSE,out_file="foobar.txt") 
  
  LEDoptim(poped.db,opt_xt=T,rsit=10,d_switch=TRUE,
           footer_flag=FALSE, header_flag=FALSE,out_file="") 

} # }
```
