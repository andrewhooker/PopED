# Header function for optimization routines

Create some output to the screen and a text file that summarizes the
problem you are tying to solve.

## Usage

``` r
blockheader(
  poped.db,
  name = "Default",
  iter = NULL,
  e_flag = !(poped.db$settings$d_switch),
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5],
  fmf = 0,
  dmf = 0,
  bpop = NULL,
  d = NULL,
  docc = NULL,
  sigma = NULL,
  name_header = poped.db$settings$strOutputFileName,
  file_path = poped.db$settings$strOutputFilePath,
  out_file = NULL,
  compute_inv = TRUE,
  trflag = TRUE,
  header_flag = TRUE,
  ...
)
```

## Arguments

- poped.db:

  A PopED database.

- name:

  The name used for the output file. Combined with `name_header` and
  `iter`. If `""` then output is to the screen.

- iter:

  The last number in the name printed to the output file, combined with
  `name`.

- e_flag:

  Should output be with uncertainty around parameters?

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

- fmf:

  The initial value of the FIM. If set to zero then it is computed.

- dmf:

  The initial OFV. If set to zero then it is computed.

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

- name_header:

  The initial portion of the file name.

- file_path:

  The path to where the file should be created.

- out_file:

  Which file should the output be directed to? A string, a file handle
  using [`file`](https://rdrr.io/r/base/connections.html) or `""` will
  output to the screen.

- compute_inv:

  should the inverse of the FIM be used to compute expected RSE values?
  Often not needed except for diagnostic purposes.

- trflag:

  Should the optimization be output to the screen and to a file?

- header_flag:

  Should the header text be printed out?

- ...:

  Additional arguments passed to further functions.

## Value

fn A file handle (or `''` if `name=''`)

## See also

Other Helper:
[`blockexp()`](https://andrewhooker.github.io/PopED/dev/reference/blockexp.md),
[`blockfinal()`](https://andrewhooker.github.io/PopED/dev/reference/blockfinal.md),
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

blockheader(poped.db,name="")
#> ==============================================================================
#> Optimization of design parameters
#> 
#> 
#> [1] ""

blockheader(name="",iter=1,poped.db)
#> ==============================================================================
#> Optimization of design parameters
#> 
#> 
#> [1] ""


blockheader(name='',
              iter=1,
              poped.db,
              e_flag=FALSE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$settings$optsw[4],
              opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$parameters$param.pt.val$bpop,
              d=poped.db$parameters$param.pt.val$d,
              docc=poped.db$parameters$docc,sigma=poped.db$parameters$param.pt.val$sigma)
#> ===============================================================================
#> Initial design evaluation
#> 
#> Initial OFV = 1.14386e+24
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
#> [1] ""



blockheader(name='',
              iter=1,
              poped.db,
              e_flag=TRUE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$settings$optsw[4],
              opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$parameters$param.pt.val$bpop,
              d=poped.db$parameters$param.pt.val$d,
              docc=poped.db$parameters$docc,sigma=poped.db$parameters$param.pt.val$sigma)
#> ===============================================================================
#> Initial design evaluation
#> 
#> Initial OFV = 1.14386e+24
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
#> [1] ""
  
  
poped.db.1 <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.add.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=rbind(c( 0.5,1,2,6,24,36,72,120),
                                           c( 0.5,1.1,2,6,24,36,72,120)),
                                  minxt=rbind(c(0,1,1.5,3,20,30,70,118),
                                              c(0.1,1.1,1.6,3.1,20.1,30.1,70.1,118.1)),
                                  maxxt=c(12,13,14,15,26,44,78,120),
                                  a=70,
                                  mina=0,
                                  maxa=100)


blockheader(poped.db.1,name="",trflag=2,opt_xt=TRUE)
#> ==============================================================================
#> Model description : PopED model 
#> 
#> Model Sizes : 
#> Number of individual model parameters                  g[j]    : Ng    = 5
#> Number of population model fixed parameters            bpop[j] : Nbpop = 4
#> Number of population model random effects parameters   b[j]    : Nb    = 3
#> 
#> Typical Population Parameters:
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> bpop[1]:  0.15 
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> bpop[2]:     8 
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> bpop[3]:     1 
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> bpop[4]:     1 
#> 
#> Between Subject Variability matrix D (variance units) 
#> 0.07 0.00 0.00
#> 0.00 0.02 0.00
#> 0.00 0.00 0.60
#> 
#> Diagonal Elements of D [sqrt(param)]:
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> D[1,1]:  0.07 [0.2646] 
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> D[2,2]:  0.02 [0.1414] 
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> D[3,3]:   0.6 [0.7746] 
#> 
#> Residual Unexplained Variability matrix SIGMA (variance units) : 
#> 0.01 0.00
#> 0.00 0.25
#> 
#> Diagonal Elements of SIGMA [sqrt(param)]:
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> SIGMA[1,1]:  0.01 [  0.1] 
#> Warning: 5 arguments not used by format '%s[%g%s]: %5.4g %s
#> '
#> SIGMA[2,2]:  0.25 [  0.5] 
#> 
#> ==============================================================================
#> Experiment description (design and design space)
#> 
#> Warning: 2 arguments not used by format 'Number of individuals: %g
#> '
#> Number of individuals: 64
#> Number of groups (individuals with same design): 2
#> Number of individuals per group:
#>  
#> Warning: 2 arguments not used by format '    Group %g: %g
#> '
#>     Group 1: 32
#>      Group 2: 32
#> Number of samples per group:
#>  Number of discrete experimental variables: 0
#> Number of model covariates: 1
#> 
#> Initial Sampling Schedule
#> Group 1:    0.5      1      2      6     24     36     72    120
#> Group 2:    0.5    1.1      2      6     24     36     72    120
#> 
#> Minimum allowed sampling values
#> Group 1:  1e-05      1    1.5      3     20     30     70    118
#> Group 2:    0.1    1.1    1.6    3.1   20.1   30.1   70.1  118.1
#> 
#> Maximum allowed sampling values
#> Group 1:     12     13     14     15     26     44     78    120
#> Group 2:     12     13     14     15     26     44     78    120
#> 
#> Covariates:
#> Group 1: 
#> Warning: 2 arguments not used by format '%g'
#> 70
#> Group 2: 
#> Warning: 2 arguments not used by format '%g'
#> 70
#> 
#> ==============================================================================
#> Criterion Specification
#> 
#> OFV calculation for FIM: 4 
#>   1=Determinant of FIM,
#>   4=log determinant of FIM,
#>   6=determinant of interesting part of FIM (Ds)
#> 
#> Approximation method: 0
#>   0=FO, 
#>   1=FOCE, 
#>   2=FOCEI, 
#>   3=FOI
#> 
#> Fisher Information Matrix type: 1
#>   0=Full FIM,
#>   1=Reduced FIM,
#>   2=weighted models,
#>   3=Loc models,
#>   4=reduced FIM with derivative of SD of sigma as pfim,
#>   5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
#>   6=Calculate one model switch at a time, good for large matrices,
#>   7=Reduced FIM parameterized with A,B,C matrices & derivative of variance
#> 
#> Design family: 1
#>   D-family design (1) or 
#>   ED-family design (0) 
#>   (with or without parameter uncertainty)
#> 
#> ==============================================================================
#> Optimization of design parameters
#> 
#> * Optimize Sampling Schedule
#> 
#> [1] ""

```
