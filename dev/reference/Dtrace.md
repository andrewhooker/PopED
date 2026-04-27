# Trace optimization routines

A helper function for writing output to the screen and files when
optimizing.

## Usage

``` r
Dtrace(
  fn,
  it,
  ni,
  xtopt,
  xopt,
  aopt,
  gxt,
  ga,
  dmf,
  diff,
  ixt,
  ia,
  itvector,
  dmfvector,
  poped.db,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5],
  rsit = poped.db$settings$rsit,
  convergence_eps = poped.db$settings$convergence_eps
)
```

## Arguments

- fn:

  A file to output information to. Can also be the screen if `''`.

- it:

  the iteration number.

- ni:

  A vector of the number of samples in each group.

- xtopt:

  The matrix defining current best sampling schedule.

- xopt:

  The cell structure defining the current best discrete design
  variables.

- aopt:

  The matrix defining the current best continuous design variables.

- gxt:

  The matrix defining the current gradient of the xt vector.

- ga:

  The matrix defining the current gradient for the continuous design
  variables.

- dmf:

  The current OFV.

- diff:

  The difference from the previous iteration.

- ixt:

  If xt Gradient Inversion occurred or not.

- ia:

  If a Gradient Inversion occurred or not.

- itvector:

  The iteration vector. Not currently used.

- dmfvector:

  The dmf vector. Not currently used.

- poped.db:

  A PopED database.

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

- rsit:

  Number of Random search iterations

- convergence_eps:

  Stochastic Gradient convergence value, (difference in OFV for
  D-optimal, difference in gradient for ED-optimal)

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


FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)

Dtrace(fn="",
       it=1,
       ni=poped.db$design$ni,
       xtopt=poped.db$design$xt,
       xopt=poped.db$design$x,
       aopt=poped.db$design$a,
       gxt=0,ga=0,
       dmf=dmf,diff=3,
       ixt=FALSE,
       ia=FALSE, 
       itvector=NULL,
       dmfvector=NULL,
       poped.db,
       opt_xt=poped.db$settings$optsw[2],
       opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
       opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
       rsit=200)
#> RS - It. : 1   OFV : 1.14386e+24

Dtrace(fn="",
       it=1,
       ni=poped.db$design$ni,
       xtopt=poped.db$design$xt,
       xopt=poped.db$design$x,
       aopt=poped.db$design$a,
       gxt=0,ga=0,
       dmf=dmf,diff=3,
       ixt=FALSE,
       ia=FALSE, 
       itvector=NULL,
       dmfvector=NULL,
       poped.db,
       opt_xt=poped.db$settings$optsw[2],
       opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
       opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
       rsit=0)
#> SG - It. : 1  OFV : 1.144e+24   Diff. :     3

Dtrace(fn="",
       it=1,
       ni=poped.db$design$ni,
       xtopt=poped.db$design$xt,
       xopt=poped.db$design$x,
       aopt=poped.db$design$a,
       gxt=0,ga=0,
       dmf=dmf,diff=3,
       ixt=FALSE,
       ia=FALSE, 
       itvector=NULL,
       dmfvector=NULL,
       poped.db,
       opt_xt=poped.db$settings$optsw[2],
       opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
       opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
       rsit=1)
#> RS - It. : 1   OFV : 1.14386e+24
#> 
#> *******************************
#> RS Results
#>  OFV(mf) = 1.14386e+24
#> 
#> *********************************
#> 

Dtrace(fn="",
       it=1,
       ni=poped.db$design$ni,
       xtopt=poped.db$design$xt,
       xopt=poped.db$design$x,
       aopt=poped.db$design$a,
       gxt=0,ga=0,
       dmf=dmf,
       diff=0,
       ixt=FALSE,
       ia=FALSE, 
       itvector=NULL,
       dmfvector=NULL,
       poped.db,
       opt_xt=poped.db$settings$optsw[2],
       opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
       opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
       rsit=1)
#> RS - It. : 1   OFV : 1.14386e+24
#> 
#> *******************************
#> RS Results
#>  OFV(mf) = 1.14386e+24
#> 
#> *********************************
#> 

Dtrace(fn="",
       it=5,
       ni=poped.db$design$ni,
       xtopt=poped.db$design$xt,
       xopt=poped.db$design$x,
       aopt=poped.db$design$a,
       gxt=0,ga=0,
       dmf=dmf,
       diff=0,
       ixt=FALSE,
       ia=FALSE, 
       itvector=NULL,
       dmfvector=NULL,
       poped.db,
       opt_xt=poped.db$settings$optsw[2],
       opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
       opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
       rsit=1)
#> SG - It. : 4  OFV : 1.144e+24   Diff. :     0
#> 
#> SG - Iteration 4 --------- FINAL -------------------------
#> OFV(mf)    : 1.14386e+24
#> diff       : 0
#> *************************************************************
```
