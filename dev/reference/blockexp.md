# Summarize your experiment for optimization routines

Create some output to the screen and a text file that summarizes the
initial design and the design space you will use to optimize.

## Usage

``` r
blockexp(
  fn,
  poped.db,
  e_flag = FALSE,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[4],
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5]
)
```

## Arguments

- fn:

  The file handle to write to.

- poped.db:

  A PopED database.

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

## See also

Other Helper:
[`blockfinal()`](https://andrewhooker.github.io/PopED/dev/reference/blockfinal.md),
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


blockexp("",poped.db, opt_xt=TRUE)
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
#> Number of individuals: 32
#> Number of groups (individuals with same design): 1
#> Number of individuals per group:
#>  
#> Warning: 2 arguments not used by format '    Group %g: %g
#> '
#>     Group 1: 32
#> Number of samples per group:
#>  Number of discrete experimental variables: 0
#> Number of model covariates: 1
#> 
#> Initial Sampling Schedule
#> Group 1:    0.5      1      2      6     24     36     72    120
#> 
#> Minimum allowed sampling values
#> Group 1:   0.01   0.01   0.01   0.01   0.01   0.01   0.01   0.01
#> 
#> Maximum allowed sampling values
#> Group 1:    120    120    120    120    120    120    120    120
#> 
#> Covariates:
#> Group 1: 
#> Warning: 2 arguments not used by format '%g'
#> 70
#> 
#> NULL
```
