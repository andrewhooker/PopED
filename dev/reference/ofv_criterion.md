# Normalize an objective function by the size of the FIM matrix

Compute a normalized OFV based on the size of the FIM matrix. This value
can then be used in efficiency calculations. This is NOT the OFV used in
optimization, see
[`ofv_fim`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md).

## Usage

``` r
ofv_criterion(
  ofv_f,
  num_parameters,
  poped.db,
  ofv_calc_type = poped.db$settings$ofv_calc_type
)
```

## Arguments

- ofv_f:

  An objective function

- num_parameters:

  The number of parameters to use for normalization

- poped.db:

  a poped database

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

## Value

The specified criterion value.

## See also

Other FIM:
[`LinMatrixH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixH.md),
[`LinMatrixLH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixLH.md),
[`LinMatrixL_occ()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixL_occ.md),
[`calc_ofv_and_fim()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_fim.md),
[`ed_laplace_ofv()`](https://andrewhooker.github.io/PopED/dev/reference/ed_laplace_ofv.md),
[`ed_mftot()`](https://andrewhooker.github.io/PopED/dev/reference/ed_mftot.md),
[`efficiency()`](https://andrewhooker.github.io/PopED/dev/reference/efficiency.md),
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md),
[`evaluate.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md),
[`gradf_eps()`](https://andrewhooker.github.io/PopED/dev/reference/gradf_eps.md),
[`mf3()`](https://andrewhooker.github.io/PopED/dev/reference/mf3.md),
[`mf7()`](https://andrewhooker.github.io/PopED/dev/reference/mf7.md),
[`mftot()`](https://andrewhooker.github.io/PopED/dev/reference/mftot.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)

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


## evaluate initial design 
FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
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
get_rse(FIM,poped.db)
#>        CL         V        KA      d_CL       d_V      d_KA  sig_prop   sig_add 
#>  5.096246  3.031164 14.260384 29.761226 36.681388 26.748640 32.011719 25.637971 

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=1),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=1) # det(FIM)
#> [1] 1016.943

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=2),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=2) 
#> [1] 1.140916

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=4),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=4)
#> [1] 1016.943

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=6),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=6)
#> [1] 1.75168

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=7),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=7) 
#> [1] 0
```
