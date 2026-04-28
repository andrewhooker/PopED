# Evaluate the Fisher Information Matrix (FIM)

Compute the FIM given the model, parameters, design and methods defined
in the PopED database. Some of the arguments coming from the PopED
database can be overwritten; by default these arguments are `NULL` in
the function, if they are supplied then they are used instead of the
arguments from the PopED database.

## Usage

``` r
evaluate.fim(
  poped.db,
  fim.calc.type = NULL,
  approx.method = NULL,
  FOCE.num = NULL,
  bpop.val = NULL,
  d_full = NULL,
  docc_full = NULL,
  sigma_full = NULL,
  model_switch = NULL,
  ni = NULL,
  xt = NULL,
  x = NULL,
  a = NULL,
  groupsize = NULL,
  deriv.type = NULL,
  ...
)
```

## Arguments

- poped.db:

  A PopED database.

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

- approx.method:

  Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI

- FOCE.num:

  Number individuals in each step of FOCE approximation method

- bpop.val:

  The fixed effects parameter values. Supplied as a vector.

- d_full:

  A between subject variability matrix (OMEGA in NONMEM).

- docc_full:

  A between occasion variability matrix.

- sigma_full:

  A residual unexplained variability matrix (SIGMA in NONMEM).

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

- groupsize:

  A vector of the number of individuals in each group.

- deriv.type:

  A number indicating the type of derivative to use:

  - 0=Complex difference

  - 1=Central difference

  - 20=Analytic derivative (placeholder)

  - 30=Automatic differentiation (placeholder)

- ...:

  Other arguments passed to the function.

## Value

The FIM.

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
[`gradf_eps()`](https://andrewhooker.github.io/PopED/dev/reference/gradf_eps.md),
[`mf3()`](https://andrewhooker.github.io/PopED/dev/reference/mf3.md),
[`mf7()`](https://andrewhooker.github.io/PopED/dev/reference/mf7.md),
[`mftot()`](https://andrewhooker.github.io/PopED/dev/reference/mftot.md),
[`ofv_criterion()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_criterion.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)

Other evaluate_design:
[`evaluate_design()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate_design.md),
[`evaluate_power()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate_power.md),
[`get_rse()`](https://andrewhooker.github.io/PopED/dev/reference/get_rse.md),
[`model_prediction()`](https://andrewhooker.github.io/PopED/dev/reference/model_prediction.md),
[`plot_efficiency_of_windows()`](https://andrewhooker.github.io/PopED/dev/reference/plot_efficiency_of_windows.md),
[`plot_model_prediction()`](https://andrewhooker.github.io/PopED/dev/reference/plot_model_prediction.md)

Other evaluate_FIM:
[`calc_ofv_and_fim()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_fim.md),
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)

## Examples

``` r
## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.md.CL
#> function (model_switch, xt, parameters, poped.db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         N = floor(xt/TAU) + 1
#>         y = (DOSE * Favail/V) * (KA/(KA - CL/V)) * (exp(-CL/V * 
#>             (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - 
#>             exp(-CL/V * TAU)) - exp(-KA * (xt - (N - 1) * TAU)) * 
#>             (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
#>         return(list(y = y, poped.db = poped.db))
#>     })
#> }
#> <bytecode: 0x55bedb304e30>
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
poped.db <- create.poped.database(ff_fun = ff.PK.1.comp.oral.sd.CL,
                                  fg_fun = sfg,
                                  fError_fun = feps.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  # notfixed_bpop=c(1,1,1,0),
                                  notfixed_bpop=c(CL=1,V=1,KA=1,Favail=0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=0.01,
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70)


## evaluate initial design with the reduced FIM
FIM.1 <- evaluate.fim(poped.db) 
FIM.1
#>             [,1]       [,2]      [,3]         [,4]         [,5]        [,6]
#> [1,] 19821.28445 -21.836551 -8.622140 0.000000e+00     0.000000  0.00000000
#> [2,]   -21.83655  20.656071 -1.807099 0.000000e+00     0.000000  0.00000000
#> [3,]    -8.62214  -1.807099 51.729039 0.000000e+00     0.000000  0.00000000
#> [4,]     0.00000   0.000000  0.000000 3.107768e+03    10.728786  0.02613561
#> [5,]     0.00000   0.000000  0.000000 1.072879e+01 27307.089308  3.26560786
#> [6,]     0.00000   0.000000  0.000000 2.613561e-02     3.265608 41.81083599
#> [7,]     0.00000   0.000000  0.000000 5.215403e+02 11214.210707 71.08763902
#>              [,7]
#> [1,]      0.00000
#> [2,]      0.00000
#> [3,]      0.00000
#> [4,]    521.54030
#> [5,]  11214.21071
#> [6,]     71.08764
#> [7,] 806176.95068
det(FIM.1)
#> [1] 5.996147e+22
det(FIM.1)^(1/7)
#> [1] 1794.658
get_rse(FIM.1,poped.db)
#>         CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] 
#>   4.738266   2.756206  13.925829  25.627205  30.344316  25.777327  11.170784 

## evaluate initial design with the full FIM
FIM.0 <- evaluate.fim(poped.db,fim.calc.type=0) 
FIM.0
#>               [,1]         [,2]        [,3]          [,4]         [,5]
#> [1,]  47625.234145  -341.996566   35.504624 -2.073844e+03 -5899.486674
#> [2,]   -341.996566    30.887205  -12.589615 -1.686280e+01   -54.629529
#> [3,]     35.504624   -12.589615  452.758773 -8.336530e-01   -43.619195
#> [4,]  -2073.844369   -16.862802   -0.833653  3.107768e+03    10.728786
#> [5,]  -5899.486674   -54.629529  -43.619195  1.072879e+01 27307.089308
#> [6,]      4.490538    -6.550313   18.653863  2.613561e-02     3.265608
#> [7,] -54419.723543 -1070.933661 2955.924225  5.215403e+02 11214.210707
#>             [,6]         [,7]
#> [1,]  4.49053810 -54419.72354
#> [2,] -6.55031322  -1070.93366
#> [3,] 18.65386273   2955.92422
#> [4,]  0.02613561    521.54030
#> [5,]  3.26560786  11214.21071
#> [6,] 41.81083599     71.08764
#> [7,] 71.08763902 806176.95068
det(FIM.0)
#> [1] 1.220371e+24
det(FIM.0)^(1/7)
#> [1] 2760.117
get_rse(FIM.0,poped.db)
#>         CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] 
#>   3.560994   2.560413   4.811952  26.270324  30.901555  26.503936  12.409516 

## evaluate initial design with the reduced FIM 
## computing all derivatives with respect to the 
## standard deviation of the residual unexplained variation 
FIM.4 <- evaluate.fim(poped.db,fim.calc.type=4) 
FIM.4
#>             [,1]       [,2]      [,3]         [,4]         [,5]        [,6]
#> [1,] 19821.28445 -21.836551 -8.622140 0.000000e+00     0.000000  0.00000000
#> [2,]   -21.83655  20.656071 -1.807099 0.000000e+00     0.000000  0.00000000
#> [3,]    -8.62214  -1.807099 51.729039 0.000000e+00     0.000000  0.00000000
#> [4,]     0.00000   0.000000  0.000000 3.107768e+03    10.728786  0.02613561
#> [5,]     0.00000   0.000000  0.000000 1.072879e+01 27307.089308  3.26560786
#> [6,]     0.00000   0.000000  0.000000 2.613561e-02     3.265608 41.81083599
#> [7,]     0.00000   0.000000  0.000000 1.043081e+02  2242.842141 14.21752780
#>             [,7]
#> [1,]     0.00000
#> [2,]     0.00000
#> [3,]     0.00000
#> [4,]   104.30806
#> [5,]  2242.84214
#> [6,]    14.21753
#> [7,] 32247.07803
det(FIM.4)
#> [1] 2.398459e+21
get_rse(FIM.4,poped.db,fim.calc.type=4)
#>         CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] 
#>   4.738266   2.756206  13.925829  25.627205  30.344316  25.777327   5.585392 

## evaluate initial design with the full FIM with A,B,C matricies
## should give same answer as fim.calc.type=0
FIM.5 <- evaluate.fim(poped.db,fim.calc.type=5) 
FIM.5
#>               [,1]         [,2]        [,3]          [,4]         [,5]
#> [1,]  47625.234145  -341.996566   35.504624 -2.073844e+03 -5899.486674
#> [2,]   -341.996566    30.887205  -12.589615 -1.686280e+01   -54.629529
#> [3,]     35.504624   -12.589615  452.758773 -8.336530e-01   -43.619195
#> [4,]  -2073.844369   -16.862802   -0.833653  3.107768e+03    10.728786
#> [5,]  -5899.486674   -54.629529  -43.619195  1.072879e+01 27307.089308
#> [6,]      4.490538    -6.550313   18.653863  2.613561e-02     3.265608
#> [7,] -54419.723543 -1070.933661 2955.924225  5.215403e+02 11214.210707
#>             [,6]         [,7]
#> [1,]  4.49053810 -54419.72354
#> [2,] -6.55031322  -1070.93366
#> [3,] 18.65386273   2955.92422
#> [4,]  0.02613561    521.54030
#> [5,]  3.26560786  11214.21071
#> [6,] 41.81083599     71.08764
#> [7,] 71.08763902 806176.95068
det(FIM.5)
#> [1] 1.220371e+24
get_rse(FIM.5,poped.db,fim.calc.type=5)
#>         CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] 
#>   3.560994   2.560413   4.811952  26.270324  30.901555  26.503936  12.409516 

## evaluate initial design with the reduced FIM with 
## A,B,C matricies and derivative of variance
## should give same answer as fim.calc.type=1 (default)
FIM.7 <- evaluate.fim(poped.db,fim.calc.type=7) 
FIM.7
#>             [,1]       [,2]      [,3]         [,4]         [,5]        [,6]
#> [1,] 19821.28445 -21.836551 -8.622140 0.000000e+00     0.000000  0.00000000
#> [2,]   -21.83655  20.656071 -1.807099 0.000000e+00     0.000000  0.00000000
#> [3,]    -8.62214  -1.807099 51.729039 0.000000e+00     0.000000  0.00000000
#> [4,]     0.00000   0.000000  0.000000 3.107768e+03    10.728786  0.02613561
#> [5,]     0.00000   0.000000  0.000000 1.072879e+01 27307.089308  3.26560786
#> [6,]     0.00000   0.000000  0.000000 2.613561e-02     3.265608 41.81083599
#> [7,]     0.00000   0.000000  0.000000 5.215403e+02 11214.210707 71.08763902
#>              [,7]
#> [1,]      0.00000
#> [2,]      0.00000
#> [3,]      0.00000
#> [4,]    521.54030
#> [5,]  11214.21071
#> [6,]     71.08764
#> [7,] 806176.95068
det(FIM.7)
#> [1] 5.996147e+22
get_rse(FIM.7,poped.db,fim.calc.type=7)
#>         CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] 
#>   4.738266   2.756206  13.925829  25.627205  30.344316  25.777327  11.170784 

## evaluate FIM and rse with prior FIM.1
poped.db.prior = create.poped.database(poped.db, prior_fim = FIM.1)
FIM.1.prior <- evaluate.fim(poped.db.prior)
all.equal(FIM.1.prior,FIM.1) # the FIM is only computed from the design in the poped.db
#> [1] TRUE
get_rse(FIM.1.prior,poped.db.prior) # the RSE is computed with the prior information
#>         CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] 
#>   3.350460   1.948932   9.847048  18.121170  21.456671  18.227322   7.898937 
```
