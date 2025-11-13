# Optimization main module for PopED

Optimize the objective function. The function works for both discrete
and continuous optimization variables. If more than one optimization
method is specified then the methods are run in series. If
`loop_methods=TRUE` then the series of optimization methods will be run
for `iter_max` iterations, or until the efficiency of the design after
the current series (compared to the start of the series) is less than
`stop_crit_eff`.

## Usage

``` r
poped_optim_3(
  poped.db,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5],
  method = c("ARS", "BFGS", "LS"),
  control = list(),
  trace = TRUE,
  fim.calc.type = poped.db$settings$iFIMCalculationType,
  ofv_calc_type = poped.db$settings$ofv_calc_type,
  ds_index = poped.db$parameters$ds_index,
  approx_type = poped.db$settings$iApproximationMethod,
  d_switch = poped.db$settings$d_switch,
  ED_samp_size = poped.db$settings$ED_samp_size,
  bLHS = poped.db$settings$bLHS,
  use_laplace = poped.db$settings$iEDCalculationType,
  out_file = "",
  parallel = F,
  parallel_type = NULL,
  num_cores = NULL,
  loop_methods = ifelse(length(method) > 1, TRUE, FALSE),
  iter_max = 10,
  stop_crit_eff = 1.001,
  stop_crit_diff = NULL,
  stop_crit_rel = NULL,
  ofv_fun = poped.db$settings$ofv_fun,
  maximize = T,
  allow_replicates = TRUE,
  allow_replicates_xt = TRUE,
  allow_replicates_a = TRUE,
  ...
)
```

## Arguments

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

- method:

  A vector of optimization methods to use in a sequential fashion.
  Options are `c("ARS","BFGS","LS","GA")`. `c("ARS")` is for Adaptive
  Random Search
  [`optim_ARS`](https://andrewhooker.github.io/PopED/dev/reference/optim_ARS.md).
  `c("LS")` is for Line Search
  [`optim_LS`](https://andrewhooker.github.io/PopED/dev/reference/optim_LS.md).
  `c("BFGS")` is for Method "L-BFGS-B" from
  [`optim`](https://rdrr.io/r/stats/optim.html). `c("GA")` is for the
  genetic algorithm from
  [`ga`](https://github.com/luca-scr/GA/reference/ga.html). If
  `opt_inds=TRUE` then this optimization is always added to the end of
  the sequential optimization.

- control:

  Contains control arguments specified for each method separately.

- trace:

  Should the algorithm output results intermittently.

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

- ds_index:

  Ds_index is a vector set to 1 if a parameter is uninteresting,
  otherwise 0. size=(1,num unfixed parameters). First unfixed bpop, then
  unfixed d, then unfixed docc and last unfixed sigma. Default is the
  fixed effects being important, everything else not important. Used in
  conjunction with `ofv_calc_type=6`.

- approx_type:

  Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI.

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

- out_file:

  Save output from the optimization to a file.

- parallel:

  Should we use parallel computations?

- parallel_type:

  Which type of parallelization should be used? Can be "snow" or
  "multicore". "snow" works on Linux-like systems & Windows. "multicore"
  works only on Linux-like systems. By default this is chosen for you
  depending on your operating system. See
  [`start_parallel`](https://andrewhooker.github.io/PopED/dev/reference/start_parallel.md).

- num_cores:

  The number of cores to use in the parallelization. By default is set
  to the number output from
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).
  See
  [`start_parallel`](https://andrewhooker.github.io/PopED/dev/reference/start_parallel.md).

- loop_methods:

  Should the optimization methods be looped for `iter_max` iterations,
  or until the efficiency of the design after the current series
  (compared to the start of the series) is less than, or equal to,
  `stop_crit_eff`?

- iter_max:

  If line search is used then the algorithm tests if line search (always
  run at the end of the optimization iteration) changes the design in
  any way. If not, the algorithm stops. If yes, then a new iteration is
  run unless `iter_max` iterations have already been run.

- stop_crit_eff:

  If `loop_methods==TRUE`, the looping will stop if the efficiency of
  the design after the current series (compared to the start of the
  series) is less than, or equal to, `stop_crit_eff` (if
  `maximize==FALSE` then 1/stop_crit_eff is the cut off and the
  efficiency must be greater than or equal to this value to stop the
  looping).

- stop_crit_diff:

  If `loop_methods==TRUE`, the looping will stop if the difference in
  criterion value of the design after the current series (compared to
  the start of the series) is less than, or equal to, `stop_crit_diff`
  (if `maximize==FALSE` then -stop_crit_diff is the cut off and the
  difference in criterion value must be greater than or equal to this
  value to stop the looping).

- stop_crit_rel:

  If `loop_methods==TRUE`, the looping will stop if the relative
  difference in criterion value of the design after the current series
  (compared to the start of the series) is less than, or equal to,
  `stop_crit_rel` (if `maximize==FALSE` then -stop_crit_rel is the cut
  off and the relative difference in criterion value must be greater
  than or equal to this value to stop the looping).

- ofv_fun:

  User defined function used to compute the objective function. The
  function must have a poped database object as its first argument and
  have "..." in its argument list. Can be referenced as a function or as
  a file name where the function defined in the file has the same name
  as the file. e.g. "cost.txt" has a function named "cost" in it.

- maximize:

  Should the objective function be maximized or minimized?

- allow_replicates:

  Should the algorithm allow parameters to have the same value?

- ...:

  arguments passed to other functions.

## Details

This function takes information from the PopED database supplied as an
argument. The PopED database supplies information about the the model,
parameters, design and methods to use. Some of the arguments coming from
the PopED database can be overwritten; if they are supplied then they
are used instead of the arguments from the PopED database.

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
[`poped_optimize()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optimize.md)
