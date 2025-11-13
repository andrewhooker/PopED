# Model linearization with respect to epsilon.

The function performs a linearization of the model with respect to the
residual variability. Derivative of model w.r.t. eps evaluated at eps=0

## Usage

``` r
LinMatrixH(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, poped.db)
```

## Arguments

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

- xt_ind:

  A vector of the individual/group sample times

- x:

  A matrix for the discrete design variables. Each row is a group.

- a:

  A matrix of covariates. Each row is a group.

- bpop:

  The fixed effects parameter values. Supplied as a vector.

- b_ind:

  vector of individual realization of the BSV terms b

- bocc_ind:

  Vector of individual realizations of the BOV terms bocc

- poped.db:

  A PopED database.

## Value

A matrix of size (samples per individual x number of epsilons)

## See also

Other FIM:
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
[`ofv_criterion()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_criterion.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)
