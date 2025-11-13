# The Fisher Information Matrix (FIM) for one individual

Compute the FIM for one individual given specific model(s), parameters,
design and methods.

## Usage

``` r
mf3(model_switch, xt, x, a, bpop, d, sigma, docc, poped.db)
```

## Arguments

- model_switch:

  A vector that is the same size as xt, specifying which model each
  sample belongs to.

- xt:

  A vector of sample times.

- x:

  A vector for the discrete design variables.

- a:

  A vector of covariates.

- bpop:

  The fixed effects parameter values. Supplied as a vector.

- d:

  A between subject variability matrix (OMEGA in NONMEM).

- sigma:

  A residual unexplained variability matrix (SIGMA in NONMEM).

- docc:

  A between occasion variability matrix.

- poped.db:

  A PopED database.

## Value

As a list:

- ret:

  The FIM for one individual

- poped.db:

  A PopED database

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
[`mf7()`](https://andrewhooker.github.io/PopED/dev/reference/mf7.md),
[`mftot()`](https://andrewhooker.github.io/PopED/dev/reference/mftot.md),
[`ofv_criterion()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_criterion.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)
