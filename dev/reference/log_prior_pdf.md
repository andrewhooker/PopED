# Compute the natural log of the PDF for the parameters in an E-family design

Compute the natural log of the PDF for the parameters in an E-family
design

## Usage

``` r
log_prior_pdf(
  alpha,
  bpopdescr,
  ddescr,
  return_gradient = F,
  return_hessian = F
)
```

## Arguments

- alpha:

  A parameter vector.

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

- return_gradient:

  Should the gradient be returned.

- return_hessian:

  Should the hessian be returned?
