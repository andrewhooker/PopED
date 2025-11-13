# Return all the unfixed parameters

all = vector of all unfixed params var derivative is a vector of 1 and
0, 1 means derivative of parameter is taken w.r.t. variance otherwise
w.r.t. sd If params is supplied then the parameter is taken from this
vector instead of poped.db

## Usage

``` r
get_unfixed_params(poped.db, params = NULL)
```

## Arguments

- poped.db:

  a PopED database.

- params:

  If params is supplied then the parameters are taken from this vector.

## Value

A list with the parameters. All unfixed parameters are also returned in
the "`all` output with the specified order
(bpop,d,covd,docc,covdocc,sigma,covsigma). `var_derivative` is a vector
of 1's or 0's, 1 means derivative of parameter is taken with respect to
the variance otherwise with respect to standard deviation.
