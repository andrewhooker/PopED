# tryCatch both warnings (with value) and errors

Catch \*and\* save both errors and warnings, and in the case of a
warning, also keep the computed result.

## Usage

``` r
tryCatch.W.E(expr)
```

## Arguments

- expr:

  an R expression to evaluate

## Value

a list with 'value' and 'warning', where 'value' may be an error caught.

## Author

Martin Maechler, The R Core Team
