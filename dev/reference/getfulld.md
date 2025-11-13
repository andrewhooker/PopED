# Create a full D (between subject variability) matrix given a vector of variances and covariances. Note, this does not test matching vector lengths.

Create a full D (between subject variability) matrix given a vector of
variances and covariances. Note, this does not test matching vector
lengths.

## Usage

``` r
getfulld(variance_vector, covariance_vector = NULL)
```

## Arguments

- variance_vector:

  The vector of the variances.

- covariance_vector:

  A vector of the covariances. Written in column major order for the
  lower triangular matrix.

## Value

The full matrix of variances for the between subject variances

## Examples

``` r
getfulld(c(1,2,3))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    2    0
#> [3,]    0    0    3

getfulld(c(1,2,3),c(7,6,5))
#>      [,1] [,2] [,3]
#> [1,]    1    7    6
#> [2,]    7    2    5
#> [3,]    6    5    3
```
