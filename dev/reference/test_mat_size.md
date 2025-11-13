# Test to make sure that matricies are the right size

Test to make sure that matricies are the right size

## Usage

``` r
test_mat_size(correct_size, mat, name)
```

## Arguments

- correct_size:

  the correct size of a matrix

- mat:

  The matrix to test.

- name:

  The name of the matrix as a string.

## Examples

``` r
test_mat_size(c(2,3),zeros(2,3),"foo")
#> [1] 1

if (FALSE) { # \dontrun{
  test_mat_size(c(2,3),zeros(2,6),"foo")
} # }

test_mat_size(c(1,3),c(2,6,7),"foo")
#> [1] 1


```
