context("Misc functions")

test_that("inv works", {
  mat <- diag(c(2,2,2))
  expect_equal(diag(c(1/2,1/2,1/2)),inv(mat))
  expect_equal(diag(c(1/2,1/2,1/2)),inv(mat,method=2))
  expect_equal(diag(c(1/2,1/2,1/2)),inv(mat,method=3))
  
  mat <- diag(c(2,2,2,0))
  expect_equal(diag(c(1/2,1/2,1/2,0)),inv(mat,pseudo_on_fail = T))
  expect_equal(diag(c(1/2,1/2,1/2,0)),inv(mat,method=3))
  expect_error(inv(mat,pseudo_on_fail = F))
})
