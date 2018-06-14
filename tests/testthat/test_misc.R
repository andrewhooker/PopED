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

test_that("bounding and unbound of parameters works",{
  expect_equal(unbound_par(0.5,0,1),0)
  expect_equal(unbound_par(0,0,1,tol=0),-Inf)
  expect_equal(unbound_par(1,0,1,tol=0),Inf)
  expect_equal(unbound_par(0.025,0,1,logit = F),stats::qnorm(0.025))
  
  expect_equal(bound_par(0,0,1),0.5)
  expect_equal(bound_par(-Inf,0,1),0)
  expect_equal(bound_par(Inf,0,1),1)
  
  expect_equal(bound_par(unbound_par(0,0,1),0,1),1e-8)
  expect_equal(bound_par(unbound_par(0,0,1,tol = 0),0,1),0)
  
  
})