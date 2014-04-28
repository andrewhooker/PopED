context("Parameter generation")

test_that("Means and variances are as expected", {
  
  source("examples_fcn_doc/warfarin_optimize.R")
  source("examples_fcn_doc/examples_pargen.R")
  
  expect_that(all(mean.diff.ln<5),is_true())
  expect_that(all(var.diff.ln[1:3]<10),is_true())
  expect_that(all(mean.diff.n<5),is_true())
  expect_that(all(var.diff.n[1:3]<10),is_true())
  expect_that(all(mean.diff.u<5),is_true())
  expect_that(all(range.diff.u[1:3]<10),is_true())
    
})

