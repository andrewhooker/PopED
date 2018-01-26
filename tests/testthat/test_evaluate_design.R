context("Design Evaluation")

test_that("Design evaluation works", {
  
  source("examples_fcn_doc/warfarin_basic.R")
  source("examples_fcn_doc/examples_evaluate_design.R")

  expect_output(str(evaluate_design(poped.db)),"List of 3")
})

