context("Design Evaluation")

test_that("Design evaluation works", {
  
  source("examples_fcn_doc/warfarin_basic.R")
  source("examples_fcn_doc/examples_evaluate_design.R")

  expect_output(str(evaluate_design(poped.db)),"List of 3")
})

test_that("Power evaluation works", {
  
  source("examples_fcn_doc/examples_evaluate_power.R")

  expect_error(evaluate_power(poped.db,bpopIdx = 4))

  expect_equal(evaluate_power(poped.db_2,bpopIdx = 4)$power$predPower,91.5)
    
})