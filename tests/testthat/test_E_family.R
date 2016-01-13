context("E-family")

test_that("ed_laplace_ofv works", {
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_ed_laplace_ofv.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string))
  sink()
  file.remove("tmp.txt")

  expect_gt(output$E_ofv,1.0e+24)
    
})