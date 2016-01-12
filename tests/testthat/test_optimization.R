context("Optimization")

test_that("optim_ARS works", {
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_optim_ARS.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string))
  sink()
  file.remove("tmp.txt")
  
  expect_true(res1$ofv <= 159.0012) # check that things are minimizing
  
  # check that box constraints hold
  expect_true(all(res_box$par<=4)) 
  expect_true(all(res_box$par >= 2))
              
  # check that there are only these allowed values in the solution
  expect_true(res_int$par %in% seq(-50,100,by=1))
  
  # check that function in maximizing
  expect_true(res_max$ofv >= 0)
  
})

test_that("a_line_search, mfea, poped_optimize and poped_optim work", {
  
  ex_string_1 <- ex_to_string("examples_fcn_doc/warfarin_optimize.R",comment_dontrun=comment_dontrun)
  eval(parse(text=ex_string_1))
  
  expect_warning(a_line_search(poped.db))
  
  ex_string_2 <- ex_to_string("examples_fcn_doc/examples_a_line_search.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string_2))
  sink()
  file.remove("tmp.txt")
  
  expect_true(output$best_changed)
  
  ex_string_3 <- ex_to_string("examples_fcn_doc/examples_mfea.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string_3))
  sink()
  file.remove("tmp.txt")
  
  expect_equivalent(out_1$poped.db$design$a[1,1],100)
  
  ex_string_4 <- ex_to_string("examples_fcn_doc/examples_poped_optimize.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string_4))
  sink()
  file.remove("tmp.txt")
  
  expect_equivalent(out_1$poped.db$design$a[1,1],100)
  
  ex_string_5 <- ex_to_string("examples_fcn_doc/examples_poped_optim.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string_5))
  sink()
  file.remove("tmp.txt")
  
  expect_equivalent(out_1$poped.db$design$a[1,1],100)
  
})