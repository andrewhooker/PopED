context("Optimization")

test_that("optim_ARS works", {
  
  if(skip_optim) skip("Optimization skipped")
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_optim_ARS.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")  
  sink()
  file.remove("tmp.txt","temp.R")
  
  
  expect_true(res1$ofv <= 159.0012) # check that things are minimizing
  
  # check that box constraints hold
  expect_true(all(res_box$par<=4)) 
  expect_true(all(res_box$par >= 2))
              
  # check that there are only these allowed values in the solution
  # expect_true(res_int$par %in% seq(-50,100,by=1))
  
  # check that function in maximizing
  expect_true(res_max$ofv >= 0)
  
})

test_that("a_line_search, mfea, poped_optimize, poped_optim and RS_opt_gen work", {
  
  if(skip_optim) skip("Optimization skipped")
  
  source("examples_fcn_doc/warfarin_optimize.R")

  expect_warning(a_line_search(poped.db))
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_a_line_search.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")  
  sink()
  file.remove("tmp.txt","temp.R")
  
  expect_true(output$best_changed)
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_mfea.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")  
  sink()
  file.remove("tmp.txt","temp.R")
  
  expect_equivalent(out_1$poped.db$design$a[1,1],100)
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_poped_optimize.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")  
  sink()
  file.remove("tmp.txt","temp.R")
  
  expect_equivalent(out_1$poped.db$design$a[1,1],100)
  #expect_equivalent(out_2$poped.db$design$a[1,1],100)
  #expect_equivalent(out_3$poped.db$design$a[1,1],0.01)
  
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_poped_optim.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")  
  sink()
  file.remove("tmp.txt","temp.R")
  
  expect_equivalent(out_1$poped.db$design$a[1,1],100)
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_RS_opt.R",comment_dontrun=comment_dontrun)
  writeLines(ex_string, "temp.R")
  sink("tmp.txt")
  source("temp.R")
  sink()
  file.remove("tmp.txt","temp.R")

  expect_output(str(out_1),"List of 6")
  
})