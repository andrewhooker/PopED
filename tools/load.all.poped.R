library(devtools)
load_all_poped <- function(...){
  load_all("/Users/ahooker/Documents/_PROJECTS/PopED/repos/PopED/",...)
}

library(testthat)
test_all_poped  <- function(...){
  test_dir("/Users/ahooker/Documents/_PROJECTS/PopED/repos/PopED/tests/testthat/",...)
}