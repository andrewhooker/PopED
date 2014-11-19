## for translation
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/tools/matlab.to.r.R")
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/tools/function.to.sub.matrix.R")

foo <- matlab.to.r("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped/ofv/ofv_criterion.m")

foo <- matlab.to.r("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/PopED_Examples_matlab/Model Templates/function_input.m")
foo <- matlab.to.r("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/PopED_Examples_matlab/Model Templates/ff.m")
foo <- matlab.to.r("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/PopED_Examples_matlab/Model Templates/sfg.m")
foo <- matlab.to.r("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/PopED_Examples_matlab/Model Templates/feps.m")




source("function_input_fo_reduced_eval_pfim_way.R")

##tmp <- function.to.sub.matrix(file="../../pargen.R",mat="\\Wret")
source("../../pargen.R")
