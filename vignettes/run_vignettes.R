# Date: February 2018
# Author: Giulia
# Task: Run "render" function to transform intro-poped.Rmd into a docunment "html_vignette"
# Output: intro-poped.html
# Note: to run this code it is necessary to Set the directory to Source file location (i.e ~/PopED/vignettes)

if(interactive()) setwd(here::here("vignettes"))
rm(list=ls(all=TRUE))
library(rmarkdown)
devtools::load_all("../")
viewer <- getOption("viewer")

render("intro-poped.Rmd")
viewer("intro-poped.html")

render("examples.Rmd")
viewer("examples.html")

#render("model_def_other_pkgs.Rmd")
#viewer("model_def_other_pkgs.html")
