# Date: February 2018
# Author: Giulia
# Task: Run "render" function to transform intro-poped.Rmd into a docunment "html_vignette"
# Output: intro-poped.html
# Note: to run this code it is necessary to Set the directory to Source file location (i.e ~/PopED/vignettes)
rm(list=ls(all=TRUE))
library(rmarkdown)
devtools::load_all("../")

render("intro-poped.Rmd")
viewer <- getOption("viewer")
viewer("intro-poped.html")

render("examples.Rmd")
viewer <- getOption("viewer")
viewer("examples.html")
