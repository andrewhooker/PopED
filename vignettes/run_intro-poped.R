# Date: February 2018
# Author: Giulia
# Task: Run "render" function to transform intro-poped.Rmd into a docunment "html_vignette"
# Output: intro-poped.html
# Note: to run this code it is necessary to be in the same directory of the Rmd file
# that we want to execute, in this case it is: ~/PopED/vignettes
rm(list=ls(all=TRUE))
library(rmarkdown)

## Ensure that no library from the user home will be loaded such that
## we run with the production packages only
.libPaths(grep("home", .libPaths(), value=TRUE, invert=TRUE))
render("intro-poped.Rmd")
