#---------- functions
dots <- function(...) {
  eval(substitute(alist(...)))
}