## Function written to match MATLAB function
## Author: Andrew Hooker

#tic <- function(...){
#    return(matlab::tic())
#}

tic <- function (gcFirst = FALSE,name=".poped_savedTime") 
{
  if (gcFirst == TRUE) {
    gc(verbose = FALSE)
  }
  #if(!exists("savedTime",envir=.PopEDNamespaceEnv)){
  #.PopEDNamespaceEnv <- new.env(parent=baseenv())
  #}
  assign(name, proc.time()[3], pos = 1)
  invisible()
}