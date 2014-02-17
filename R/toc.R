## Function written to match MATLAB function
## Author: Andrew Hooker

#toc <- function(...){
#    return(matlab::toc(echo=FALSE))
#}

toc <- function (echo = TRUE,name=".poped_savedTime") 
{
  prevTime <- get(name, pos = 1)
  diffTimeSecs <- proc.time()[3] - prevTime
  if (echo) {
    cat(sprintf("elapsed time is %f seconds", diffTimeSecs), 
        "\n")
    return(invisible(diffTimeSecs))
  }
  else {
    return(diffTimeSecs)
  }
}
