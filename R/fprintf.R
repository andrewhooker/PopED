## Function written to match MATLAB function
## Author: Andrew Hooker

fprintf <- function(file="",...,both=FALSE){
  if(any(class(file)=="file")){
    cat(sprintf(...),file=file)
  } else {
    if(file==""){cat(sprintf(...))} else {
      cat(sprintf(file,...))
    }
  }
  if(both && file!=""){
    cat(sprintf(...))
  }
}
