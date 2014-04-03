#' Summarize your optimization settings for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the optimization settings
#' you will use to optimize.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams blockexp
#' 
#' @family Helper
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockopt <- function(fn,poped.db){
  
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'                     Optimization Settings : %s\n',poped.db$opttit)
  fprintf(fn,'------------------------------------------------------------------------------\n')
  fprintf(fn,'Random Search :\n')
  fprintf(fn,'Number of cycles : %g\n',poped.db$rsit)
  fprintf(fn,'Locality factor for xt : %g\n',poped.db$rslxt)
  fprintf(fn,'Locality factor for a  : %g\n',poped.db$rsla)
  fprintf(fn,'------------------------------------------------------------------------------\n')
  fprintf(fn,'Stochastic Gradient :\n')
  if((poped.db$convergence_eps!=0)){
    fprintf(fn,'Maximum number of cycles : %g\n',poped.db$sgit)
    fprintf(fn,'Epsilon for termination : %g\n',poped.db$convergence_eps)
  } else {
    fprintf(fn,'Number of cycles : %g\n',poped.db$sgit)
  }
  fprintf(fn,'First step factor for xt: %g\n', poped.db$cfaxt)
  fprintf(fn,'First step factor for a: %g\n', poped.db$cfaa)
  fprintf(fn,'RS m0it: %g\n',poped.db$maxrsnullit)
  fprintf(fn,'------------------------------------------------------------------------------\n')
  fprintf(fn,'Discrete Optimization  :\n')
  fprintf(fn,'RS int it: %g\n',poped.db$intrsit)
  fprintf(fn,'SG int it: %g\n',poped.db$intsgit)
  fprintf(fn,'==============================================================================\n')
return( ) 
}
