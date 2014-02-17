## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockopt <- function(fn,globalStructure){
  
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'                     Optimization Settings : %s\n',globalStructure$opttit)
  fprintf(fn,'------------------------------------------------------------------------------\n')
  fprintf(fn,'Random Search :\n')
  fprintf(fn,'Number of cycles : %g\n',globalStructure$rsit)
  fprintf(fn,'Locality factor for xt : %g\n',globalStructure$rslxt)
  fprintf(fn,'Locality factor for a  : %g\n',globalStructure$rsla)
  fprintf(fn,'------------------------------------------------------------------------------\n')
  fprintf(fn,'Stochastic Gradient :\n')
  if((globalStructure$convergence_eps!=0)){
    fprintf(fn,'Maximum number of cycles : %g\n',globalStructure$sgit)
    fprintf(fn,'Epsilon for termination : %g\n',globalStructure$convergence_eps)
  } else {
    fprintf(fn,'Number of cycles : %g\n',globalStructure$sgit)
  }
  fprintf(fn,'First step factor for xt: %g\n', globalStructure$cfaxt)
  fprintf(fn,'First step factor for a: %g\n', globalStructure$cfaa)
  fprintf(fn,'RS m0it: %g\n',globalStructure$maxrsnullit)
  fprintf(fn,'------------------------------------------------------------------------------\n')
  fprintf(fn,'Discrete Optimization  :\n')
  fprintf(fn,'RS int it: %g\n',globalStructure$intrsit)
  fprintf(fn,'SG int it: %g\n',globalStructure$intsgit)
  fprintf(fn,'==============================================================================\n')
return( ) 
}
