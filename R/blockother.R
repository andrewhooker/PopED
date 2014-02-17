## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockother <- function(fn,globalStructure){
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'Internal linearizations h :\n')
  fprintf(fn,'Model :\n')
  fprintf(fn,'hlf = %g\n',globalStructure$hlf)
  fprintf(fn,'hlg = %g\n',globalStructure$hlg)
  fprintf(fn,'hle = %g\n',globalStructure$hle)
  fprintf(fn,'F$I.M. :\n')
  fprintf(fn,'hm1 = %g\n',globalStructure$hm1)
  fprintf(fn,'hm2 = %g\n',globalStructure$hm2)
  fprintf(fn,'\n')
  fprintf(fn,'Gradients for optimization h :\n')
  fprintf(fn,'hgd = %g\n',globalStructure$hgd)
  fprintf(fn,'==============================================================================\n')
return( ) 
}
