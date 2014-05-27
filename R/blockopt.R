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
#' @param opt_method If "RS" (random search), "SG" (stochastic gradient) or "DO" (discrete optimization) then specifc output is produced.
#' 
#' @family Helper
#' @export
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_blockopt.R
# @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockopt <- function(fn,poped.db,opt_method=""){
  
  if(any(opt_method==c("RS","SG","DO"))){
    fprintf(fn,'==============================================================================\n')
    fprintf(fn,'Optimization Settings\n\n')
    if(opt_method=="RS"){
      fprintf(fn,'Random Search :\n')
      fprintf(fn,'Number of cycles : %g\n',poped.db$rsit)
      fprintf(fn,'Locality factor for xt : %g\n',poped.db$rslxt)
      fprintf(fn,'Locality factor for a  : %g\n',poped.db$rsla)
    }
    if(opt_method=="SG"){
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
    }
    if(opt_method=="DO"){    
      fprintf(fn,'Discrete Optimization  :\n')
      fprintf(fn,'RS int it: %g\n',poped.db$intrsit)
      fprintf(fn,'SG int it: %g\n',poped.db$intsgit)
    }
    fprintf(fn,"\n")
  }
  return( ) 
}