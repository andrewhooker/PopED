#' Result function for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the problem you solved.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams blockexp
#' 
#' @family Helper
# @export
# @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockfinal <- function(fn,fmf,dmf,groupsize,ni,xt,x,a,bpop,d,docc,sigma,m,poped.db){

fprintf(fn,'FINAL:===============================================================================\n')
fprintf(fn,'Individuals per group :\n')
fprintf(fn,'%g ',groupsize)
fprintf(fn,'\n')
if((poped.db$optsw[2]==TRUE)){
    fprintf(fn,'xt :\n')
    #gxt(1:m,1:maxni)=xt
    for(i in 1:m){
        for(j in 1:ni[i]){
            fprintf(fn,'%g ',xt[i,j])
        }
        fprintf(fn,'\n')
    }
}
if((poped.db$optsw[3]==TRUE)){
    fprintf(fn,'x :\n')
    fprintf(fn,'%g \n',x)
}
if((poped.db$optsw[4]==TRUE)){
    fprintf(fn,'a :\n')
    fprintf(fn,'%g \n',a)
}
if((poped.db$d_switch==TRUE)){
    fprintf(fn,'mf : \n')
    write_matrix(fn,fmf)
    fprintf(fn,'\nInverse(mf):\n')
    write_matrix(fn,inv(fmf))
}
fprintf(fn,'\ndet(mf) = %g\n',dmf)
fprintf(fn,'========================================================================================\n')
fprintf(fn,'Benchmark section\n')
time_value=toc()
fprintf(fn,'Total running time: %g seconds\n',time_value)
fprintf(fn,'========================================================================================\n')

fprintf(fn,'\n')
fprintf(fn,'========================================================================================\n')
fprintf(fn,'========================================================================================\n')
fprintf(fn,'\n')
fprintf(fn,'Variances and CVs\n')
fprintf(fn,'\n')

fprintf(fn,'========================================================================================\n')
fprintf(fn,'\n')
fprintf(fn,'Parameters\n')
param_vars=diag_matlab(inv(fmf))
 returnArgs <-  get_cv(param_vars,bpop,d,docc,sigma,poped.db) 
params <- returnArgs[[1]]
 param_cvs <- returnArgs[[2]]
if((poped.db$Engine$Type==1) ){#Matlab
   fprintf(fn,'%g\n',params)
} else {
   for(i in 1:length(params)){
      fprintf(fn,'%g ',params[i])
   }
   fprintf(fn,'\n')
}
fprintf(fn,'\n')

fprintf(fn,'========================================================================================\n')
fprintf(fn,'\n')
fprintf(fn,'Parameter variances\n')
param_vars=diag_matlab(inv(fmf))
if((poped.db$Engine$Type==1) ){#Matlab
   fprintf(fn,'%g\n',param_vars)
} else {
   for(i in 1:length(param_vars)){
      fprintf(fn,'%g ',param_vars[i])
   }
   fprintf(fn,'\n')
}
fprintf(fn,'\n')

fprintf(fn,'========================================================================================\n')
fprintf(fn,'\n')
fprintf(fn,'Parameter CVs\n')
if((poped.db$Engine$Type==1) ){#Matlab
   fprintf(fn,'%g\n',param_cvs)
} else {
  for(i in 1:length(param_cvs)){
    fprintf(fn,'%g ',param_cvs(i))
  }
  fprintf(fn,'\n')
}
fprintf(fn,'\n')
fprintf(fn,'========================================================================================\n')


return( ) 
}
