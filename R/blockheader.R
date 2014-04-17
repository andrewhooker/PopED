#' Header function for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the problem you are tying to solve.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams blockheader_2
#' @inheritParams create.poped.database 
#' @param bDiscreteOpt is discrete optimization being performed?
#' @family Helper
#' @export
# @keywords internal
#' 
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockheader <- function(bDiscreteOpt,iter,poped.db){

if((bDiscreteOpt)){
    tmpfile=sprintf('%s_Discrete_%g%s',poped.db$strOutputFileName,iter,poped.db$strOutputFileExtension)
} else {
    tmpfile=sprintf('%s_RS_SG_%g%s',poped.db$strOutputFileName,iter,poped.db$strOutputFileExtension)
}
tmpfile = fullfile(poped.db$strOutputFilePath,tmpfile)

if(!(tmpfile=='')){
    fn=file(tmpfile,'w')
    if((fn==-1)){
        stop(sprintf('output file could not be opened'))
    }
} else {
  filename= readline("File to open for output: ")
  fn = file(filename, 'w')
  if((fn == -1)){
    stop(sprintf('output file could not be opened'))
  }
}

tic()

# -------------- LOG FILE: initial status
fprintf(fn,'        PopED          Optimization Results \n\n')
fprintf(fn,'        ')
fprintf(fn,datestr_poped(poped.db$Engine$Type))
fprintf(fn,'\n\n')

blockexp(fn,poped.db)
blockopt(fn,poped.db)
blockother(fn,poped.db)

blockoptwrt(fn,poped.db$optsw)

fprintf(fn,'=====================================================================\n')


return( fn) 
}
