## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockheader <- function(bDiscreteOpt,iter,globalStructure){

if((bDiscreteOpt)){
    tmpfile=sprintf('%s_Discrete_%g%s',globalStructure$strOutputFileName,iter,globalStructure$strOutputFileExtension)
} else {
    tmpfile=sprintf('%s_RS_SG_%g%s',globalStructure$strOutputFileName,iter,globalStructure$strOutputFileExtension)
}
tmpfile = fullfile(globalStructure$strOutputFilePath,tmpfile)

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
fprintf(fn,datestr_poped(globalStructure$Engine$Type))
fprintf(fn,'\n\n')

blockexp(fn,globalStructure)
blockopt(fn,globalStructure)
blockother(fn,globalStructure)

blockoptwrt(fn,globalStructure$optsw)

fprintf(fn,'=====================================================================\n')


return( fn) 
}
