## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

blockheader_2 <- function(name,iter,globalStructure,
                          e_flag=FALSE,opt_xt=globalStructure$optsw[2],
                          opt_a=globalStructure$optsw[4],opt_x=globalStructure$optsw[4],
                          opt_samps=globalStructure$optsw[1],opt_inds=globalStructure$optsw[5],
                          fmf=0,dmf=0,bpop=NULL,d=NULL,docc=NULL,sigma=NULL,...)
{
  # BLOCKHEADER_2
  #   filename to write to is 
  #   globalStructure$strOutputFilePath,globalStructure$strOutputFileName,NAME,iter,globalStructure$strOutputFileExtension
  
  #   if((bDiscreteOpt)){
  #     tmpfile=sprintf('%s_Discrete_%g%s',globalStructure$strOutputFileName,iter,globalStructure$strOutputFileExtension)
  #   } else {
  #     tmpfile=sprintf('%s_RS_SG_%g%s',globalStructure$strOutputFileName,iter,globalStructure$strOutputFileExtension)
  #   }
  
  #tmpfile=sprintf('%s_%s_%g%s',globalStructure$strOutputFileName,name,iter,globalStructure$strOutputFileExtension)
  tmpfile=sprintf('%s_%s_%g.txt',globalStructure$strOutputFileName,name,iter)
  tmpfile = fullfile(globalStructure$strOutputFilePath,tmpfile)
  
  if(!(tmpfile=='')){
    fn=file(tmpfile,'w')
    if((fn==-1)){
      stop(sprintf('output file could not be opened'))
    }
  } else {
    filename=readline("File to open for output: ")
    fn = file(filename, 'w')
    if((fn == -1)){
      stop(sprintf('output file could not be opened'))
    }
  }
  
  tic()
  tic(name=".poped_total_time")
  
  # -------------- LOG FILE: initial status
  if(name=="RS"){
    alg_name <- "Adaptive Random Search"
    fprintf(fn,'PopED Optimization Results for the %s Algorithm \n\n',alg_name)
  }  else {
    fprintf(fn,'PopED Results \n\n')
  }
  fprintf(fn,'        ')
  fprintf(fn,datestr_poped(globalStructure$Engine$Type))
  fprintf(fn,'\n\n')
  
  blockexp(fn,globalStructure,
           e_flag=e_flag,opt_xt=opt_xt,
           opt_a=opt_a,opt_x=opt_x,
           opt_samps=opt_samps,opt_inds=opt_inds)
  
  if(dmf!=0 || fmf != 0) fprintf(fn,'===============================================================================\nInitial design evaluation\n')
  if(dmf!=0) fprintf(fn,'\ndet(FIM) = %g\n',dmf)
  
  if(any(fmf!=0)){
    param_vars=diag_matlab(inv(fmf))
    returnArgs <-  get_cv(param_vars,bpop,d,docc,sigma,globalStructure) 
    params <- returnArgs[[1]]
    param_cvs <- returnArgs[[2]]
    
    fprintf(fn,'\nEfficiency criterion det(FIM)^(1/npar) = %g\n',dmf^(1/length(params)))
    
    parnam <- get_parnam(globalStructure)
    fprintf(fn,'\nInitial design expected parameter variance and relative standard error (%sRSE)\n','%')
    fprintf('\nInitial design expected parameter variance and relative standard error (%sRSE)\n','%')
    df <- data.frame("Parameter"=parnam,"Values"=params, "Variance"=param_vars, "RSE"=t(param_cvs*100))
    print(df,digits=3, print.gap=3,row.names=F)
    capture.output(print(df,digits=3, print.gap=3,row.names=F),file=fn)
    fprintf('\n')
    fprintf(fn,'\n')
    
  }
  
  blockopt_2(fn,globalStructure,opt_method=name)
  blockother_2(fn,globalStructure)
  
  blockoptwrt(fn,globalStructure$optsw, opt_xt=opt_xt,
              opt_a=opt_a,opt_x=opt_x,
              opt_samps=opt_samps,opt_inds=opt_inds)
  
  return( fn) 
}
