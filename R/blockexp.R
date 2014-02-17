## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockexp <- function(fn,globalStructure,e_flag=FALSE,
                     opt_xt=globalStructure$optsw[2],opt_a=globalStructure$optsw[4],opt_x=globalStructure$optsw[4],
                     opt_samps=globalStructure$optsw[1],opt_inds=globalStructure$optsw[5]){
  
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'Model description : %s \n',globalStructure$modtit)
  fprintf(fn,'\n')
  fprintf(fn,'Model Sizes : \n')
  fprintf(fn,'Number of individual model parameters                  g[j]    : Ng    = %g\n',globalStructure$ng)
  fprintf(fn,'Number of population model fixed parameters            bpop[j] : Nbpop = %g\n',globalStructure$nbpop)
  fprintf(fn,'Number of population model random effects parameters   b[j]    : Nb    = %g\n',globalStructure$NumRanEff)
  fprintf(fn,'\n')
  print_params(globalStructure$gbpop,"bpop",fn=fn,globalStructure=globalStructure, 
               head_txt="Typical Population Parameters",e_flag=e_flag)
  fprintf(fn,'\n')
  if((globalStructure$NumRanEff!=0)){
    fprintf(fn,"Between Subject Variability matrix D (variance units) \n")
    d=getfulld(globalStructure$gd[,2,drop=F],globalStructure$covd)
    MASS::write.matrix(d,file=fn)
    fprintf(fn,'\n')
    print_params(globalStructure$gd,"D",fn=fn,globalStructure=globalStructure,param_sqrt=TRUE, matrix_elements=T,
                 head_txt="Diagonal Elements of D",e_flag=e_flag)
    fprintf(fn,'\n')
  }
  
  docc_full = getfulld(poped.db$docc[,2,drop=F],poped.db$covdocc)
  
  fprintf(fn,'Residual Unexplained Variability matrix SIGMA (variance units) : \n')
  sigma = globalStructure$sigma
  MASS::write.matrix(sigma,file=fn)
  fprintf(fn,'\n')
  
  #sigma_d = diag_matlab(globalStructure$sigma)
  sigma_d <- cbind(c(0,0),diag_matlab(globalStructure$sigma),c(1,1))
  print_params(sigma_d,"SIGMA",fn=fn,globalStructure=globalStructure,param_sqrt=TRUE, matrix_elements=T,
               head_txt="Diagonal Elements of SIGMA",e_flag=e_flag)
  fprintf(fn,'\n')
  
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'Experiment description (design and design space)\n')
  fprintf(fn,'\n')
  
  tmp_txt <- "Numer of individuals"
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(min, max)',sep=" ")
  tmp_txt <- paste(tmp_txt,': %g',sep="")
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
  tmp_txt <- paste(tmp_txt,'\n',sep="")
  fprintf(fn,tmp_txt,sum(globalStructure$groupsize),globalStructure$design$mintotgroupsize,globalStructure$design$maxtotgroupsize)
  
  fprintf(fn,'Number of groups (individuals with same design): %g\n',globalStructure$m)
  
  tmp_txt <- "Numer of individuals per group"
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(min, max)',sep=" ")
  tmp_txt <- paste(tmp_txt,':\n',sep="")
  fprintf(fn,tmp_txt)
  
  fprintf(fn," ")
  tmp_txt <- '    Group %g: %g'
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
  tmp_txt <- paste(tmp_txt,'\n',sep="")
  fprintf(fn,tmp_txt,1:globalStructure$m,globalStructure$groupsize,globalStructure$mingroupsize, globalStructure$maxgroupsize)
  
  tmp_txt <- "Numer of samples per group"
  if(opt_samps) tmp_txt <- paste(tmp_txt,'(min, max)',sep=" ")
  tmp_txt <- paste(tmp_txt,':\n',sep="")
  fprintf(fn,tmp_txt)
  
  fprintf(fn," ")
  tmp_txt <- '    Group %g: %g'
  if(opt_samps) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
  tmp_txt <- paste(tmp_txt,'\n',sep="")
  fprintf(fn,tmp_txt,1:globalStructure$m,globalStructure$gni,globalStructure$minni, globalStructure$maxni)
  
  fprintf(fn,'Number of discrete experimental variables: %g\n',globalStructure$nx)
  fprintf(fn,'Number of model covariates: %g\n',globalStructure$na)
  
  fprintf(fn,'\n')
  
  print_xt(globalStructure$gxt,globalStructure$gni,poped.db$global_model_switch,fn,
           head_txt="Initial Sampling Schedule\n")
  fprintf(fn,'\n')
  if(opt_xt){
    print_xt(globalStructure$gxt,globalStructure$gni,poped.db$global_model_switch,fn,
          head_txt="Minimum allowed sampling values\n",xt_other=globalStructure$gminxt)
    fprintf(fn,'\n')
    print_xt(globalStructure$gxt,globalStructure$gni,poped.db$global_model_switch,fn,
             head_txt="Maximum allowed sampling values\n",xt_other=globalStructure$gmaxxt)
    fprintf(fn,'\n')
  }  
  
  
  if((globalStructure$nx!=0)){
    tmp_txt <- "Discrete Variables"
    if(opt_x) tmp_txt <- paste(tmp_txt,' (possible vales)',sep=" ")
    tmp_txt <- paste(tmp_txt,':\n',sep="")
    fprintf(fn,tmp_txt)
    for(ct1 in 1:globalStructure$m){
      fprintf(fn,'Group %g: ', ct1)
      for(ct2 in 1:globalStructure$nx){
        tmp_txt <- '%g'
        if(opt_x) tmp_txt <- paste(tmp_txt,'(%s)',sep=" ")
        if(ct2<globalStructure$nx) tmp_txt <- paste(tmp_txt,' : ',sep="")        
        discrete_val = globalStructure$discrete_x[[ct1,ct2]]  
        fprintf(fn,tmp_txt,globalStructure$gx[ct1,ct2],get_vector_str(discrete_val))
      }
      fprintf(fn,'\n')
    }
    fprintf(fn,'\n')
  }
  
  
  if((globalStructure$na!=0)){   
    tmp_txt <- "Covariates"
    if(opt_a) tmp_txt <- paste(tmp_txt,' (min, max)',sep=" ")
    tmp_txt <- paste(tmp_txt,':\n',sep="")
    fprintf(fn,tmp_txt)
    for(ct1 in 1:globalStructure$m){
      fprintf(fn,'Group %g: ', ct1)
      for(ct2 in 1:globalStructure$na){
        tmp_txt <- '%g'
        if(opt_a) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
        if(ct2<globalStructure$na) tmp_txt <- paste(tmp_txt,' : ',sep="")
        fprintf(fn,tmp_txt,globalStructure$ga[ct1,ct2],globalStructure$gmina[ct1,ct2],globalStructure$gmaxa[ct1,ct2])
      }
      fprintf(fn,'\n')
    }
    fprintf(fn,'\n')
  }
  
  return( ) 
}

print_params <- function (params,name_str, fn, globalStructure, param_sqrt=FALSE,head_txt=NULL,matrix_elements=F,e_flag=FALSE) {
  if(is.null(head_txt)) head_txt <- "Parameter Values"
  uncer_txt <- ""
  if(e_flag) uncer_txt <- " (Uncertainty Distribution)"
  sqrt_txt <- ""
  if(param_sqrt) sqrt_txt <- " [sqrt(param)]"
  fprintf(fn,paste(head_txt,sqrt_txt,uncer_txt,":\n",sep=""))
  for(ct in 1:size(params,1)){
    
    par_val <- params[ct,2]
    
    uncer_val <- params[ct,3]
    
    cv_uncer <- sqrt(uncer_val)/par_val*100
    
    cv_str <- ", %CV="
    uncer_str <- ", Var="
    
    par_val_sqrt <- ""
    if(param_sqrt) par_val_sqrt =sqrt(par_val)
    #if(param_sqrt) cv_uncer =cv_uncer/2
    
    if((params[ct,1]==0)){ 
      dist_str <- "Point value" 
      uncer_val <- ""
      cv_uncer <- ""
      cv_str <- ""
      uncer_str <- ""
      
    }
    
    if((params[ct,1]==2)){
      dist_str <- "Uniform" 
      cv_uncer <- ""
      cv_str <- ""
      uncer_str <- ", Max-Min" 
    }
    
    if((params[ct,1]==4)){
      dist_str <- "Log-Normal" 
    }
    if(params[ct,1]==1){
      dist_str <- "Normal" 
    }
    
    if((params[ct,1]==3)){
      dist_str <- "User Defined"
      cv_uncer <- ""
      cv_str <- ""
    }
    
    if((params[ct,1]==5)){
      dist_str <- "Zero-Truncated Normal" 
    }
    
    if(!is.character(uncer_val)) uncer_val <- sprintf("%5.4g",uncer_val)
    if(!is.character(cv_uncer)) cv_uncer <- sprintf("%5.4g",cv_uncer)
    if(!is.character(par_val_sqrt)) par_val_sqrt <- sprintf("[%5.4g] ",par_val_sqrt)
    mat_str <- ""
    if(matrix_elements) mat_str <- sprintf(",%g",ct)

    
    if(e_flag){ 
      fprintf(fn,'%s[%g%s]: %5.4g %s(%s%s%s%s%s)\n', name_str,ct,mat_str,par_val, par_val_sqrt,dist_str, uncer_str, uncer_val, cv_str, cv_uncer)
    } else {
      fprintf(fn,'%s[%g%s]: %5.4g %s\n', name_str,ct,mat_str,par_val, par_val_sqrt,dist_str, uncer_str, uncer_val, cv_str, cv_uncer)
    }
  }
}
