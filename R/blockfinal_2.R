#' Result function for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the problem you solved.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams blockexp
#' @inheritParams blockheader_2
#' @param fmf_init Initial FIM.
#' @param dmf_init Initial OFV.
#' @param param_cvs_init The inital design parameter RSE values.
#' 
#' @family Helper
#' @export
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_blockfinal_2.R
# @keywords internal
#' 
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

blockfinal_2 <- function(fn,fmf,dmf,groupsize,ni,xt,x,a,model_switch,bpop,d,docc,sigma,m,poped.db,
                         opt_xt=poped.db$optsw[2],opt_a=poped.db$optsw[4],opt_x=poped.db$optsw[4],
                         fmf_init=NULL,dmf_init=NULL,param_cvs_init=NULL){
  
  fprintf(fn,'===============================================================================\nFINAL RESULTS\n\n')
  time_value = toc(echo=FALSE,name=".poped_total_time")
  if((opt_xt==TRUE)){
    print_xt(xt,ni,model_switch,fn,head_txt="Optimized Sampling Schedule\n")
    if(fn!="") print_xt(xt,ni,model_switch,head_txt="\nOptimized Sampling Schedule\n")
  }
  if((opt_x==TRUE)){
    #     fprintf(fn,'x :\n')
    #     fprintf(fn,'%g \n',x)
    #     cat("Optimized x values:\n")
    #     print(x)
    tmp_txt <- "\nOptimized Discrete Variables"
    tmp_txt <- paste(tmp_txt,':\n',sep="")
    fprintf(fn,tmp_txt)
    fprintf(tmp_txt)
    for(ct1 in 1:poped.db$m){
      fprintf(fn,'Group %g: ', ct1)
      fprintf('Group %g: ', ct1)
      for(ct2 in 1:poped.db$nx){
        tmp_txt <- '%g'
        if(ct2<poped.db$nx) tmp_txt <- paste(tmp_txt,' : ',sep="")        
        fprintf(fn,tmp_txt,x[ct1,ct2])
        fprintf(tmp_txt,x[ct1,ct2])
      }
      fprintf(fn,'\n')
      fprintf('\n')
    }
#     fprintf(fn,'\n')
#     fprintf('\n')
  }
  if((opt_a==TRUE)){
    tmp_txt <- "\nOptimized Covariates"
    tmp_txt <- paste(tmp_txt,':\n',sep="")
    fprintf(fn,tmp_txt)
    fprintf(tmp_txt)
    for(ct1 in 1:poped.db$m){
      fprintf(fn,'Group %g: ', ct1)
      fprintf('Group %g: ', ct1)
      for(ct2 in 1:poped.db$na){
        tmp_txt <- '%g'
        if(ct2<poped.db$na) tmp_txt <- paste(tmp_txt,' : ',sep="")
        fprintf(fn,tmp_txt,a[ct1,ct2])
        fprintf(tmp_txt,a[ct1,ct2])
      }
      fprintf(fn,'\n')
      fprintf('\n')
    }
#     fprintf(fn,'\n')
#     fprintf('\n')
    #     
    #     fprintf(fn,'a :\n')
    #     fprintf(fn,'%g \n',a)
    #     cat("Optimized a values:\n")
    #     print(a)
  }
  if((poped.db$d_switch==TRUE)){
    fprintf(fn,'\n FIM: \n')
    #write_matrix(fn,fmf)
    MASS::write.matrix(fmf,file=fn)
    fprintf(fn,'\n\nInverse(FIM):\n')
    #write_matrix(fn,inv(fmf))
    MASS::write.matrix(inv(fmf),file=fn)
  }
  fprintf(fn,'\ndet(FIM) = %g\n',dmf)
  
  param_vars=diag_matlab(inv(fmf))
  returnArgs <-  get_cv(param_vars,bpop,d,docc,sigma,poped.db) 
  params <- returnArgs[[1]]
  param_cvs <- returnArgs[[2]]

  fprintf(fn,'\nEfficiency criterion: det(FIM)^(1/npar) = %g\n',dmf^(1/length(params)))
  fprintf(fn,'\nEfficiency (final_design/initial_design): %g\n',(dmf^(1/length(params)))/(dmf_init^(1/length(params))))
  if(fn!="") fprintf('\nEfficiency (final_design/initial_design): %g\n',(dmf^(1/length(params)))/(dmf_init^(1/length(params))))
  
  
  parnam <- get_parnam(poped.db)
  fprintf(fn,'\nExpected parameter variance \nand relative standard error (%sRSE):\n','%')
  if(fn!="") fprintf('\nExpected parameter variance \nand relative standard error (%sRSE):\n','%')
  df <- data.frame("Parameter"=parnam,"Values"=params, "Variance"=param_vars, "RSE"=t(param_cvs*100),"RSE_initial_design"=t(param_cvs_init*100))
  print(df,digits=3, print.gap=3,row.names=F)
  if(fn!="") capture.output(print(df,digits=3, print.gap=3,row.names=F),file=fn)
  
  fprintf(fn,'\nTotal running time: %g seconds\n',time_value)
  if(fn!="") fprintf('\nTotal running time: %g seconds\n',time_value)
  
  return( ) 
}

print_xt <- function (xtopt, ni, model_switch,fn="",head_txt="Optimized xt values:\n",xt_other=NULL) {
  cat(head_txt,file=fn)
  for(j in 1:size(xtopt,1)){
    xtopt_i = xtopt[j,1:ni[j]]
    model_switch_i = model_switch[j,1:ni[j]]
    if(!is.null(xt_other)) xt_other_i = xt_other[j,1:ni[j]]
    for(i in unique(as.vector(model_switch_i))){
      xtopt_i_sort = sort(xtopt_i[model_switch_i==i])
      if(!is.null(xt_other)) xt_other_i_sort = xt_other_i[order(xtopt_i[model_switch_i==i])]
      if(length(unique(as.vector(model_switch_i)))>1) cat(sprintf("Model %g : ", i),file=fn)
      if(size(xtopt,1)>1) cat(sprintf("Group %g : ", j),file=fn)
      if(!is.null(xt_other)) {
        cat(sprintf("%6.4g", xt_other_i_sort),file=fn)
      } else {
        cat(sprintf("%6.4g", xtopt_i_sort),file=fn)
      }
      cat("\n",file=fn)
    }
  }
  invisible()
}

get_parnam <- function (poped.db) {
  nbpop = length(poped.db$notfixed_bpop)
  nd = length(poped.db$notfixed_d)
  ncovd = length(poped.db$notfixed_covd)
  ndocc = length(poped.db$notfixed_docc)
  ncovdocc = length(poped.db$notfixed_covdocc)
  nsigma = length(poped.db$notfixed_sigma)
  ncovsigma = length(poped.db$notfixed_covsigma)
  
  not_fixed <- list("bpop"=poped.db$notfixed_bpop,
                    "D"=poped.db$notfixed_d,
                    "D_cov"=poped.db$notfixed_covd,
                    "D.occ"=poped.db$notfixed_docc,
                    "D.occ_cov"=poped.db$notfixed_covdocc,
                    "SIGMA"=poped.db$notfixed_sigma,
                    "SIGMA_cov"=poped.db$notfixed_covsigma)
  parnam <- c()
  for(i in 1:size(not_fixed,2)){
    if(length(not_fixed[[i]])==0) next
    for(j in 1:length(not_fixed[[i]])){
      #     if(grep("_cov",names(not_fixed[i]))){
      #       k <- c()
      #       l <- 1
      #       while(length(k)<length(not_fixed[[i]])){
      #         k <- c(k,rep(l,l))
      #         l <- l+1
      #       }
      #     }
      if(not_fixed[[i]][j]==1){ 
        if(names(not_fixed[i])=="bpop") parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,"]",sep=""))  
        if(any(names(not_fixed[i])==c("D","SIGMA"))) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,",",j,"]",sep=""))
        if(length(grep("_cov",names(not_fixed[i])))!=0) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,"]",sep="")) 
      }
    }
  }
  return(parnam)
}
