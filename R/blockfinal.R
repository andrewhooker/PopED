#' Result function for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the problem you solved.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams blockexp
#' @inheritParams blockheader
#' @param fmf_init Initial FIM.
#' @param dmf_init Initial OFV.
#' @param param_cvs_init The initial design parameter RSE values in percent.
#' 
#' @family Helper
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_blockfinal.R
#' @export
#' @keywords internal
#' 
#' 
# @importFrom MASS write.matrix
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

blockfinal <- function(fn,fmf,dmf,groupsize,ni,xt,x,a,model_switch,bpop,d,docc,sigma,poped.db,
                       opt_xt=poped.db$settings$optsw[2],opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
                       opt_inds=poped.db$settings$optsw[5],
                       fmf_init=NULL,dmf_init=NULL,param_cvs_init=NULL,
                       compute_inv=TRUE,out_file=NULL,trflag=TRUE,footer_flag=TRUE,
                       run_time = NULL,
                       ...){
  time_value <- NULL
  
  if(!trflag) return(invisible() ) 
  if(footer_flag){
    if(is.null(fmf)) compute_inv <- FALSE
    if(!is.matrix(fmf)) compute_inv <- FALSE
    
    
    
    fprintf(fn,'===============================================================================\nFINAL RESULTS\n')
    if(fn!="") fprintf('===============================================================================\nFINAL RESULTS\n')
    
    time_value <- run_time
    if(is.null(time_value) & exists(".poped_total_time", envir=.PopedNamespaceEnv)) time_value = toc(echo=FALSE,name=".poped_total_time")
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
      if(fn!="") fprintf(tmp_txt)
      for(ct1 in 1:poped.db$design$m){
        fprintf(fn,'Group %g: ', ct1)
        if(fn!="") fprintf('Group %g: ', ct1)
        for(ct2 in 1:size(poped.db$design$x,2)){
          tmp_txt <- '%g'
          if(ct2<size(poped.db$design$x,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")        
          fprintf(fn,tmp_txt,x[ct1,ct2])
          if(fn!="") fprintf(tmp_txt,x[ct1,ct2])
        }
        fprintf(fn,'\n')
        if(fn!="") fprintf('\n')
      }
      #     fprintf(fn,'\n')
      #     fprintf('\n')
    }
    if((opt_a==TRUE)){
      tmp_txt <- "\nOptimized Covariates"
      tmp_txt <- paste(tmp_txt,':\n',sep="")
      fprintf(fn,tmp_txt)
      if(fn!="") fprintf(tmp_txt)
      for(ct1 in 1:poped.db$design$m){
        fprintf(fn,'Group %g: ', ct1)
        if(fn!="") fprintf('Group %g: ', ct1)
        for(ct2 in 1:size(poped.db$design$a,2)){
          tmp_txt <- '%g'
          if(ct2<size(poped.db$design$a,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")
          fprintf(fn,tmp_txt,a[ct1,ct2])
          if(fn!="") fprintf(tmp_txt,a[ct1,ct2])
        }
        fprintf(fn,'\n')
        if(fn!="") fprintf('\n')
      }
      #     fprintf(fn,'\n')
      #     fprintf('\n')
      #     
      #     fprintf(fn,'a :\n')
      #     fprintf(fn,'%g \n',a)
      #     cat("Optimized a values:\n")
      #     print(a)
    }
    if((opt_inds==TRUE)){
      tmp_txt <- "\nOptimized groupsize"
      tmp_txt <- paste(tmp_txt,':\n',sep="")
      fprintf(fn,tmp_txt)
      if(fn!="") fprintf(tmp_txt)
      for(ct1 in 1:poped.db$design$m){
        fprintf(fn,'Group %g: ', ct1)
        if(fn!="") fprintf('Group %g: ', ct1)
        for(ct2 in 1:size(poped.db$design$groupsize,2)){
          tmp_txt <- '%g'
          if(ct2<size(poped.db$design$groupsize,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")
          fprintf(fn,tmp_txt,groupsize[ct1,ct2])
          if(fn!="") fprintf(tmp_txt,groupsize[ct1,ct2])
        }
        fprintf(fn,'\n')
        if(fn!="") fprintf('\n')
      }
      #     fprintf(fn,'\n')
      #     fprintf('\n')
      #     
      #     fprintf(fn,'a :\n')
      #     fprintf(fn,'%g \n',a)
      #     cat("Optimized a values:\n")
      #     print(a)
    }
    if((poped.db$settings$d_switch==TRUE && (fn!="" || trflag>1))){
      fprintf(fn,'\n FIM: \n')
      #write_matrix(fn,fmf)
      MASS::write.matrix(fmf,file=fn)
      fprintf(fn,'\n\nInverse(FIM):\n')
      #write_matrix(fn,inv(fmf))
      if(compute_inv) MASS::write.matrix(inv(fmf),file=fn)
    }
    fprintf(fn,'\nOFV = %g\n',dmf)
    if(fn!="") fprintf('\nOFV = %g\n',dmf)
    
    if(compute_inv){
      #param_vars=diag_matlab(inv(fmf))
      #returnArgs <-  get_cv(param_vars,bpop,d,docc,sigma,poped.db) 
      #params <- returnArgs[[1]]
      #param_cvs <- returnArgs[[2]]
      param_cvs <- get_rse(fim=fmf,poped.db,bpop,diag(d),docc,sigma)
    }
    
    output <- get_unfixed_params(poped.db)
    params <- output$all
    npar <- length(params)
    
    if(fn!="" || trflag>1) fprintf(fn,'\nEfficiency criterion [usually defined as det(FIM)^(1/npar)]  = %g\n',
            ofv_criterion(dmf,npar,poped.db))
    
    # fprintf(fn,'\nEfficiency [typically: (OFV_final/OFV_initial)^(1/npar)]: %g\n',
    #         ofv_criterion(dmf,npar,poped.db)/ofv_criterion(dmf_init,npar,poped.db))
    # if(fn!=""){
    #   fprintf('\nEfficiency [typically: (OFV_final/OFV_initial)^(1/npar)]: %g\n',
    #           ofv_criterion(dmf,npar,poped.db)/ofv_criterion(dmf_init,npar,poped.db))
    # }
    
    if(!is.null(dmf_init)){
      eff <- efficiency(dmf_init, dmf, poped.db,...)
      fprintf(fn,"\nEfficiency: \n  (%s) = %.5g\n",attr(eff,"description"),eff,both=TRUE)
    }
    # fprintf(fn,'\nEfficiency (Final/Initial): %0.5g\n',
    #         ofv_criterion(dmf,npar,poped.db)/ofv_criterion(dmf_init,npar,poped.db),both=TRUE)
    # 
    #fprintf(fn,'\nEfficiency criterion: det(FIM)^(1/npar) = %g\n',dmf^(1/length(params)))
    #fprintf(fn,'\nEfficiency (final_design/initial_design): %g\n',(dmf^(1/length(params)))/(dmf_init^(1/length(params))))
    #if(fn!="") fprintf('\nEfficiency (final_design/initial_design): %g\n',(dmf^(1/length(params)))/(dmf_init^(1/length(params))))
    
    # if((poped.db$settings$ofv_calc_type==4) && !is.null(fmf_init) && !is.null(fmf) ){#D-Optimal Design
    #   fprintf(fn,'\nD-Efficiency [(det(FIM_final)/det(FIM_initial))^(1/npar)]: %g\n',
    #           (det(fmf)/det(fmf_init))^(1/npar))
    #   if(fn!=""){
    #     fprintf(fn,'\nD-Efficiency [(det(FIM_final)/det(FIM_initial))^(1/npar)]: %g\n',
    #             (det(fmf)/det(fmf_init))^(1/npar))
    #   }
    # }
    # 
    
    if(is.null(param_cvs_init) && !is.null(fmf_init) && is.matrix(fmf_init) && compute_inv){
      if(is.finite(dmf_init)) {
        #param_vars_init=diag_matlab(inv(fmf_init))
        #returnArgs <-  get_cv(param_vars_init,bpop,d,docc,sigma,poped.db) 
        #params_init <- returnArgs[[1]]
        #param_cvs_init <- returnArgs[[2]]
        param_cvs_init <- get_rse(fim=fmf_init,poped.db,bpop,diag(d),docc,sigma)
      } else {
        param_cvs_init <- suppressMessages(suppressWarnings(get_rse(fim=fmf_init,poped.db,bpop,diag(d),docc,sigma)))
      }
    }
    
    if(compute_inv){
      parnam <- get_parnam(poped.db)
      fprintf(fn,'\nExpected relative standard error\n(%sRSE, rounded to nearest integer):\n','%')
      if(fn!="") fprintf('\nExpected relative standard error\n(%sRSE, rounded to nearest integer):\n','%')
      df <- data.frame("Parameter"=parnam,"Values"=sprintf("%6.3g",params), #"Variance"=param_vars, 
                       "RSE_0"=round(param_cvs_init),"RSE"=round(param_cvs))
                       #"RSE_0"=sprintf("%6.3g",param_cvs_init),"RSE"=sprintf("%6.3g",param_cvs))
      print(df,digits=3, print.gap=3,row.names=F)
      if(fn!="") capture.output(print(df,digits=3, print.gap=3,row.names=F),file=fn)
    }
    
    if(!is.null(time_value)){
      fprintf(fn,'\nTotal running time: %g seconds\n',time_value)
      if(fn!="") fprintf('\nTotal running time: %g seconds\n',time_value)
    }
  } # end footer_flag
  if(!any(class(out_file)=="file") && (fn != '')) close(fn)
  return(invisible(time_value)) 
}

print_xt <- function (xtopt, ni, model_switch,fn="",head_txt="Optimized sample times:\n",xt_other=NULL,
                      digits=4) {
  cat(head_txt,file=fn)
  for(j in 1:size(xtopt,1)){
    xtopt_i = xtopt[j,1:ni[j]]
    model_switch_i = model_switch[j,1:ni[j]]
    if(!is.null(xt_other)) xt_other_i = xt_other[j,1:ni[j]]
    for(i in unique(as.vector(model_switch_i))){
      xtopt_i_sort = sort(xtopt_i[model_switch_i==i])
      if(!is.null(xt_other)) xt_other_i_sort = xt_other_i[order(xtopt_i[model_switch_i==i])]
      # if(size(xtopt,1)>1) cat(sprintf("Group %g : ", j),file=fn)
      cat(sprintf("Group %g: ", j),file=fn)
      if(length(unique(as.vector(model_switch_i)))>1) cat(sprintf("Model %g: ", i),file=fn)
      if(!is.null(xt_other)) {
        cat(sprintf(paste0("%",digits+2,".",digits,"g"), xt_other_i_sort),file=fn)
      } else {
        cat(sprintf(paste0("%",digits+2,".",digits,"g"), xtopt_i_sort),file=fn)
        # cat(sprintf("%6.4g", xtopt_i_sort),file=fn)
      }
      cat("\n",file=fn)
    }
  }
  invisible()
}

# print_a <- function (aopt,fn="",head_txt="Optimized covariates:\n") {
#   cat(head_txt,file=fn)
#   for(j in 1:size(a,1)){
#     aopt_i = aopt[j,1:ni[j]]
#     model_switch_i = model_switch[j,1:ni[j]]
#     if(!is.null(xt_other)) xt_other_i = xt_other[j,1:ni[j]]
#     for(i in unique(as.vector(model_switch_i))){
#       xtopt_i_sort = sort(xtopt_i[model_switch_i==i])
#       if(!is.null(xt_other)) xt_other_i_sort = xt_other_i[order(xtopt_i[model_switch_i==i])]
#       if(size(xtopt,1)>1) cat(sprintf("Group %g : ", j),file=fn)
#       if(length(unique(as.vector(model_switch_i)))>1) cat(sprintf("Model %g : ", i),file=fn)
#       if(!is.null(xt_other)) {
#         cat(sprintf("%6.4g", xt_other_i_sort),file=fn)
#       } else {
#         cat(sprintf("%6.4g", xtopt_i_sort),file=fn)
#       }
#       cat("\n",file=fn)
#     }
#   }
#   invisible()
# }

get_parnam <- function (poped.db) {
  nbpop = length(poped.db$parameters$notfixed_bpop)
  nd = length(poped.db$parameters$notfixed_d)
  ncovd = length(poped.db$parameters$notfixed_covd)
  ndocc = length(poped.db$parameters$notfixed_docc)
  ncovdocc = length(poped.db$parameters$notfixed_covdocc)
  nsigma = length(poped.db$parameters$notfixed_sigma)
  ncovsigma = length(poped.db$parameters$notfixed_covsigma)
  
  not_fixed <- list("bpop"=poped.db$parameters$notfixed_bpop,
                    "D"=poped.db$parameters$notfixed_d,
                    "D_cov"=poped.db$parameters$notfixed_covd,
                    "D.occ"=poped.db$parameters$notfixed_docc,
                    "D.occ_cov"=poped.db$parameters$notfixed_covdocc,
                    "SIGMA"=poped.db$parameters$notfixed_sigma,
                    "SIGMA_cov"=poped.db$parameters$notfixed_covsigma)
  
  bpop_names <- rownames(poped.db$parameters$bpop)
  d_names <- rownames(poped.db$parameters$d)
  sig_names <- rownames(poped.db$parameters$sigma)
  
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
        if(names(not_fixed[i])=="bpop"){
          default_name <- TRUE
          if(!is.null(bpop_names)){
            if(bpop_names[j]!="") {
              default_name <- FALSE
              parnam <- c(parnam,bpop_names[j])
            }
          } 
          if(default_name)  parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,"]",sep=""))    
        } 
        if(any(names(not_fixed[i])==c("D"))){
          default_name <- TRUE
          if(!is.null(d_names)){
            if(d_names[j]!="") {
              default_name <- FALSE
              parnam <- c(parnam,paste0("d_",d_names[j]))
            }
          } 
          if(default_name) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,",",j,"]",sep=""))
        }
        if(any(names(not_fixed[i])==c("SIGMA"))){
          default_name <- TRUE
          if(!is.null(sig_names)){
            if(sig_names[j]!="") {
              default_name <- FALSE
              parnam <- c(parnam,paste0("sig_",sig_names[j]))
            }
          } 
          if(default_name) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,",",j,"]",sep=""))
        }
        if(any(names(not_fixed[i])==c("D_cov"))){
          mat_ind <- which(lower.tri(poped.db$parameters$param.pt.val$d, diag = FALSE), arr.ind=T)[j,]
          parnam <- c(parnam,paste("D","[",mat_ind[1],",",mat_ind[2],"]",sep=""))
        }
        
        if(any(names(not_fixed[i])==c("D.occ"))) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,",",j,"]",sep=""))
        if((length(grep("_cov",names(not_fixed[i])))!=0) && (!any(names(not_fixed[i])==c("D_cov")))) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,"]",sep="")) 
      }
    }
  }
 
  
  return(parnam)
}
