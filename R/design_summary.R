#' Display a summary of output from poped_db
#'
#' @param poped_db An object returned from \code{\link{create.poped.database}} to summarize. 
#' @param file A file handle to write to.  Default is to the R console.
#' @param ... Additional arguments. Passed to \code{\link{blockfinal}}.
#'
#' @return NULL
# @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
# @example tests/testthat/examples_fcn_doc/examples_summary.poped_optim.R
#' @export

design_summary <- function(poped_db,
                           file="",
                           ...){
  
  
  print_xt(poped_db$design$xt,poped_db$design$ni,poped_db$design$model_switch,file,
           head_txt="Sampling Schedule\n------------\n",...)#digits=digits)
  
  if(ncol(poped_db$design$a)>0){
    cat("\n\nCovariates\n------------\n",file=file)
    tmp_mat <- poped_db$design$a
    names_tmp <- gsub("grp_","Group ",rownames(tmp_mat))
    names_tmp <- paste0(names_tmp,":")
    rownames(tmp_mat) <- names_tmp                                    
    print(tmp_mat,
          ...)#digits=digits)
  }
  
  if(ncol(poped_db$design$x)>0){
    cat("\n\nDiscrete design variables\n------------\n",file=file)
    tmp_mat <- poped_db$design$x
    names_tmp <- gsub("grp_","Group ",rownames(tmp_mat))
    names_tmp <- paste0(names_tmp,":")
    rownames(tmp_mat) <- names_tmp                                    
    print(tmp_mat,
          ...)#digits=digits)  
  }
 
  
  if(ncol(poped_db$design$groupsize)>0){
    cat("\n\nGroup Size\n------------\n",file=file)
    tmp_mat <- poped_db$design$groupsize
    names_tmp <- gsub("grp_","Group ",rownames(tmp_mat))
    names_tmp <- paste0(names_tmp,":")
    rownames(tmp_mat) <- names_tmp    
    colnames(tmp_mat) <- "N" 
    print(tmp_mat,
          ...)#digits=digits)  
  }
  
  return(invisible())
  
}