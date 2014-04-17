#' Evaluate a criterion of the Fisher Information Matrix (FIM)
#' 
#' Compute a criterion of the FIM given the model, parameters, design and methods defined in the 
#' PopED database. 
#' 
#' @param fmf The FIM
#' @param poped.db A poped database
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' 
#' @return The specified criterion value.
#' 
#' @family FIM
#' @family evaluate_FIM
#' 
#' @examples 
#' \dontrun{
#' evaluate.fim(poped.db)
#' }  
#' ## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

ofv_fim <- function(fmf,poped.db,
                    ofv_calc_type = poped.db$ofv_calc_type,
                    ds_index=poped.db$ds_index){
  
  #Input: the FIM
  #Return the single value that should be maximized
  
  ofv_value = 0
  
  if((!isempty(poped.db$prior_fim) && all(size(poped.db$prior_fim)==size(fmf)))){
    fmf = fmf + poped.db$prior_fim
  }
  
  if((ofv_calc_type==1) ){#determinant of FIM
    if((isempty(fmf))){
      ofv_value = 0
    } else {
      ofv_value = det(fmf)
    }
    return
  }
  
  if((ofv_calc_type==2) ){#trace of the inverse of FIM
    imf = inv(fmf)
    ofv_value = trace_matrix(imf)
    ofv_value = 1/ofv_value #Make it a max-problem
    return
  }
  
  if((ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('S-optimal design not implemented yet!'))
  }
  
  if((ofv_calc_type==4) ){#log determinant of FIM
    #ofv_value = sum(log(svd(fmf)))
    ofv_value = log(det(fmf))
  }
  
  if((ofv_calc_type==5) ){# C-optimal design
    stop(sprintf('C-optimal design is not implemented yet!'))
  }
  if((ofv_calc_type==6) ){#Ds-optimal design
    tmp = fmf
    tmp <- tmp[c(col(ds_index)[ds_index==1]),,drop=F]
    tmp <- tmp[,c(col(ds_index)[ds_index==1]),drop=F]
    ofv_value = det(fmf)/det(tmp)
  }
  if((ofv_calc_type==7) ){#sum of CV
    if((sum(sum(fmf))!=0 && !isnan(sum(sum(fmf))))){
      imf = inv(fmf)
      returnArgs <-  get_cv(diag_matlab(imf),poped.db$gbpop,poped.db$gd,poped.db$docc,poped.db$sigma,poped.db) 
      params <- returnArgs[[1]]
      params_cvs <- returnArgs[[2]]
      if((isnan(sum(diag_matlab(imf))))){
        ofv_value = 0
      } else {
        ofv_value=1/sum(abs(params_cvs))
      }
    } else {
      ofv_value = 0
    }
    return
  }
  return( ofv_value ) 
}
