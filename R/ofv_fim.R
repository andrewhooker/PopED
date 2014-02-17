## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

ofv_fim <- function(fmf,globalStructure,
                    ofv_calc_type = globalStructure$ofv_calc_type,
                    ds_index=globalStructure$ds_index){
  
  #Input: the FIM
  #Return the single value that should be maximized
  
  ofv_value = 0
  
  if((!isempty(globalStructure$prior_fim) && all(size(globalStructure$prior_fim)==size(fmf)))){
    fmf = fmf + globalStructure$prior_fim
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
    ofv_value = sum(log(svd(fmf)))
    #ofv_value = log(det(fmf))
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
      returnArgs <-  get_cv(diag_matlab(imf),globalStructure$gbpop,globalStructure$gd,globalStructure$docc,globalStructure$sigma,globalStructure) 
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
