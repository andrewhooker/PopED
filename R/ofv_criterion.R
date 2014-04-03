#' Normalize an objective function by the size of the FIM matrix
#' 
#' Compute a criterion of the FIM given the non-normalized FIM, model, parameters, design and methods defined in the 
#' PopED database. 
#' 
#' @param ofv_f An objective function
#' @param num_parameters The number of parameters to use for normalization
#' @param poped.db a poped database
#' 
#' @return The specified criterion value.
#' 
#' @family FIM
#' 
#' @examples 
#' \dontrun{
#' evaluate.fim(poped.db)
#' }  
#' 
ofv_criterion <- function(ofv_f,num_parameters,poped.db){
  
  #Input: the ofv
  #Return the single value that should be maximized
  
  criterion_value = 0
  
  if((poped.db$ofv_calc_type==1 || poped.db$ofv_calc_type==4) ){#D-Optimal Design
    criterion_value = ofv_f^(1/num_parameters)
    return
  }
  
  if((poped.db$ofv_calc_type==2) ){#A-Optimal Design
    criterion_value=ofv_f/num_parameters
  }
  
  if((poped.db$ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('Criterion for S-optimal design not implemented yet'))
  }
  
  if((poped.db$ofv_calc_type==6) ){#Ds-Optimal design
    criterion_value = ofv_f^(1/sum(poped.db$ds_index))
  }   
  
  return( criterion_value ) 
}
