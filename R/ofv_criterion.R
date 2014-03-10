ofv_criterion <- function(ofv_f,num_parameters,globalStructure){
  
  #Input: the ofv
  #Return the single value that should be maximized
  
  criterion_value = 0
  
  if((globalStructure$ofv_calc_type==1 || globalStructure$ofv_calc_type==4) ){#D-Optimal Design
    criterion_value = ofv_f^(1/num_parameters)
    return
  }
  
  if((globalStructure$ofv_calc_type==2) ){#A-Optimal Design
    criterion_value=ofv_f/num_parameters
  }
  
  if((globalStructure$ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('Criterion for S-optimal design not implemented yet'))
  }
  
  if((globalStructure$ofv_calc_type==6) ){#Ds-Optimal design
    criterion_value = ofv_f^(1/sum(globalStructure$ds_index))
  }   
  
  return( criterion_value ) 
}
