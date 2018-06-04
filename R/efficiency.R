#' Compute efficiency.
#' 
#' Efficiency calculation between two designs.
#' 
#' 
#' @param ofv_init An initial objective function
#' @param ofv_final A final objective function.
#' @param npar The number of parameters to use for normalization.
#' @param poped_db a poped database
#' @inheritParams ofv_fim
#' @inheritParams poped_optim
#' @inheritParams create.poped.database
#' 
#' @return The specified efficiency value depending on the ofv_calc_type.  
#' The attribute "description" tells you how the calculation was made 
#' \code{attr(return_vale,"description")}
#' 
#' @family FIM
#' 
#' 
## @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
## @example tests/testthat/examples_fcn_doc/examples_ofv_criterion.R
#' 
#' @export

efficiency <- function(ofv_init, ofv_final, poped_db,
                       npar = get_fim_size(poped_db),
                       ofv_calc_type=poped_db$settings$ofv_calc_type,
                       ds_index = poped_db$parameters$ds_index) {
  
  
  eff = ofv_final/ofv_init
  attr(eff,"description") <- "ofv_final / ofv_init"
  
  if((ofv_calc_type==1) ){#D-Optimal Design
    eff = eff^(1/npar)
    attr(eff,"description") <- "(ofv_final / ofv_init)^(1/n_parameters)"
  }
  if((ofv_calc_type==4) ){#lnD-Optimal Design
    eff = (exp(ofv_final)/exp(ofv_init))^(1/npar)
    attr(eff,"description") <- "(exp(ofv_final) / exp(ofv_init))^(1/n_parameters)"
    
  }
  # if((ofv_calc_type==2) ){#A-Optimal Design
  #   eff=ofv_f/npar
  # }
  # 
  # if((ofv_calc_type==3) ){#S-Optimal Design
  #   stop(sprintf('Criterion for S-optimal design not implemented yet'))
  # }
  # 
  if((ofv_calc_type==6) ){#Ds-Optimal design
    eff = eff^(1/sum(!ds_index))
    attr(eff,"description") <- "(ofv_final / ofv_init)^(1/sum(interesting_parameters))"
    
  }   
  
  return( eff ) 
}

