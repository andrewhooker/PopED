#' Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).
#' 
#' Compute the expectation of the FIM given the model, parameters, distributions of parameter uncertainty, design and methods defined in the 
#' PopED database. Some of the arguments coming from the PopED database can be overwritten;  
#' by default these arguments are \code{NULL} in the 
#' function, if they are supplied then they are used instead of the arguments from the PopED database.
#' 
#' @inheritParams evaluate.fim
#' 
#' @return The FIM.
#' 
#' @family FIM
#' @family E-family
#'  
#' @examples 
#' \dontrun{
#' evaluate.ed.fim(poped.db)
#' }  
#' 


evaluate.e.ofv.fim <- function(poped.db,
                         calc.type=NULL,
                         bpop=poped.db$gbpop,
                         d=poped.db$gd,
                         covd=poped.db$covd,
                         docc=poped.db$docc,
                         sigma=poped.db$sigma,
                         model_switch=NULL,
                         ni=NULL,
                         xt=NULL,
                         x=NULL,
                         a=NULL,
                         groupsize=poped.db$design$groupsize,
                         deriv.type = NULL,
                         bLHS=poped.db$bLHS,
                         ofv_calc_type = poped.db$ofv_calc_type,
                         ED_samp_size = poped.db$ED_samp_size,
                         use.laplace=FALSE, # not working
                         ...){
  ## update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
    }
  }
  
  downsize.list <- downsizing_general_design(poped.db)
  if(is.null(ni)) ni <- downsize.list$ni
  if(is.null(xt)) xt <- downsize.list$xt
  if(is.null(model_switch)) model_switch <- downsize.list$model_switch
  if(is.null(x)) x <- downsize.list$x
  if(is.null(a)) a <- downsize.list$a    
 
  if(is.null(groupsize)) groupsize <- poped.db$downsized.design$groupsize
  
  if(!is.null(calc.type)) poped.db$iFIMCalculationType=calc.type
  
  if(!is.null(deriv.type)){ 
    poped.db$m1_switch=deriv.type
    poped.db$m2_switch=deriv.type
    poped.db$hle_switch=deriv.type
    poped.db$gradff_switch=deriv.type
    poped.db$gradfg_switch=deriv.type
  }
  
  E_fim <- NULL
  E_ofv <- NULL
  
  if(!use.laplace){
    output <- ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped.db)
    E_fim <- output$ED_fim
    E_ofv <- output$ED_ofv
    poped.db=output$globalStructure
  } else {  
    E_ofv  <- ed_laplace_ofv(c(),0, 0, model_switch,groupsize,ni,xt,x,a,bpop,d,covd,sigma,docc,poped.db)      
  }    
  return(list(E_ofv=E_ofv,E_fim= E_fim, poped.db=output$globalStructure))
}