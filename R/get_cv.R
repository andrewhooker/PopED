get_cv <- function(param_vars,poped.db){
  #Return the RSE,CV of parameters
  ## Author: Andrew Hooker
  params_all <-  get_all_params(poped.db)[[8]] 
 
  returnArgs <-  get_unfixed_params(poped.db,params_all) 
  params <- returnArgs[[8]]
  var_derivative <- returnArgs[[9]]
  
  if((length(param_vars)!=length(params))){
    stop("Number of unfixed parameters not the same as the size of the FIM,\nno RSE can be computed!\n")
  }
  
  params_cv = zeros(size(param_vars))
  for(i in 1:length(params)){
    if((params[i]!=0)){
      if((var_derivative[i]==1)){
        params_cv[i] = sqrt(param_vars[i])/abs(params[i])
      } else { #Derivative w$r.t to SD instead of var
        params_cv[i] = sqrt(param_vars[i])/sqrt(params[i])
      }
    } else {
      params_cv[i] = sqrt(param_vars[i])
    }
  }
  
  return(list( params= params, params_cv = params_cv )) 
}

#' Compute the expected parameter relative standard errors 
#' 
#' This function  computes the expected relative standard errors of a model given a design and a previously computed
#' FIM.
#' 
#' @param fim A Fisher Information Matrix (FIM).
#' @param bpop A vector containing the values of the fixed effects used to compute the \code{fim}. 
#' @param d A vector containing the values of the diagonals of the between subject variability matrix.
#' @param use_percent Should RSE be reported as percent? 
#' @param prior_fim A prior FIM to be added to the \code{fim}. Should be the same size as the \code{fim}.
#' @param ... Additional arguments passed to \code{\link{inv}}. 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' 
#' @return A named list of RSE values.  If the estimated parameter is assumed to be zero then for that 
#'   parameter the standard error is returned.
#' 
#' @family evaluate_design
#' 
# @example inst/examples_fcn_doc/examples_evaluate.fim.R
#' @example tests/testthat/examples_fcn_doc/examples_evaluate.fim.R
#' @export
get_rse <- function (fim, poped.db,
                     bpop=poped.db$parameters$bpop[,2],
                     #bpop=poped.db$parameters$bpop[,2,drop=F],
                     d=poped.db$parameters$d[,2],
                     # d=poped.db$parameters$d[,2,drop=F],
                     docc=poped.db$parameters$docc,
                     sigma=poped.db$parameters$sigma,
                     use_percent=TRUE,
                     fim.calc.type=poped.db$settings$iFIMCalculationType,
                     prior_fim = poped.db$settings$prior_fim,
                     #pseudo_on_fail = FALSE,
                     ...) {
  
  ## update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      # if (i %in% c('bpop','d')) {
      #   if (eval(parse(text=paste("dim(",i,")[2]>1"))))
      #     (eval(parse(text=paste(i, "<-",i,"[,2]"))))
      # }
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
    }
  }

  ## if prior is given in poped.db then add it to the given fim
  if((!isempty(prior_fim) && all(size(prior_fim)==size(fim)))){
    fim = fim + prior_fim
  }
  
  inv_fim <- tryCatch({
    inv(fim,...)
  }, error=function(e){
    warning(e)
    return(NULL)
  })
   
  if(is.null(inv_fim)){
    mess <- paste0("\n  Could not invert the FIM.",
                   "\n  Is the design adequate to estimate all parameters?")
    eig <- eigen(fim)[["values"]]
    names(eig) <- get_parnam(poped.db)
    neg.vals <- eig[eig< 0]
    num.neg <- length(neg.vals)
    if(num.neg>0){
      mess <- paste0(mess,"\n  Potentially problematic parameters and associated eigenvalues:")
      for(i in 1:num.neg){
        mess <- paste0(mess,sprintf("\n %12s  %8.7e",names(neg.vals[i]),neg.vals[i]))
      }
    }
    #warning(simpleWarning(mess,call="get_rse()"))
    warning(mess)
    return(rep(NA,length(get_parnam(poped.db))))
  }  
  param_vars=diag_matlab(inv_fim)
  returnArgs <-  get_cv(param_vars,poped.db) 
  params <- returnArgs[[1]]
  params_rse <- returnArgs[[2]]
  parnam <- get_parnam(poped.db)
  ret <- params_rse[,,drop=T]
  if(use_percent) ret[params!=0]=ret[params!=0]*100
  names(ret) <- parnam
  if(any(ret==0)){
    zero_ret <- names(ret[ret==0])
    mess <- paste0("  The following parameters are not estimable:\n  ",
                   paste0(zero_ret,collapse = ", "),
                   "\n  Is the design adequate to estimate all parameters?")
    warning(mess, call. = FALSE)
    ret[ret==0] <- NA
  }
  return(ret)
}
