#' Calculate the Fisher Information Matrix (FIM) and the OFV(FIM) for either point values or parameters or distributions.
#' 
#' This function computes the expectation of the FIM and OFV(FIM) for either point values of parameter estimates
#' or parameter distributions given the model, parameters, 
#' distributions of parameter uncertainty, design and methods defined in the 
#' PopED database.
#' 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param ofv The current ofv.  If other than zero then this value is simply returned unchanged.
#' @param fim The current FIM. If other than zero then this value is simply returned unchanged.
#' @param use_laplace Should the Laplace method be used in calculating the expectation of the OFV?  
#' @param laplace.fim Should an E(FIM) be calculated when computing the Laplace approximated E(OFV).  Typically
#' the FIM does not need to be computed and, if desired,  this calculation
#' is done using the standard MC integration technique, so can be slow. 
#' @param evaluate_fim Should the FIM be calculated?
#' 
#' @return A list containing the FIM and OFV(FIM) or the E(FIM) and E(OFV(FIM)) according  to the function arguments.
#' 
#' @family FIM
#' @family E-family
#' @family evaluate_FIM
#'  
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_ed.R
#' @example tests/testthat/examples_fcn_doc/examples_calc_ofv_and_fim.R
#' @export

calc_ofv_and_fim <- function (poped.db,
                              ofv=0,
                              fim=0, 
                              d_switch=poped.db$settings$d_switch,  
                              bpopdescr=poped.db$parameters$bpop, 
                              ddescr=poped.db$parameters$d,
                              bpop=bpopdescr[,2,drop=F], 
                              d=getfulld(ddescr[,2,drop=F],poped.db$parameters$covd), 
                              docc_full= getfulld(poped.db$parameters$docc[,2,drop=F],poped.db$parameters$covdocc), 
                              model_switch=poped.db$design$model_switch, 
                              ni=poped.db$design$ni, 
                              xt=poped.db$design$xt, 
                              x=poped.db$design$x, 
                              a=poped.db$design$a, 
                              fim.calc.type=poped.db$settings$iFIMCalculationType, 
                              use_laplace=poped.db$settings$iEDCalculationType, 
                              laplace.fim=FALSE, 
                              ofv_fun = poped.db$settings$ofv_fun,
                              evaluate_fim = TRUE,
                              ...) {
  
  ## compute the OFV
  if((ofv==0)){
    if(d_switch){ 
      if(is.null(ofv_fun)){
        if(!is.matrix(fim)){ 
          fmf <- evaluate.fim(poped.db,
                              bpop.val=bpop,
                              d_full=d,
                              docc_full=docc_full,
                              sigma_full=poped.db$parameters$sigma,
                              model_switch=model_switch,
                              ni=ni,
                              xt=xt,
                              x=x,
                              a=a,
                              groupsize=poped.db$design$groupsize,
                              fim.calc.type=fim.calc.type,
                              ...)
          
          #     returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          #     fmf <- returnArgs[[1]]
          #     poped.db <- returnArgs[[2]]
        }
        dmf=ofv_fim(fmf,poped.db,...)
      } else {
        ## update poped.db with options supplied in function
        called_args <- match.call()
        default_args <- formals()
        for(i in names(called_args)[-1]){
          if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
            #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
            if(!is.null(eval(parse(text=paste(i))))) eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
          }
        }
        out_tmp <- do.call(ofv_fun,list(poped.db,...))
        dmf <- out_tmp[[1]]
        fmf <- NULL
        if(length(out_tmp)>1) fmf <- out_tmp[[2]]
      }
    } else { # e-family
      if(is.null(ofv_fun)){
      output <- evaluate.e.ofv.fim(poped.db,
                                   fim.calc.type=fim.calc.type,
                                   bpop=bpopdescr,
                                   d=ddescr,
                                   covd=poped.db$parameters$covd,
                                   docc=poped.db$parameters$docc,
                                   sigma=poped.db$parameters$sigma,
                                   model_switch=model_switch,
                                   ni=ni,
                                   xt=xt,
                                   x=x,
                                   a=a,
                                   groupsize=poped.db$design$groupsize,
                                   use_laplace=use_laplace,
                                   laplace.fim=laplace.fim, 
                                   ...)
      dmf <- output$E_ofv
      fmf <- output$E_fim 
      } else {
        ## update poped.db with options supplied in function
        called_args <- match.call()
        default_args <- formals()
        for(i in names(called_args)[-1]){
          if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
            #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
            if(!is.null(eval(parse(text=paste(i))))) eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
          }
        }
        
        dmf <- mc_mean(ofv_fun,poped.db,...)
        fmf <- NULL
      }
    }
  
    ofv <- dmf
    if(!is.matrix(fim)) fim <- fmf
  }
  
  ## Should we compute the FIM?
  calc_fim <- FALSE
  if(is.null(fim)){
    calc_fim <- TRUE
  } else if(!is.matrix(fim)){ 
    calc_fim <- TRUE
  }    
  if(use_laplace && !laplace.fim) calc_fim <- FALSE
  if(!evaluate_fim) calc_fim <- FALSE
  
  if(calc_fim){
    if(d_switch){ 
      fmf <- evaluate.fim(poped.db,
                          bpop.val=bpop,
                          d_full=d,
                          docc_full=docc_full,
                          sigma_full=poped.db$parameters$sigma,
                          model_switch=model_switch,
                          ni=ni,
                          xt=xt,
                          x=x,
                          a=a,
                          groupsize=poped.db$design$groupsize,
                          fim.calc.type=fim.calc.type,
                          ...)
      
      #     returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
      #     fmf <- returnArgs[[1]]
      #     poped.db <- returnArgs[[2]]
    } else {
      output <- evaluate.e.ofv.fim(poped.db,
                                   fim.calc.type=fim.calc.type,
                                   bpop=bpopdescr,
                                   d=ddescr,
                                   covd=poped.db$parameters$covd,
                                   docc=poped.db$parameters$docc,
                                   sigma=poped.db$parameters$sigma,
                                   model_switch=model_switch,
                                   ni=ni,
                                   xt=xt,
                                   x=x,
                                   a=a,
                                   groupsize=poped.db$design$groupsize,
                                   use_laplace=FALSE,
                                   ...)  
      fmf<- output$E_fim
    }
    fim <- fmf
  }
  return(list(ofv=ofv,fim=fim)) 
}