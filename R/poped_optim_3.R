#' Optimization main module for PopED
#' 
#' Optimize the objective function. The function works for both discrete and 
#' continuous optimization variables. If more than one optimization method is 
#' specified then the methods are run in series.  If \code{loop_methods=TRUE} 
#' then the series of optimization methods will be run for \code{iter_max} 
#' iterations, or until the efficiency of the design after the current series 
#' (compared to the start of the series) is less than \code{stop_crit_eff}.
#' 
#' This function takes information from the PopED database supplied as an 
#' argument. The PopED database supplies information about the the model, 
#' parameters, design and methods to use. Some of the arguments coming from the 
#' PopED database can be overwritten; if they are supplied then they are used 
#' instead of the arguments from the PopED database.
#' 
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams Dtrace
#' @inheritParams calc_ofv_and_fim
#' @inheritParams optim_LS
#' @param ... arguments passed to other functions.
#' @param control Contains control arguments specified for each method separately. 
#' @param method A vector of optimization methods to use in a sequential 
#'   fashion.  Options are \code{c("ARS","BFGS","LS","GA")}. \code{c("ARS")} is 
#'   for Adaptive Random Search \code{\link{optim_ARS}}.  \code{c("LS")} is for 
#'   Line Search \code{\link{optim_LS}}. \code{c("BFGS")} is for Method 
#'   "L-BFGS-B" from \code{\link[stats]{optim}}. \code{c("GA")} is for the 
#'   genetic algorithm from \code{\link[GA]{ga}}. If \code{opt_inds=TRUE} then
#'   this optimization is always added to the end of the sequential optimization.
#' @param out_file Save output from the optimization to a file.
#' @param loop_methods Should the optimization methods be looped for
#'   \code{iter_max} iterations, or until the efficiency of the design after the
#'   current series (compared to the start of the series) is less than, or equal to,
#'   \code{stop_crit_eff}?
#' @param stop_crit_eff If \code{loop_methods==TRUE}, the looping will stop if the
#'   efficiency of the design after the current series (compared to the start of
#'   the series) is less than, or equal to, \code{stop_crit_eff} (if \code{maximize==FALSE} then 1/stop_crit_eff is the cut
#'   off and the efficiency must be greater than or equal to this value to stop the looping).
#' @param stop_crit_diff If \code{loop_methods==TRUE}, the looping will stop if the
#'   difference in criterion value of the design after the current series (compared to the start of
#'   the series) is less than, or equal to, \code{stop_crit_diff} (if \code{maximize==FALSE} then -stop_crit_diff is the cut
#'   off and the difference in criterion value must be greater than or equal to this value to stop the looping).
#' @param stop_crit_rel If \code{loop_methods==TRUE}, the looping will stop if the
#'   relative difference in criterion value of the design after the current series (compared to the start of
#'   the series) is less than, or equal to, \code{stop_crit_rel} (if \code{maximize==FALSE} then -stop_crit_rel is the cut
#'   off and the relative difference in criterion value must be greater than or equal to this value to stop the looping).
#' @param maximize Should the objective function be maximized or minimized?
#'   
#'   
#'   
#' @references \enumerate{ \item M. Foracchia, A.C. Hooker, P. Vicini and A. 
#'   Ruggeri, "PopED, a software fir optimal experimental design in population 
#'   kinetics", Computer Methods and Programs in Biomedicine, 74, 2004. \item J.
#'   Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. 
#'   Hooker, "PopED: An extended, parallelized, nonlinear mixed effects models 
#'   optimal design tool", Computer Methods and Programs in Biomedicine, 108, 
#'   2012. }
#'   
#' @family Optimize
#'   
#' @keywords internal
# @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
# @example tests/testthat/examples_fcn_doc/examples_poped_optim.R
# @export

poped_optim_3 <- function(poped.db,
                          opt_xt=poped.db$settings$optsw[2],
                          opt_a=poped.db$settings$optsw[4],
                          opt_x=poped.db$settings$optsw[3],
                          opt_samps=poped.db$settings$optsw[1],
                          opt_inds=poped.db$settings$optsw[5],
                          method=c("ARS","BFGS","LS"),
                          control=list(),
                          trace = TRUE,
                          fim.calc.type=poped.db$settings$iFIMCalculationType,
                          ofv_calc_type=poped.db$settings$ofv_calc_type,
                          ds_index=poped.db$parameters$ds_index,
                          approx_type=poped.db$settings$iApproximationMethod,
                          d_switch=poped.db$settings$d_switch,
                          ED_samp_size = poped.db$settings$ED_samp_size,
                          bLHS=poped.db$settings$bLHS,
                          use_laplace=poped.db$settings$iEDCalculationType,
                          out_file="",
                          parallel=F,
                          parallel_type=NULL,
                          num_cores = NULL,
                          loop_methods=ifelse(length(method)>1,TRUE,FALSE),
                          iter_max = 10,
                          stop_crit_eff = 1.001,
                          stop_crit_diff = NULL,
                          stop_crit_rel = NULL,
                          ofv_fun = poped.db$settings$ofv_fun,
                          maximize=T,
                          allow_replicates=TRUE,
                          allow_replicates_xt=TRUE,
                          allow_replicates_a=TRUE,
                          ...){
  
  #------------ update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
    }
  }
  

  
  #------------- initialization
  fmf = 0 #The best FIM so far
  dmf = 0 #The best ofv of FIM  so far
  #output <-calc_ofv_and_fim(poped.db,...)
  
  #--------------------- Header information
  ## only for header information
  if(is.null(ofv_fun) || is.function(ofv_fun)){
    ofv_fun_user <- ofv_fun 
  } else {
    # source explicit file
    # here I assume that function in file has same name as filename minus .txt and pathnames
    if(file.exists(as.character(ofv_fun))){
      source(as.character(ofv_fun))
      ofv_fun_user <- eval(parse(text=fileparts(ofv_fun)[["filename"]]))
    } else {
      stop("ofv_fun is not a function or NULL, and no file with that name was found")
    }
    
  }
  if(!is.null(ofv_fun)){
    poped.db$settings$ofv_calc_type = 0
    ofv_calc_type = 0
  }
  
  output <-calc_ofv_and_fim(poped.db,d_switch=d_switch,
                            ED_samp_size=ED_samp_size,
                            bLHS=bLHS,
                            use_laplace=use_laplace,
                            ofv_calc_type=ofv_calc_type,
                            fim.calc.type=fim.calc.type,
                            ofv_fun = ofv_fun_user,
                            ...)
  
  fmf <- output$fim
  dmf <- output$ofv
  fmf_init <- fmf
  dmf_init <- dmf
  poped.db_init <- poped.db
  
  if(is.nan(dmf_init)) stop("Objective function of initial design is NaN")
  
  fn=blockheader(poped.db,name="optim",e_flag=!d_switch,
                 fmf=fmf_init,dmf=dmf_init,
                 out_file=out_file,
                 trflag=trace,
                 ...)
  
  # handle optimization of n distribution
  if(opt_inds){
    if(poped.db$design$m==1){
      message("There is only one group, so proportions of individuals in each group cannot be optimized\n")
      opt_inds <- FALSE
    } else{
      method <- c(method,"opt_n_dist")
      if(missing(loop_methods) & (length(method)>1)) loop_methods <- TRUE
      opt_inds <- FALSE
    }
  } 
  
  #------------ optimize
  if(!(fn=="")) sink(fn, append=TRUE, split=TRUE)
  
  iter <- 0
  stop_crit <- FALSE
  output$ofv <- dmf_init
  while(stop_crit==FALSE && iter < iter_max){
    ofv_init <- output$ofv
    iter=iter+1
    method_loop <- method
    if(loop_methods){
      cat("************* Iteration",iter," for all optimization methods***********************\n\n") 
    }
    
    while(length(method_loop)>0){
      cur_meth <- method_loop[1]
      method_loop <- method_loop[-1]
      if(!(cur_meth %in% c("ARS","LS","GA","BFGS","opt_n_dist"))) 
        stop("Current method ", cur_meth, " not defined")
      if(cur_meth=="ARS"){
        cat("*******************************************\n")
        cat("Running Adaptive Random Search Optimization\n")
        cat("*******************************************\n")
        
        # Get parameters and space
        ps_tbl <- get_par_and_space_optim(poped.db,
                                          opt_xt=opt_xt,
                                          opt_a=opt_a,
                                          opt_x=opt_x,
                                          opt_samps=opt_samps,
                                          opt_inds=opt_inds,
                                          transform_parameters=F)
        
        # handle control arguments
        con <- list(trace = trace, 
                    parallel=parallel,
                    parallel_type=parallel_type,
                    num_cores = num_cores)
        con[names(control$ARS)] <- control$ARS

        # handle replicates
        if(!allow_replicates)
          con$replicates_index=ps_tbl$type
        if(!allow_replicates_xt | !allow_replicates_a)            
          con$replicates_index <- seq(1,length(ps_tbl$par))
        if(!allow_replicates_xt)
          con$replicates_index[ps_tbl$type=="xt"]="xt"
        if(!allow_replicates_a)
          con$replicates_index[ps_tbl$type=="a"]="a"
              
        tmp_ofv_fun <- function(par,...){
          ofv_optim(par,ps_tbl,poped.db,
                    d_switch=d_switch,
                    ED_samp_size=ED_samp_size,
                    bLHS=bLHS,
                    use_laplace=use_laplace,
                    ofv_calc_type=ofv_calc_type,
                    fim.calc.type=fim.calc.type,
                    ofv_fun = ofv_fun_user,
                    ...)}
        
        output <- 
          do.call(optim_ARS,
                  c(list(par=ps_tbl$par,
                         fn=tmp_ofv_fun,
                         lower=ps_tbl$lower,
                         upper=ps_tbl$upper,
                         allowed_values = ps_tbl$allowed_values,
                         maximize=maximize),
                    con,
                    list(...)))
        
        poped.db <- add_to_poped_db(poped.db, ps_tbl, par = output$par)

      }
      if(cur_meth=="LS"){
        cat("*******************************************\n")
        cat("Running Line Search Optimization\n")
        cat("*******************************************\n")
        
        # Get parameters and space
        ps_tbl <- get_par_and_space_optim(poped.db,
                                          opt_xt=opt_xt,
                                          opt_a=opt_a,
                                          opt_x=opt_x,
                                          opt_samps=opt_samps,
                                          opt_inds=opt_inds,
                                          transform_parameters=F)
        
        
        
        # handle control arguments
        con <- list(trace = trace, 
                    parallel=parallel,
                    parallel_type=parallel_type,
                    num_cores = num_cores)
        con[names(control$LS)] <- control$LS
        
        # handle replicates
        if(!allow_replicates)
          con$replicates_index=ps_tbl$type
        if(!allow_replicates_xt | !allow_replicates_a)            
          con$replicates_index <- seq(1,length(ps_tbl$par))
        if(!allow_replicates_xt)
          con$replicates_index[ps_tbl$type=="xt"]="xt"
        if(!allow_replicates_a)
          con$replicates_index[ps_tbl$type=="a"]="a"

        tmp_ofv_fun <- function(par,...){ofv_optim(par,ps_tbl,poped.db,
                                                   d_switch=d_switch,
                                                   ED_samp_size=ED_samp_size,
                                                   bLHS=bLHS,
                                                   use_laplace=use_laplace,
                                                   ofv_calc_type=ofv_calc_type,
                                                   fim.calc.type=fim.calc.type,
                                                   ofv_fun = ofv_fun_user,
                                                   ...)}
        output <- 
          do.call(optim_LS,
                  c(list(par=ps_tbl$par,
                         fn=tmp_ofv_fun,
                         lower=ps_tbl$lower,
                         upper=ps_tbl$upper,
                         allowed_values = ps_tbl$allowed_values,
                         maximize=maximize
                  ),
                  con,
                  list(...)))
        
        poped.db <- add_to_poped_db(poped.db, ps_tbl, par = output$par)
      }
      if(cur_meth=="BFGS"){
        
        cat("*******************************************\n")
        cat("Running BFGS Optimization\n")
        cat("*******************************************\n")
        
        
        
        # Get parameters and space
        ps_tbl <- get_par_and_space_optim(poped.db,
                                          opt_xt=opt_xt,
                                          opt_a=opt_a,
                                          opt_x=opt_x,
                                          opt_samps=opt_samps,
                                          opt_inds=opt_inds,
                                          transform_parameters=T,
                                          cont_cat = "cont",
                                          warn_when_none=F)
        
       
        if(nrow(ps_tbl)==0){
          cat("\nNo continuous variables to optimize, BFGS Optimization skipped\n\n")
          next
        }
        
        if(!allow_replicates_xt | !allow_replicates_a | !allow_replicates){
          msg <- stringr::str_glue(
            'One or more of allow_replicates, allow_replicates_a, allow_replicates_xt are set to FALSE ', 
            'with design parameters that will be optimized accross a continuous space. ',
            'These options will be ignored for this optimization method and the ',
            'continuous optimization parameters.'
          )
          message(msg)
        }
          
        if(trace) trace_optim=3
        if(is.numeric(trace)) trace_optim = trace
        #if(trace==2) trace_optim = 4
        #if(trace==3) trace_optim = 5
        #if(trace==4) trace_optim = 6
        
        # handle control arguments
        con <- list(trace=trace_optim)
        nmsC <- names(con)
        con[(namc <- names(control$BFGS))] <- control$BFGS
        fnscale=-1
        if(!maximize) fnscale=1
        if(is.null(con[["fnscale"]])) con$fnscale <- fnscale
        #if (length(noNms <- namc[!namc %in% nmsC])) warning("unknown names in control: ", paste(noNms, collapse = ", "))
        
        tmp_ofv_fun <- function(par,...){ofv_optim(par,ps_tbl,poped.db,
                                                   d_switch=d_switch,
                                                   ED_samp_size=ED_samp_size,
                                                   bLHS=bLHS,
                                                   use_laplace=use_laplace,
                                                   ofv_calc_type=ofv_calc_type,
                                                   fim.calc.type=fim.calc.type,
                                                   ofv_fun = ofv_fun_user,
                                                   ...)}
        output <- optim(par=ps_tbl$par,
                        fn=tmp_ofv_fun,
                        gr=NULL,
                        ...,
                        #par_full=par_full,
                        # only_cont=T,
                        method = "BFGS",
                        control=con)
        
        output$ofv <- output$value
        poped.db <- add_to_poped_db(poped.db, ps_tbl, par=output$par)
        
        fprintf('\n')
        if(fn!="") fprintf(fn,'\n')
      }
      
      if(cur_meth=="GA"){
        
        cat("*******************************************\n")
        cat("Running Genetic Algorithm (GA) Optimization\n")
        cat("*******************************************\n")
        
        if (!requireNamespace("GA", quietly = TRUE)) {
          stop("GA package needed for this function to work. Please install it.",
               call. = FALSE)
        }
        
        # Get parameters and space
        ps_tbl <- get_par_and_space_optim(poped.db,
                                          opt_xt=opt_xt,
                                          opt_a=opt_a,
                                          opt_x=opt_x,
                                          opt_samps=opt_samps,
                                          opt_inds=opt_inds,
                                          transform_parameters=F,
                                          cont_cat = "cont",
                                          warn_when_none=F)
        
        
        if(nrow(ps_tbl)==0){
          cat("\nNo continuous variables to optimize, GA Optimization skipped\n\n")
          next
        }
        
        if(!allow_replicates_xt | !allow_replicates_a | !allow_replicates){
          msg <- stringr::str_glue(
            'One or more of allow_replicates, allow_replicates_a, allow_replicates_xt are set to FALSE ', 
            'with design parameters that will be optimized accross a continuous space. ',
            'These options will be ignored for this optimization method and the ',
            'continuous optimization parameters.'
          )
          message(msg)
        }
        
        # handle control arguments
        parallel_ga <- parallel
        if(!is.null(num_cores))  parallel_ga <- num_cores
        if(!is.null(parallel_type))  parallel_ga <- parallel_type
        
        con <- list(parallel=parallel_ga)
        #dot_vals <- dots(...)
        #if(is.null(dot_vals[["monitor"]])){
          #if(packageVersion("GA")>="3.0.2" && packageVersion("GA")<"3.1.1") con$monitor <- GA::gaMonitor2
          #if(packageVersion("GA")>="3.1.1") con$monitor <- GA::gaMonitor
        #}
        nmsC <- names(con)
        con[(namc <- names(control$GA))] <- control$GA
        #if (length(noNms <- namc[!namc %in% nmsC])) warning("unknown names in control: ", paste(noNms, collapse = ", "))
        
        tmp_ofv_fun <- function(par,...){ofv_optim(par,ps_tbl,poped.db,
                                                   d_switch=d_switch,
                                                   ED_samp_size=ED_samp_size,
                                                   bLHS=bLHS,
                                                   use_laplace=use_laplace,
                                                   ofv_calc_type=ofv_calc_type,
                                                   fim.calc.type=fim.calc.type,
                                                   ofv_fun = ofv_fun_user,
                                                   ...)}
        if(!maximize) tmp_ofv_fun <- function(par,...){-ofv_optim(par,ps_tbl,poped.db,
                                                                  d_switch=d_switch,
                                                                  ED_samp_size=ED_samp_size,
                                                                  bLHS=bLHS,
                                                                  use_laplace=use_laplace,
                                                                  ofv_calc_type=ofv_calc_type,
                                                                  fim.calc.type=fim.calc.type,
                                                                  ofv_fun = ofv_fun_user,
                                                                  ...)}
        if(packageVersion("GA")<"3.1.1")
          output_ga <- do.call(GA::ga,c(list(type = "real-valued", 
                                             fitness = tmp_ofv_fun,
                                             #par_full=par_full,
                                             # only_cont=T,
                                             min=ps_tbl$lower,
                                             max=ps_tbl$upper,
                                             suggestions=ps_tbl$par),
                                        #allowed_values = allowed_values),
                                        con,
                                        list(...)))
        if(packageVersion("GA")>="3.1.1")
          output_ga <- do.call(GA::ga,c(list(type = "real-valued", 
                                             fitness = tmp_ofv_fun,
                                             #par_full=par_full,
                                             # only_cont=T,
                                             lower=ps_tbl$lower,
                                             upper=ps_tbl$upper,
                                             suggestions=ps_tbl$par),
                                        #allowed_values = allowed_values),
                                        con,
                                        list(...)))
        output$ofv <- output_ga@fitnessValue
        if(!maximize) output$ofv <- -output$ofv
        
        
        output$par <- output_ga@solution
        
        poped.db <- add_to_poped_db(poped.db, ps_tbl, par=output$par)
        
        fprintf('\n')
        if(fn!="") fprintf(fn,'\n')
      }
      if(cur_meth=="opt_n_dist"){
        
        cat("*******************************************\n")
        cat("Running distribution of individuals optimization\n")
        cat("*******************************************\n")
          
        
        # Get parameters and space
        # ps_tbl <- get_par_and_space_optim(poped.db,
        #                                   opt_xt=opt_xt,
        #                                   opt_a=opt_a,
        #                                   opt_x=opt_x,
        #                                   opt_samps=opt_samps,
        #                                   opt_inds=opt_inds,
        #                                   transform_parameters=F,
        #                                   cont_cat = "cont",
        #                                   warn_when_none=F)
        # 
        
        # handle control arguments
        con <- control$opt_n_dist
        
        output_opt_n_dist <- do.call(optimize_groupsize,
                                     c(list(poped.db=poped.db),
                                       con,
                                       list(...)))
        
        poped.db <- create.poped.database(poped.db,
                              groupsize = output_opt_n_dist$opt_n_per_group,
                              mingroupsize = output_opt_n_dist$opt_n_per_group,
                              maxgroupsize = output_opt_n_dist$opt_n_per_group)
        
        fprintf('\n')
        if(fn!="") fprintf(fn,'\n')
      }
      
    }
    
    if(!loop_methods){
      stop_crit <- TRUE
    } else {
      cat("*******************************************\n")
      cat("Stopping criteria testing\n")
      cat("(Compare between start of iteration and end of iteration)\n")
      cat("*******************************************\n")
      
      # relative difference
      rel_diff <- (output$ofv - ofv_init)/abs(ofv_init)
      abs_diff <- (output$ofv - ofv_init)
      fprintf("Difference in OFV:  %.3g\n",abs_diff)
      fprintf("Relative difference in OFV:  %.3g%%\n",rel_diff*100)
      
      # efficiency
      
      
      eff <- efficiency(ofv_init, output$ofv, poped.db,...)
      fprintf("Efficiency: \n  (%s) = %.5g\n",attr(eff,"description"),eff)
      #cat("Efficiency: \n  ", attr(eff,"description"), sprintf("%.5g",eff), "\n")
      #if(eff<=stop_crit_eff) stop_crit <- TRUE
      
      #cat("eff: ",sprintf("%.3g",(output$ofv - ofv_init)/p), "\n")
      #cat("eff: ",sprintf("%.3g",(exp(output$ofv)/exp(ofv_init))^(1/p)), "\n")
      
      
      compare <-function(crit,crit_stop,maximize,inv=FALSE,neg=FALSE,text=""){
        
        if(is.null(crit_stop)){
          #cat("  Stopping criteria not defined\n")
          return(FALSE)
        }
        cat("\n",text,"\n")
        if(is.nan(crit)){
          fprintf("  Stopping criteria using 'NaN' as a comparitor cannot be used\n")
          return(FALSE)
        } 
        
        comparitor <- "<="
        if(!maximize) comparitor <- ">="
        if(inv) crit_stop <- 1/crit_stop
        if(neg) crit_stop <- -crit_stop
        fprintf("  Is (%0.5g %s %0.5g)? ",crit,comparitor,crit_stop)
        res <- do.call(comparitor,list(crit,crit_stop))
        if(res) cat("  Yes.\n  Stopping criteria achieved.\n")
        if(!res) cat("  No.\n  Stopping criteria NOT achieved.\n")
        return(res)
        #if(maximize) cat("Efficiency stopping criteria (lower limit) = ",crit_stop, "\n")
        #if(!maximize) cat("Efficiency stopping criteria (upper limit) = ",1/crit_stop, "\n")
        
        #if(maximize) return(crit <= crit_stop)
        #if(!maximize) return(crit >= 1/crit_stop)
      } 
      
      if(all(is.null(c(stop_crit_eff,stop_crit_rel,stop_crit_diff)))){
        cat("No stopping criteria defined")
      } else {
        
        
        stop_eff <- compare(eff,stop_crit_eff,maximize,inv=!maximize,
                            text="Efficiency stopping criteria:")
        
        stop_abs <- compare(abs_diff,stop_crit_diff,maximize,neg=!maximize,
                            text="OFV difference stopping criteria:")
        
        stop_rel <- compare(rel_diff,stop_crit_rel,maximize,neg=!maximize,
                            text="Relative OFV difference stopping criteria:")
        
        if(stop_eff || stop_rel || stop_abs) stop_crit <- TRUE
        
        if(stop_crit){
          cat("\nStopping criteria achieved.\n")
        } else {
          cat("\nStopping criteria NOT achieved.\n")
        }
        cat("\n")
      }
    }
    
  } # end of total loop 
  
  if(!(fn=="")) sink()
  
  FIM <-calc_ofv_and_fim(poped.db,
                         ofv=output$ofv,
                         fim=0, 
                         d_switch=d_switch,
                         ED_samp_size=ED_samp_size,
                         bLHS=bLHS,
                         use_laplace=use_laplace,
                         ofv_calc_type=ofv_calc_type,
                         fim.calc.type=fim.calc.type,
                         ofv_fun = ofv_fun_user,
                         ...)[["fim"]]
  
  if(use_laplace) poped.db$settings$ofv_calc_type <- 1
  time_value <- 
    blockfinal(fn=fn,fmf=FIM,
               dmf=output$ofv,
               groupsize=poped.db$design$groupsize,
               ni=poped.db$design$ni,
               xt=poped.db$design$xt,
               x=poped.db$design$x,
               a=poped.db$design$a,
               model_switch=poped.db$design$model_switch,
               poped.db$parameters$param.pt.val$bpop,
               poped.db$parameters$param.pt.val$d,
               poped.db$parameters$docc,
               poped.db$parameters$param.pt.val$sigma,
               poped.db,
               fmf_init=fmf_init,
               dmf_init=dmf_init,
               ...)
  
  
  #  }
  #}
  results <- list( ofv= output$ofv, FIM=FIM, initial=list(ofv=dmf_init, FIM=fmf_init, poped.db=poped.db_init), 
                   run_time=time_value, poped.db = poped.db )
  class(results) <- "poped_optim"
  return(invisible(results)) 
}

