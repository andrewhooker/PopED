#' Optimize a design defined in a PopED database
#'  
#' 
#' Optimize a design defined in a PopED database using the objective function
#' described in the database (or in the arguments to this function). The 
#' function works for both discrete and 
#' continuous optimization variables. 
#' 
#' This function takes information from the PopED database supplied as an 
#' argument. The PopED database supplies information about the the model, 
#' parameters, design and methods to use. Some of the arguments coming from the 
#' PopED database can be overwritten; if they are supplied then they are used 
#' instead of the arguments from the PopED database.
#' 
#' If more than one optimization method is 
#' specified then the methods are run in series.  If \code{loop_methods=TRUE} 
#' then the series of optimization methods will be run for \code{iter_max} 
#' iterations, or until the efficiency of the design after the current series 
#' (compared to the start of the series) is less than \code{stop_crit_eff}.
#' 
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams Dtrace
#' @inheritParams calc_ofv_and_fim
#' @inheritParams optim_LS
#' @param ... arguments passed to other functions.
#' @param control Contains control arguments for each method specified.
#' @param method A vector of optimization methods to use in a sequential 
#'   fashion.  Options are \code{c("ARS","BFGS","LS","GA")}. \code{c("ARS")} is 
#'   for Adaptive Random Search \code{\link{optim_ARS}}.  \code{c("LS")} is for 
#'   Line Search \code{\link{optim_LS}}. \code{c("BFGS")} is for Method 
#'   "L-BFGS-B" from \code{\link[stats]{optim}}. \code{c("GA")} is for the 
#'   genetic algorithm from \code{\link[GA]{ga}}.
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
#' @param allow_replicates Should the algorithm allow optimized design components to have the same value? If FALSE then
#' all discrete optimizations will not allow replicates within variable types 
#' (equivalent to \code{allow_replicates_xt=FALSE} and \code{allow_replicates_a=FALSE}). 
#' @param allow_replicates_xt Should the algorithm allow optimized \code{xt} design components to have the same value? If FALSE then
#' all discrete optimizations will not allow replicates.
#' @param allow_replicates_a Should the algorithm allow optimized \code{a} design components to have the same value? If FALSE then
#' all discrete optimizations will not allow replicates.
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
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_poped_optim.R
#' @export

poped_optim <- function(poped.db,
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
                        approx_type=poped.db$settings$iApproximationMethod,
                        d_switch=poped.db$settings$d_switch,
                        ED_samp_size = poped.db$settings$ED_samp_size,
                        bLHS=poped.db$settings$bLHS,
                        use_laplace=poped.db$settings$iEDCalculationType,
                        out_file="",
                        parallel=F,
                        parallel_type=NULL,
                        num_cores = NULL,
                        mrgsolve_model = NULL,
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
  
  
 
  #------------ update argument list with called arguments
  arg.list <- formals()
  arg.names <- names(arg.list)
  arg.names <- arg.names[-length(arg.names)]
  new.arg.list <- vector("list",length(arg.names))
  names(new.arg.list) <- arg.names
  for (argnam in arg.names){
    tmp <- get(argnam)
    if (!is.null(tmp)){
      new.arg.list[[argnam]]=tmp
    }
  }

  dot_vals <- dots(...)
  #dot_vals <- list(...)
  optim_ver <- dot_vals[["optim_ver"]]
  dot_vals[["optim_ver"]] <- NULL
  if(is.null(optim_ver)) optim_ver <- 3
  if(optim_ver!=3 & (!allow_replicates | 
                     !allow_replicates_xt | 
                     !allow_replicates_a)){
    msg <- stringr::str_glue(
      'One or more of allow_replicates, allow_replicates_a, allow_replicates_xt are set to FALSE ', 
      'This will only work with optim_ver=3'
    )
    message(msg)
  } 
  
  if(optim_ver==1) results <- do.call(poped_optim_1,c(new.arg.list,dot_vals))
  if(optim_ver==2) results <- do.call(poped_optim_2,c(new.arg.list,dot_vals))
  if(optim_ver==3) results <- do.call(poped_optim_3,c(new.arg.list,dot_vals))
  
  return(invisible(results)) 
}

