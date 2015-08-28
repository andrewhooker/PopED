#' Optimization Using a Line Search Algorithm. 
#' 
#' \code{optim_LS} performs sequential grid search optimization of an arbitrary function with respect 
#' to each of the parameters to be optimized over. 
#' The function works for both discrete and continuous optimization parameters 
#' and allows for box-constraints and sets of allowed values. 
#' 
#' @references \enumerate{
#' \item M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a software fir optimal 
#' experimental design in population kinetics", Computer Methods and Programs in Biomedicine, 74, 2004.
#' \item J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. Hooker, "PopED: An extended, 
#' parallelized, nonlinear mixed effects models optimal design tool",  
#' Computer Methods and Programs in Biomedicine, 108, 2012.
#' }
#' 
#' @family Optimize
#' 
## @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
## @example tests/testthat/examples_fcn_doc/examples_optim_LS.R

optim_LS <- function(par,
                     fn,
                     lower,
                     upper,
                     allowed_values=NULL,
                     #constraints=NULL,
                     #control=NULL,
                     iter_chunk = NULL,
                     trace = TRUE,
                     trace_iter=5,
                     par_grouping=NULL,
                     minimize=F, # default is to maximize
                     parallel=F,
                     parallel_type=NULL,
                     num_cores = NULL,
                     seed=round(runif(1,0,10000000)),
                     n_steps=50,
                     ...){
  
  # constratints

  #---------------- start trace
  if((trace)){
    tic(name=".ls_savedTime")
  }
  
  #--------------- checks
  if((is.null(lower) || is.null(upper)) && is.null(allowed_values)){
    stop("At least 'lower' and 'upper' or 'allowed_values' should be supplied.")
  }   
  
  #----------------- initialization 
  par0=par
  fn1 <- function(par) fn(par, ...)  
  npar <- max(c(length(lower),length(upper),length(allowed_values)))
  ofv_opt <- NULL
  par_opt <- par
  if(!is.null(seed)) set.seed(seed)
  
  if(parallel){
    parallel <- start_parallel(parallel,seed=seed,parallel_type=parallel_type,num_cores=num_cores,...) 
    on.exit(if(parallel && (attr(parallel,"type")=="snow")) parallel::stopCluster(attr(parallel,"cluster")))
  }   
  if(is.null(iter_chunk)) if(parallel) iter_chunk <- attr(parallel,"cores") else iter_chunk <- 1
  
  # continuous and discrete parameters
  par_type <- rep("cont",npar)
  if(!is.null(allowed_values)){
    for(k in 1:npar){
      if(!is.na(allowed_values[[k]]) && length(allowed_values[[k]]>0)){
        par_type[k] <- "cat"          
      }
    }
  }
  
  compare <-function(a,b) a > b
  if(minimize)  compare <-function(a,b) a < b
  
  
  # ------------ generate and evaluate new parameter sets
  gen_par_ofv <- function (it, par, ...) {
    need_new_par <- TRUE
    new_par_it <- 0
    
    #---------------- generate new parameter set
    while(need_new_par && (new_par_it < new_par_max_it)){
      
      if(it==1 && !is.null(par)){
        need_new_par <- FALSE
        break
      }
      # generate new continuous parameters
      if(is.null(par)) par  <- lower+(upper - lower)/2
      par[par_type=="cont"]=par[par_type=="cont"]+dpar[par_type=="cont"]/ff*rnorm(length(par[par_type=="cont"]))
      
      # generate new discrete parameters
      par[par_type=="cat"] <- vapply(allowed_values[par_type=="cat"],resample,c(0))
      
      # Group samples that should be grouped
      if(!is.null(par_grouping)){
        for(k in unique(par_grouping)){
          par[par_grouping==k] <- resample(par[par_grouping==k])
        }
      }
      
      # set to boundary if beyond the boundary
      if(!all(is.null(upper))) par=par-((par>upper)*(par-upper))
      if(!all(is.null(upper))) par=par+((par<lower)*(lower-par))
      
      # check if design has already been tested
      #need_new_par <- any(sapply(sapply(par_vec, length)!=0,function(x) all(x == par)))
      #par_vec[sapply(par_vec, length)!=0,]
      need_new_par <- any(sapply(par_vec[sapply(par_vec, length)!=0,],function(x) all(x==par)))
      #browser()
      new_par_it <- new_par_it + 1
    } # end while  
    
    #------------- evaluate new parameter sets
    if(need_new_par){
      ofv <- NULL
      par <- NULL
    } else {
      ofv <- fn1(par)
    }
    return(list(ofv=ofv,par=par,need_new_par=need_new_par))
  } # end function
  
  for(it in 1:ceiling(iter/iter_chunk)){
    start_it <- (it-1)*iter_chunk+1
    stop_it <- min(it*iter_chunk,iter)
    it_seq <- start_it:stop_it
    
    # generate new parameters and OFVs
    if(parallel && (attr(parallel,"type")=="multicore")){
      if(!is.null(seed)) set.seed(seed+it)
      res <- mapply(c,parallel::mclapply(it_seq,gen_par_ofv,par_opt,mc.cores=attr(parallel, "cores")))
    } else if(parallel && (attr(parallel,"type")=="snow")){
      res <- mapply(c,parallel::parLapply(attr(parallel, "cluster"),it_seq,gen_par_ofv,par_opt))
    } else {
      res <- mapply(c,lapply(it_seq,gen_par_ofv,par_opt))  
    }
    
    res2 <- res[,!sapply(res["ofv",],is.null),drop=F]
    out <- res2[,which.max(res2["ofv",])]
    
    # check if a unique new parameter vector was generated
    if(is.null(out$ofv)){
      cat(paste0("Maximum number of duplicate parameter samples reached (new_par_max_it=",new_par_max_it,"), optimization stopped.\n"))
      break
    }
    
    ofv <- out$ofv
    par <- out$par
    
    par_vec[it_seq] <- res["par",] # save new parameter vectors in a list
    
    if((compare(ofv,ofv_opt) || is.null(ofv_opt))){
      par_opt <- par 
      ofv_opt <- ofv
      nullit=1
      runs <- 1
      #ff=1
    } else {
      nullit=nullit+(length(it_seq))
      runs <- runs+(length(it_seq))
    }
    if((nullit>=iter_adapt) ){# when to make the search area smaller
      ff=ff+1
      nullit=1
    }
    
    if((trace && any(rem(it_seq,trace_iter)==0))){
      if(length(it_seq)==1){ 
        cat(sprintf(paste0("It. %",wd_iter,"i"),start_it))
      } else {
        cat(sprintf(paste0("It. %",wd_iter,"i to %",wd_iter,"i"),start_it,stop_it))
      }
      cat(sprintf(" | OFV = %g",ofv_opt))
      if(trace==2) cat(" | opt. par. = ",par_opt)
      #cat(" | runs = ",runs)
      #cat(" | nullit = ",nullit)
      if(trace==3) cat(" | par tried = ",par)
      cat("\n")
    }
    
    if(runs>=max_run){
      cat(paste0("Maximum number of identical optimal values reached (max_run=",max_run,"), optimization stopped.\n"))
      break
    }     
  }
  
  #--------- Write results
  if((trace)){
    cat("\nTotal iterations:",stop_it,"\n")
    toc(name=".ls_savedTime")
    cat("\nFinal OFV = ", ofv_opt, "\n") 
    cat("Parameters:",par_opt, "\n\n")
  }
  
  return(list( par=par_opt,ofv=ofv_opt)) 
}
