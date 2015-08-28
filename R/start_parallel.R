start_parallel <- function(parallel=TRUE,
                           num_cores=NULL,
                           parallel_type=NULL,
                           seed=NULL,
                           dlls=NULL,
                           ...)
{
  # Start parallel computing for poped package
  # edited from GA package version startParallel.R
  
  if(is.null(num_cores)) num_cores <- parallel::detectCores()
  if(is.null(parallel_type)) parallel_type <- if(.Platform$OS.type == "windows") 
    "snow" else "multicore"
  attr(parallel, "type") <- parallel_type
  attr(parallel, "cores") <- num_cores
  
  # start "parallel backend" if needed
  if(parallel){ 
    if(parallel_type == "snow"){ 
      # snow functionality on Unix-like systems & Windows
      cl <- parallel::makeCluster(num_cores, ...)
      attr(parallel, "cluster") <- cl
      
      # export parent environment
      varlist <- ls(envir = parent.frame(), all.names = TRUE)
      varlist <- varlist[varlist != "..."]
      list(...) #evaluate any promises
      parallel::clusterExport(cl, varlist = varlist,
                              # envir = parent.env(environment())
                              envir = parent.frame() )
      
      # export global environment (workspace)
      parallel::clusterExport(cl, 
                              varlist = ls(envir = globalenv(), 
                                           all.names = TRUE),
                              envir = globalenv())
      
      # load current packages in workers
      pkgs <- .packages()
      foo <- lapply(pkgs, function(pkg) 
        parallel::clusterCall(cl, library, package = pkg, 
                              character.only = TRUE))
      if(!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
      
      if(!is.null(dlls)){
        for(i in dlls){
          parallel::clusterCall(cl, 
                                dyn.load,
                                x=paste0(i,.Platform$dynlib.ext))
        }
      }
      
    } else if(parallel_type == "multicore") { 
      if(!is.null(seed)){
        RNGkind("L'Ecuyer-CMRG") 
        set.seed(seed)
        #set.seed(seed, "L'Ecuyer")
      } 
      # multicore functionality on Unix-like systems
    }
    else { stop("Only 'snow' and 'multicore' clusters allowed!") }
  }
  
  return(parallel)
}