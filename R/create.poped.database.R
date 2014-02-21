#' Create a PopED database from a poped input file
#' 
#' This function takes the input file supplied by the user and creates a database that can then be used to 
#' run all other PopED functions.  The function supplies default values to elements of the database that are not specified in the
#' input file.  The function also computes needed quatities.
#' 
#' @param popedInput An input file to PopED 
#' @return A PopED database
#' @family poped_input
#' 
#' @examples 
#' \dontrun{
#' source("warfarin.POPED.R")
#' poped.db <- create.poped.database(warfarin.input())
#' }  
#' 

create.poped.database <- function(popedInput,
                                  ...){
  
  poped.db <- convert_popedInput(popedInput,...)
  poped.db <- convert_variables(poped.db)
  param.val <- get_all_params(poped.db)
  tmp.names <- names(param.val)
  eval(parse(text=paste(tmp.names,".val","<-","param.val$",tmp.names,sep="")))
  d.val <- d.val # for package check
  covd.val <- covd.val
  docc.val <- docc.val
  covdocc.val <- covdocc.val
  bpop.val <- bpop.val
  d_full=getfulld(d.val,covd.val)
  docc_full = getfulld(docc.val,covdocc.val)
  sigma_full = poped.db$sigma
  poped.db$param.pt.val$bpop <- bpop.val
  poped.db$param.pt.val$d <- d_full
  poped.db$param.pt.val$docc <- docc_full
  poped.db$param.pt.val$sigma <- sigma_full
  
  downsize.list <- downsizing_general_design(poped.db)
  tmp.names <- names(downsize.list)
  model_switch <- c()
  ni <- c()
  xt <- c()
  x <- c()
  a <- c()
  eval(parse(text=paste(tmp.names,"<-","downsize.list$",tmp.names,sep="")))
  poped.db$downsized.design$model_switch <- model_switch
  poped.db$downsized.design$ni <- ni
  poped.db$downsized.design$xt <- xt
  poped.db$downsized.design$x <- x
  poped.db$downsized.design$a <- a
  poped.db$downsized.design$groupsize <- poped.db$design$groupsize
  
  retargs <- fileparts(poped.db$output_file)
  poped.db$strOutputFilePath <- retargs[[1]]
  poped.db$strOutputFileName <- retargs[[2]]
  poped.db$strOutputFileExtension <- retargs[[3]]
  
  return(poped.db) 
}
