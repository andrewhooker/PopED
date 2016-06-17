#' Compute the monte-carlo mean of a function
#'
#'Function computes the monte-carlo mean of a function by varying the parameter inputs to the function
#'
#' @param ofv_fcn A function with poped.db as the first input
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param doccdescr Matrix defining the IOV.
#' per row (row number = parameter_number) we should have:
#' \itemize{
#' \item column 1 the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
#'  3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal)
#' \item column 2  defines the mean of the variance.
#' \item column 3 defines the variance of the distribution (or length of uniform distribution).
#' }
#' @param user_distribution_pointer Function name for user defined distributions for E-family designs 
#' @return The mean of the function evaluated at different parameter values.
#' @export
#'
# @examples
 
mc_mean <- function(
  ofv_fcn,
  poped.db,
  bpopdescr=poped.db$parameters$bpop, 
  ddescr=poped.db$parameters$d,
  doccdescr=poped.db$parameters$d,
  user_distribution_pointer=poped.db$model$user_distribution_pointer,
  ED_samp_size=poped.db$settings$ED_samp_size,
  bLHS=poped.db$settings$bLHS,
  ...) 
{
  ofv_sum=0
  
  
  bpop_gen  <-  pargen(bpopdescr,
                       user_distribution_pointer,
                       ED_samp_size,
                       bLHS,
                       NULL,#zeros(1,0),
                       poped.db)
  
  d_gen <- pargen(ddescr,
                  user_distribution_pointer,
                  ED_samp_size,
                  bLHS,
                  NULL,
                  poped.db)
  
  docc_gen <- pargen(doccdescr,
                     user_distribution_pointer,
                     ED_samp_size,
                     bLHS,
                     NULL,
                     poped.db)
  
  poped.db_tmp <- poped.db
  for(ct in 1:ED_samp_size){
    poped.db_tmp$parameters$bpop[,2] <- bpop_gen[ct,]
    poped.db_tmp$parameters$d[,2] <- d_gen[ct,]
    poped.db_tmp$parameters$docc[,2] <- docc_gen[ct,]
    
    dmf_tmp <- do.call(ofv_fcn,list(poped.db_tmp,...))
    
    ofv_sum=ofv_sum+dmf_tmp
  }
  ofv_mean=ofv_sum/poped.db$settings$ED_samp_size
  return(ofv_mean)
}