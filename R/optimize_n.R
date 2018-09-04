#' Extract a normalized group FIM
#'
#' Extract an individual FIM for each group in a design 
#'
#' @param poped.db A PopED database
#' @param ... Arguments passed to \code{\link{evaluate.fim}}
#'
#' @return A list of FIMs, one for each group in a design.
#' @keywords internal
#'
extract_norm_group_fim <- function (poped.db,...) {
  num_groups <- poped.db$design$m
  fim_group <- list()
  for(i in 1:num_groups){
    tmp_design <- poped.db$design
    tmp_design$m <- 1
    for(nam in names(tmp_design)[names(tmp_design)!="m"]){
      tmp_design[[nam]] <- tmp_design[[nam]][i,,drop=F]
    }
    tmp_design$groupsize[,] <- 1
    tmp_design
    
    tmp_db <- poped.db
    tmp_db$design <- tmp_design
    
    tmp_fim <- evaluate.fim(tmp_db,...)
    fim_group[[paste0("grp_",i)]] <- tmp_fim
  }
  return(fim_group)
}


combine_norm_group_fim <- function(norm_group_fim,n_per_group){
  fim <- 0
  for(i in 1:length(norm_group_fim)){
    fim <- fim+n_per_group[i]*norm_group_fim[[i]]
  }
  return(fim)
}


## optimize on where the individuals should be
optimize_n_dist <- 
  function(poped.db,
           props = c(poped.db$design$groupsize/sum(poped.db$design$groupsize)),
           ...){
  
  # need to fix:
  # limits on the proportions to account for max and min values of N in each group
  # return actual values.
  
  
  if(sum(props) != 1) stop("The sum of the proportions are not equal to 1\n") # check that the proportions add up to one
  
  # Transform the parameters
  convert_to_lcp <- function(p){
    if(sum(p)!=1) stop("Probabilities do not add up to one")
    lcp=p
    p_rem <- 1
    for(i in 1:(length(p)-1)){
      lcp[i] <- boot::logit(p[i]/p_rem)
      p_rem <- p_rem-p[i]
    }
    lcp[length(p)] <- boot::logit(1)
    return(lcp)
  }
  
  convert_from_lcp <- function(p) {
    q=p
    p_rem <- 1
    for(i in 1:(length(p)-1)){
      q[i] <- p_rem * boot::inv.logit(p[i])
      p_rem <- p_rem - q[i]
    }
    q[length(p)] <- p_rem
    return(q)
  }
  
  # convert to and from logit of conditional probability
  # (foo <- convert_to_lcp(c(0.3,0.2,0.3, 0.2)))
  # convert_from_lcp(foo)
  # 
  # convert_from_lcp(c(0,0,0,0,0,0))
  
  # Optimization
  n_per_group <- poped.db$design$groupsize
  n_tot <- sum(n_per_group)
  #props <- c(n_per_group/n_tot)
  norm_group_fim <- extract_norm_group_fim(poped.db)
  
  ofv_fun <- function(props,norm_group_fim,n_tot, poped.db, ...){
    fim_tmp <- combine_norm_group_fim(norm_group_fim,props*n_tot)
    ofv_fim(fim_tmp,poped.db,...)
  }
  
  initial_ofv <- ofv_fun(props,norm_group_fim,n_tot, poped.db, ...)
  
  # ofv_fun(poped.db$design$groupsize/sum(poped.db$design$groupsize),norm_group_fim,
  #         sum(poped.db$design$groupsize),poped.db)
  
  # same as the typical calculation?
  # evaluate_design(poped.db)$ofv
  
  
  opt_fun <- function(lcp,...){
    ofv_fun(convert_from_lcp(c(lcp,Inf)),norm_group_fim,n_tot,poped.db,...)
  } 
  
  # convert proportions to logit of conditional probabilities
  lcp <- convert_to_lcp(props)
  opt_var <- lcp[-(length(lcp))]
  
  #if(length(opt_var)==1){method = c("BFGS")} else {method = c("Nelder-Mead")}
  #if(length(opt_var)==1){method = c("BFGS")} else {method = c("BFGS")}
  
  p <- optim(opt_var, opt_fun,...,control=list(fnscale=-1),method = c("BFGS"))
  final_props <- convert_from_lcp(c(p$par,Inf)) # the proportions in each group
  if(!all.equal(sum(final_props),1)) stop("something went wrong\n the sum of the optimized proportions are not equal to 1\n") # check that the proportions add up to one
  n_per_group_opt <- final_props*n_tot # the numbers of individuals in each group
  ofv_opt <- p$value # the new OFV value
  #initial_ofv <- evaluate_design(poped.db)$ofv # the original value
  
  return(list(initial_props=props,
               initial_ofv=initial_ofv,
               opt_props=final_props,
               opt_ofv=ofv_opt,
               opt_n_per_group=n_per_group_opt) )
  
}


## optimize HOW MANY n there should be, based on current design and a single parameter
optimize_n <- function(poped.db,
                       bpop_idx,
                       need_rse,
                       allowed_values = seq(poped.db$design$m,
                                            sum(poped.db$design$groupsize)*5,
                                            by=poped.db$design$m),
                       ...){

  n_per_group = poped.db$design$groupsize
  n_tot <- sum(n_per_group)
  props = c(n_per_group/n_tot)
  norm_group_fim <- extract_norm_group_fim(poped.db)

  ofv_fun <- function(n_tot){
    fim_tmp <- combine_norm_group_fim(norm_group_fim,props*n_tot)
    ofv <- get_rse(fim_tmp,poped.db,use_percent = F)[cumsum(poped.db$parameters$notfixed_bpop==1)[bpop_idx]] - need_rse
    if(ofv>0){
      ofv <- Inf
    } else {
      ofv <- abs(ofv)
    }
    return(ofv)
  }

  result <- optim_LS(n_tot, ofv_fun, allowed_values = allowed_values,maximize = F,trace = F)
  return(result)

}

