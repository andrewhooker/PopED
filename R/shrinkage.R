#' Predict shrinkage of empirical Bayes estimates (EBEs) in a population model
#'
#' @param poped.db A PopED database
#' @param num_sim_ids If \code{use_mc=TRUE}, how many individuals should be
#'   simulated to make the computations.
#' @param use_mc Should the calculation be based on monte-carlo simulations. If
#'   not then then a first order approximation is used
#' @param use_purrr If \code{use_mc=TRUE} then should the method use the package
#'   purrr in calculations?  This may speed up computations (potentially).
#'
#' @return The shrinkage computed in variance units, standard deviation units
#'   and the relative standard errors of the EBEs.
#' @export
#'
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_shrinkage.R
#' @references \enumerate{ 
#'   \item Combes, F. P., Retout, S.,
#'   Frey, N., & Mentre, F. (2013). Prediction of shrinkage of individual
#'   parameters using the Bayesian information matrix in non-linear mixed effect
#'   models with evaluation in pharmacokinetics. Pharmaceutical Research, 30(9),
#'   2355-67. \url{https://doi.org/10.1007/s11095-013-1079-3}. 
#'   \item Hennig, S., Nyberg, J., Fanta, S., Backman, J.
#'   T., Hoppu, K., Hooker, A. C., & Karlsson, M. O. (2012). Application of the
#'   optimal design approach to improve a pretransplant drug dose finding design
#'   for ciclosporin. Journal of Clinical Pharmacology, 52(3), 347-360.
#'   \url{https://doi.org/10.1177/0091270010397731}. 
#'   }

shrinkage <- function(poped.db,
                      use_mc = FALSE,
                      num_sim_ids = 1000,
                      use_purrr = FALSE){
  
  # if (poped.db$design$m > 1) {
  #   warning("Shrinkage should only be computed for a single arm, please adjust your script accordingly.")
  # }
  
  group <- type <- NULL
  
  # tranform random effects to fixed effects 
  tmp_fg <- poped.db$model$fg_pointer
  if(is.character(tmp_fg)) tmp_fg <- eval(parse(text=tmp_fg))
  largest_bpop <- find.largest.index(tmp_fg,lab = "bpop")
  largest_b <- find.largest.index(tmp_fg,lab = "b")
  largest_eps <- find.largest.index(poped.db$model$ferror_pointer,lab = "epsi",mat=T,mat.row = F)
  replace_fun <- function(txt){
    #extract number
    num <- stringr::str_extract(txt,"\\d+") %>% as.numeric() %>% "+"(.,largest_bpop)
    stringr::str_replace(txt,"\\d+",paste0(num))
  }
  txt <- capture.output(body(tmp_fg))
  txt <- stringr::str_replace_all(txt,"b\\[(\\d+)",replace_fun)
  txt <- stringr::str_replace_all(txt,"b\\[","bpop\\[")
  env <- environment()
  body(tmp_fg,envir=env) <- parse(text=txt)
  #environment(sfg_tmp) <- env
  
  # fix population parameters
  # add IIV and IOV as fixed effects
  # change to one individual
  extra_bpop <- matrix(0,nrow = largest_b,ncol = 3)
  bpop_tmp <-rbind(poped.db$parameters$bpop,extra_bpop)
  notfixed_bpop_tmp <- c(rep(0,largest_bpop),rep(1,largest_b))
  poped.db_sh <- create.poped.database(poped.db,
                                       fg_fun=tmp_fg,
                                       bpop=bpop_tmp, 
                                       nbpop=nrow(bpop_tmp),
                                       notfixed_bpop = notfixed_bpop_tmp,
                                       d=matrix(nrow = 0,ncol = 3),
                                       notfixed_d = matrix(nrow = 0,ncol = 3),
                                       NumRanEff = 0,
                                       notfixed_sigma = c(rep(0,largest_eps)),
                                       groupsize=1,
                                       mingroupsize = 1,
                                       mintotgroupsize = 1)
  
  # Compute FIM_map 
  fulld = getfulld(poped.db$parameters$d[,2,drop=F],poped.db$parameters$covd)
  
  # get just groupwise values as well
  db_list <- list()
  #db_list[["all"]] <- poped.db_sh
  num_groups <- poped.db_sh$design$m
  for(i in 1:num_groups){
    tmp_design <- poped.db_sh$design
    tmp_design$m <- 1
    for(nam in names(tmp_design)[names(tmp_design)!="m"]){
      tmp_design[[nam]] <- tmp_design[[nam]][i,,drop=F]
    }

    tmp_db <- poped.db_sh
    tmp_db$design <- tmp_design
    
    db_list[[paste0("grp_",i)]] <- tmp_db
  }
  
  
  #browser()
  
  #out_list <- list()
  out_df <- c()
  #for(i in 1:1){
  for(i in 1:length(db_list)){
    
    poped.db_sh <- db_list[[i]]
    
    if(use_mc){
      # create list of eta values from pop model 
      if(any(size(fulld)==0)){
        b_sim_matrix = zeros(num_sim_ids,0)
      } else {
        b_sim_matrix = mvtnorm::rmvnorm(num_sim_ids,sigma=fulld)            
      }
      
      # create list of occasion eta values from pop model 
      NumOcc=poped.db$parameters$NumOcc
      if(NumOcc!=0) warning("Shrinkage not computed for occasions\n")
      fulldocc = getfulld(poped.db$parameters$docc[,2,drop=F],poped.db$parameters$covdocc)
      bocc_sim_matrix = zeros(num_sim_ids*NumOcc,length(poped.db$parameters$docc[,2,drop=F]))
      if(nrow(fulldocc)!=0) bocc_sim_matrix = mvtnorm::rmvnorm(num_sim_ids*NumOcc,sigma=fulldocc)
      
      #now use these values to compute FIMmap
      if(use_purrr){
        fim_mean <- b_sim_matrix %>% t() %>% data.frame() %>% as.list() %>% 
          purrr::map(function(x){
            x_tmp <- matrix(0,nrow = largest_b,ncol = 3)
            x_tmp[,2] <- t(x)
            bpop_tmp <-rbind(poped.db$parameters$bpop,x_tmp)
            poped.db_sh_tmp <- create.poped.database(poped.db_sh, bpop=bpop_tmp)
            evaluate.fim(poped.db_sh_tmp)
          }) %>% simplify2array() %>% apply(1:2, mean)
      } else {
        fim_tmp <- 0
        for(i in 1:num_sim_ids){
          x_tmp <- matrix(0,nrow = largest_b,ncol = 3)
          x_tmp[,2] <- t(b_sim_matrix[i,])
          bpop_tmp <-rbind(poped.db$parameters$bpop,x_tmp)
          poped.db_sh_tmp <- create.poped.database(poped.db_sh, bpop=bpop_tmp)
          fim_tmp <- fim_tmp+evaluate.fim(poped.db_sh_tmp)
        }
        fim_mean <- fim_tmp/num_sim_ids
      }    
      fim_map <- fim_mean + inv(fulld)
    } else {
      fim_map <- evaluate.fim(poped.db_sh)+ inv(fulld)
    }
    
    rse <- get_rse(fim_map,poped.db = poped.db_sh)
    names(rse) <- paste0("d[",1:length(rse),"]")
    
    # shrinkage on variance scale
    shrink <- diag(inv(fim_map)%*%inv(fulld))
    names(shrink) <-  names(rse)
    # shrinkage on SD scale
    var_n <- (1-shrink)*diag(fulld)
    shrink_sd <- 1-sqrt(var_n)/sqrt(diag(fulld))
    names(shrink_sd) <-  names(rse)
    
    out_df_tmp <- tibble::as.tibble(rbind(shrink,shrink_sd,rse))
    out_df_tmp$type <- c("shrink_var","shrink_sd","se")
    out_df_tmp$group <- c(names(db_list[i]))
    out_df <- rbind(out_df,out_df_tmp)
    #out_list[[names(db_list[i])]] <- list(shrinkage_var=shrink,shrinkage_sd=shrink_sd, rse=rse)
  }
  
  #return(list(shrinkage_var=shrink,shrinkage_sd=shrink_sd, rse=rse))
  #return(out_list)
  out_df <- dplyr::arrange(out_df,dplyr::desc(type),group)
  
  if(poped.db$design$m > 1){
    weights <- poped.db$design$groupsize/sum(poped.db$design$groupsize)
    data_tmp <- out_df 
    data_tmp[data_tmp==1] <- NA 
    data_tmp <- data_tmp %>% dplyr::group_by(type) %>% 
      dplyr::summarise_at(dplyr::vars(dplyr::starts_with('d[')),
                          dplyr::funs(stats::weighted.mean(., weights,na.rm = T)))
    data_tmp$group <- "all_groups"
    out_df <- rbind(out_df,data_tmp)
    out_df <- dplyr::arrange(out_df,dplyr::desc(type),group)
  }
  return(out_df)
  
  
}