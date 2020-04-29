
#' Build PopED parameter function
#'
#' @param model 
#' @param covariates 
#' @param par_names 
#' @param etas 
#' @param no_etas 
#' @param env 
#'
#' @return
# @export
#' @importFrom codetools findGlobals

#'
#' @examples
#' build_sfg(model="ff.PK.1.comp.oral.md.CL")
#' 
#' etas <- c(Favail="exp",KA="exp",V="add",CL="exp")
#' build_sfg(model="ff.PK.1.comp.oral.md.CL",etas = etas)

build_sfg <- function (model="ff.PK.1.comp.oral.sd.CL",covariates=c("dose","tau"),par_names=NULL,
                       etas="exp", # can be exp, prop, add, none. can be one for all or per parameter
                       no_etas=c("F","Favail"),
                       env = parent.frame()) {
  
  ## get variable of function
  parameter_names_ff <- par_names
  if(is.null(parameter_names_ff)) parameter_names_ff <- codetools::findGlobals(eval(parse(text=model)),merge=F)$variables  
  
  ## create an empty function
  sfg_tmp <- function(x,a,bpop,b,bocc){}
  
  ## get body of function
  parameters <- "parameters=c("
  bpop_num <- 1
  b_num <- 1
  a_num <- 1
 
  cov_locs <- grep(paste0("^",paste(covariates,collapse="|"),"$"),
                   parameter_names_ff,ignore.case = TRUE)
  covariate_names <- parameter_names_ff[cov_locs]
  parameter_names <- parameter_names_ff[-cov_locs]
  
  # match names
  df <- NULL
  if(!missing(etas) & !is.null(names(etas))){
     if(all(parameter_names %in% names(etas))){
       df <- cbind(parameter_names,eta_mod = etas[parameter_names])
     }
  } 
  if(is.null(df)) df <- cbind(parameter_names,eta_mod = etas)
  
  if(!missing(no_etas) | (length(etas)==1)){
    df[parameter_names %in% no_etas,"eta_mod"] <- "none"
  }
  
  for(k in 1:nrow(df)){
    ending <- ", "
    if(k==nrow(df) & (length(covariate_names)==0)) ending <- ")"
    
    parameters <- paste0(parameters,
                        df[k,"parameter_names"], 
                        "=bpop[",bpop_num,"]")
    bpop_num <- bpop_num+1
    if(df[k,"eta_mod"]=="exp") parameters <- paste0(parameters,"*exp(b[",b_num,"])")
    if(df[k,"eta_mod"]=="prop") parameters <- paste0(parameters,"*(1+b[",b_num,"])")
    if(df[k,"eta_mod"]=="add") parameters <- paste0(parameters,"+ b[",b_num,"]")  
    if(df[k,"eta_mod"]=="none") parameters <- parameters  
    if(df[k,"eta_mod"]!="none") b_num <- b_num+1
    parameters <- paste0(parameters,ending)        
  }
  
  if(length(covariate_names)!=0){
    for(k in 1:length(covariate_names)){
      ending <- ", "
      if(k==length(covariate_names)) ending <- ")"
      parameters <- paste(parameters,covariate_names[k],"=a[",a_num,"]",ending,sep="")
      a_num <- a_num+1
    }  
  }
  
    
  ## add body of funciton and set environment
  body(sfg_tmp) <- parse(text=parameters)
  environment(sfg_tmp) <- env

  return(sfg_tmp)
}