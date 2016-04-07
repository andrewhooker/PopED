# @importFrom codetools findGlobals

build_sfg <- function (model="ff.PK.1.comp.oral.sd.CL",covariates=c("DOSE"),par_names=NULL,
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
  
  covariate <- 1:length(parameter_names_ff)*FALSE
  names(covariate) <- parameter_names_ff
  covariate[parameter_names_ff %in% covariates]  <- TRUE
  
  no_eta <- 1:length(parameter_names_ff)*FALSE
  names(no_eta) <- parameter_names_ff
  if(!is.null(no_etas)) no_eta[parameter_names_ff %in% no_etas]  <- TRUE
  
  if(length(etas)==1) etas <- rep(etas,length(parameter_names_ff))
  
  for(k in 1:length(parameter_names_ff)){
    ending <- ", "
    if(k==length(parameter_names_ff)) ending <- ")"
    if(!covariate[k]){ 
      parameters <- paste(parameters,parameter_names_ff[k],"=bpop[",bpop_num,"]",sep="")
      bpop_num <- bpop_num+1
      if(!no_eta[k]){
        if(etas[k]=="exp") parameters <- paste(parameters,"*exp(b[",b_num,"])",sep="")
        if(etas[k]=="prop") parameters <- paste(parameters,"*(1+b[",b_num,"])",sep="")
        if(etas[k]=="add") parameters <- paste(parameters,"+ b[",b_num,"]",sep="")  
        if(etas[k]=="none") parameters <- parameters  
        if(etas[k]!="none") b_num <- b_num+1
      }
      parameters <- paste(parameters,ending,sep="")        
    }
    if(covariate[k]) {
      parameters <- paste(parameters,parameter_names_ff[k],"=a[",a_num,"]",ending,sep="")
      a_num <- a_num+1
    }  
  }
    
  ## add body of funciton and set environment
  body(sfg_tmp) <- parse(text=parameters)
  environment(sfg_tmp) <- env

  return(sfg_tmp)
}