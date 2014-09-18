build_sfg <- function (model="ff.PK.1.comp.oral.sd.CL",covariates=c("DOSE"),env = parent.frame()) {
  
  ## get variable of function
  parameter_names_ff <- codetools::findGlobals(eval(parse(text=model)),merge=F)$variables  
  
  ## create an empty function
  sfg_tmp <- function(x,a,bpop,b,bocc){}
  
  ## get body of function
  parameters <- "parameters=c("
  bpop_num <- 1
  b_num <- 1
  a_num <- 1
  covariate <- 1:length(parameter_names_ff)*FALSE
  names(covariate) <- parameter_names_ff
  covariate[parameter_names_ff %in% covariates]  <- 1
  for(k in 1:length(parameter_names_ff)){
    ending <- ", "
    if(k==length(parameter_names_ff)) ending <- ")"
    if(!covariate[k]){ 
      parameters <- paste(parameters,parameter_names_ff[k],"=bpop[",bpop_num,"]*exp(b[",b_num,"])",ending,sep="")
      bpop_num <- bpop_num+1
      b_num <- b_num+1
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