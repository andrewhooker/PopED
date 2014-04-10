feps.add.prop <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- do.call(poped.db$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list( y= y,poped.db =poped.db )) 
}

feps.add <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- do.call(poped.db$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y+epsi[,1]
  
  return(list( y= y,poped.db =poped.db )) 
}

feps.prop <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- do.call(poped.db$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y*(1+epsi[,1])
  
  return(list( y= y,poped.db =poped.db )) 
}
