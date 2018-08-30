
  ## following example 1.4.2 in PFIM
  sfg <- function(x,a,bpop,b,bocc){
    ## -- parameter definition function 
    parameters=c(KA=bpop[1]*exp(b[1]),
                 K=bpop[2]*exp(b[2]),
                 V=bpop[3]*exp(b[3]),
                 DOSE=a[1])
    return(parameters) 
  }
  
  ff <- function(model_switch,xt,parameters,poped.db){
    ##-- Model: One comp first order absorption
    with(as.list(parameters),{
      y<-(DOSE/V*KA/(KA-K)*(exp(-K*xt)-exp(-KA*xt)))
      return(list(y=y,poped.db=poped.db))
    })
  }
  
  ## -- Residual unexplained variablity (RUV) function
  ## -- Additive + Proportional  
  feps <- function(model_switch,xt,parameters,epsi,poped.db){
    returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
    f <- returnArgs[[1]]
    y = f + (0.5 + 0.15*f)*epsi[,1]
    
    return(list( y= y,poped.db =poped.db )) 
  }
  
  ## -- Define initial design  and design space
  poped.db <- create.poped.database(ff_fun=ff,
                                    fg_fun=sfg,
                                    fError_fun =feps,
                                    bpop=c(KA=2, K=0.25, V=15), 
                                    d=c(KA=1, V=0.25,0.1), 
                                    sigma=c(1),
                                    notfixed_sigma = c(0),
                                    groupsize=1,
                                    xt=c( 1,3,8),
                                    minxt=0,
                                    maxxt=10,
                                    a=100)
  
  
  shr_out <- shrinkage(poped.db)

  

