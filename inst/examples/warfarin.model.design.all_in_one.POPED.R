sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  return( parameters ) 
}

sfg.discrete <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=x[1])
  return( parameters ) 
}


sfg_vals <- function(){
  bpop <- rbind(c(4,8,4),
             c(4,1,0.0625),
             c(4,0.15,0.0014),
             c(0,1,0))
  
  notfixed_bpop=cbind(1,1,1,0)
  
  d=rbind(c(4,0.02,0.000025),
          c(4,0.6,0.000225),
          c(4,0.07,0.00030625))
  
  return(list(bpop=bpop,notfixed_bpop=notfixed_bpop,d=d))
}

param_test <- function(){
  V_pop = list(bpop[1])
}


ff <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp first order absorption
  ##-- Description: Parameterized by V, ka,CL and F. 
  y=xt
  with(as.list(parameters),{
    
    KE = CL/V
    y=(DOSE*Favail*KA*KE/(CL*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
    
    return(list( y= y,globalStructure=globalStructure))
  })
}

ff.ODE <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE, A2 = 0)
    times <- drop(xt)##xt[,,drop=T] ## put in solve function
    times <- sort(times) ## put in solve function
    times <- c(0,times) ## put in solve function: add extra time for start of integration
    out   <- ode(A_ini, times, one.comp.ode, parameters)#,atol=1e-13,rtol=1e-13)
    y = out[,"A2"]/(V/Favail)
    y=y[-1] # remove initial time for start of integration
    ## use order to sort output based on order of times or use extraction to get right time and y correlation
    y = cbind(y) ## must be a matrix with columns as group results
    return(list( y= y,globalStructure=globalStructure)) 
  })
}

one.comp.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    KE <- CL/V
    
    dA1  <- -KA*A1
    dA2  <- KA*A1 - KE*A2
    
    return(list(c(dA1, dA2)))
  })
}


feps.add.prop <- function(model_switch,xt,parameters,epsi,globalStructure){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- feval(globalStructure$ff_pointer,model_switch,xt,parameters,globalStructure) 
  y <- returnArgs[[1]]
  globalStructure <- returnArgs[[2]]
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list( y= y,globalStructure =globalStructure )) 
}

# -- Matrix defining the variances of the residual variability terms --
sigma_vals <- function(){
  diag(c(0.01,0.25))
} 

feps.add <- function(model_switch,xt,parameters,epsi,globalStructure){
  ## -- Residual Error function
  ## -- Additive  
  returnArgs <- feval(globalStructure$ff_pointer,model_switch,xt,parameters,globalStructure) 
  y <- returnArgs[[1]]
  globalStructure <- returnArgs[[2]]
  y = y+epsi[,1]
  
  return(list( y= y,globalStructure =globalStructure )) 
}

feps.prop <- function(model_switch,xt,parameters,epsi,globalStructure){
  ## -- Residual Error function
  ## -- Proportional 
  returnArgs <- feval(globalStructure$ff_pointer,model_switch,xt,parameters,globalStructure) 
  y <- returnArgs[[1]]
  globalStructure <- returnArgs[[2]]
  y = y*(1+epsi[,1])
  
  return(list( y= y,globalStructure =globalStructure )) 
}


warfarin.design.1.red.input <- function(){
  
  
  popedInput <- list()
  ## --------------------------
  ## ---- Labeling and file names
  ## --------------------------
  
  ## -- Filname and path of the model file --
  popedInput$ff_file='ff'
  ## -- Filname and path of the parameter file --
  popedInput$fg_file='sfg'
  ## -- Filname and path of the residual error model file --
  popedInput$fError_file='feps.add.prop'
  ## -- The model title --
  popedInput$modtit='One comp first order absorption'
  

  ## --------------------------
  ## ---- Initial Design 
  ## --------------------------
    
  ## -- Vector defining the size of the different groups (num individuals in each group) --
  #   popedInput$design$groupsize=rbind(8,8,8,8)
  popedInput$design$groupsize=8
  
  ## -- Matrix defining the initial sampling schedule --
  #     popedInput$design$xt=rbind(c(0.5,1,2,6,24,36,72,120),
  #                                c(0.5,1,2,6,24,36,72,120),
  #                                c(0.5,1,2,6,24,36,72,120),
  #                                c(0.5,1,2,6,24,36,72,120))
  #  
  popedInput$m=4
  
  popedInput$design$xt=c(0.5,1,2,6,24,36,72,120)
  
  ## -- Vector defining the initial covariate values --
  #   popedInput$design$a=rbind(70, 70, 70, 70)
  popedInput$design$a=70
  
  ## --------------------------
  ## ---- Design space
  ## --------------------------
  
  ## -- Matrix defining the max value for each sample --
  #   popedInput$design$maxxt=rbind(c(25,25,25,120,120,120,120,120),
  #                                 c(25,25,25,120,120,120,120,120),
  #                                 c(25,25,25,120,120,120,120,120),
  #                                 c(25,25,25,120,120,120,120,120))
  popedInput$design$maxxt=c(25,25,25,120,120,120,120,120)
                                
  ## -- Matrix defining the min value for each sample --
  popedInput$design$minxt=0*popedInput$design$xt
  ## -- Vector defining the max value for each covariate --
  #   popedInput$design$maxa=rbind(100, 100, 100, 100)
  popedInput$design$maxa=100
  
  ## -- Vector defining the min value for each covariate --
  #   popedInput$design$mina=rbind(1,1,1,1)
  popedInput$design$mina=1
  
  ## --------------------------
  ## ---- Model parameters and parmeter distributions and fixed or not
  ## --------------------------
  
  params <- sfg_vals()
  ## -- Matrix defining the fixed effects, per row (row number = parameter_number),
  popedInput$design$bpop=params$bpop
  
  ## -- Matrix defining the diagnonals of the IIV (same logic as for the fixed efects) --
  popedInput$design$d=params$d

  ## -- Matrix defining the variances of the residual variability terms --
  popedInput$design$sigma <- sigma_vals()

  ## -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_bpop=params$notfixed_bpop

  
  return( popedInput ) 
}
