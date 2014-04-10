library(PopED)
library(deSolve)

#defintion of parametes
#defintion of parametes
sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  return( parameters ) 
}

#definition of values?
bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
d_vals <- c(V=0.09,KA=0.09,CL=0.25^2)


###########################################################
#definition of model using odeSolve
#CODE for using odeSOLVE ipv y=closed form solution 
###########################################################

###############
##Single Dose
###############
ff.ODE.SD <- function(model_switch, xt, parameters, globalStructure){
  ##--Model: one comp, linear absorptin, single dose
  ##parametrized by CL, KA, F and V
  with(as.list(parameters),{
    A_ini <- c(A1=DOSE, A2=0)
    times <- drop(xt) #put in solve function
    times <- sort(times) #put in solve function
    times <- c(0,times) #put in solve function: add extra time for start of itegration
    out <- ode(A_ini, times, one.comp.ode, parameters) #atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V/Favail)
    y = y[-1] #remove initial time for start of integaration
    ##use order to sort output based on order of times or use extraction to get right time and y corelation
    y=cbind(y)##must be a matrix with columns as group results
    return(list(y=y,globalStructure=globalStructure))
    
  })
}


###############
##Multiple Dose
###############
ff.ODE.MD <- function(model_switch, xt, parameters, globalStructure){
  
  
  ##--Model: one comp, linear absorptin, multiple doses
  ##parametrized by CL, KA, F and V
  with(as.list(parameters),{
    
 
    
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt)#xt[,,drop=T] #put in solve function
    
    ## multiple dosing
    tau=24    #tau is dosing interval  
    dose_times = seq(from=0,to=max(times_xt),by=tau)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE), method = c("add"))
    
    times <- sort(c(times_xt,dose_times))
    
    
    out <- ode(A_ini, times, one.comp.ode, parameters, events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V/Favail)
    
    y=y[match(times_xt,out[,"time"])]
    
    y=cbind(y)
    return(list(y=y,globalStructure=globalStructure))
    
  })
}


one.comp.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    KE <- CL/V
    
    dA1 <- -KA*A1
    dA2 <- KA*A1 - KE*A2
    
    return(list(c(dA1, dA2)))
  })
}




###########################################################
#definition of model using Closed Form Solution (CSF)
###########################################################

###############
##Single Dose
###############
ff.1.comp.oral.SD <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  y=xt
  with(as.list(parameters),{
    tau=250    #tau is dosing interval
    N = floor(xt/tau)+1 #N determines the number of doses
    KE = CL/V
    y=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * tau)) * (1 - exp(-N * KE * tau))/(1 - exp(-KE * tau)) - 
         exp(-KA * (xt - (N - 1) * tau)) * (1 - exp(-N * KA * tau))/(1 - exp(-KA * tau))) 
    
    
    return(list( y= y,globalStructure=globalStructure))
  })
}

###############
##Multiple Dose
###############
ff.1.comp.oral.MD <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  y=xt
  with(as.list(parameters),{
    tau=24    #tau is dosing interval
    N = floor(xt/tau)+1 #N determines the number of doses
    
    
    KE = CL/V
    y=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * tau)) * (1 - exp(-N * KE * tau))/(1 - exp(-KE * tau)) - 
         exp(-KA * (xt - (N - 1) * tau)) * (1 - exp(-N * KA * tau))/(1 - exp(-KA * tau))) 
    
    
    return(list( y= y,globalStructure=globalStructure))
  })
}






###########################################################
#definition of residual error structure
###########################################################

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
#definition of residual error values
sigma_vals <- diag(c(0.04,0.000005))



#########################
#Question 4 Part II: MD
#########################

poped.db.4.odeSolve.MD <- create.poped.database(list(),
                                                ff_file="ff.ODE.MD",
                                                fError_file="feps.add.prop",
                                                groupsize=20,
                                                m=2,      #number of groups
                                                sigma=sigma_vals,
                                                bpop=bpop_vals,  
                                                d=d_vals, 
                                                xt=c( 1,2,8, 240,240),
                                                minxt=c(0,0,0,240,240),
                                                maxxt=c(10,10,10,248,248),
                                                bUseGrouped_xt=1,#assume both groups have same sample times (useful for pcb to avoid unblinding)
                                                a=rbind(20,40),#two different doses, each group one
                                                notfixed_bpop=c(1,1,1,0),
                                                notfixed_d=c(1,1,1),
                                                notfixed_sigma=c(1,1),
                                                mina=0,   #can make vector here eg c(0,5)
                                                maxa=100)


poped.db.4.CFS.MD <- create.poped.database(list(),
                                           ff_file="ff.1.comp.oral.MD",
                                           fError_file="feps.add.prop",
                                           groupsize=20,
                                           m=2,      #number of groups
                                           sigma=sigma_vals,
                                           bpop=bpop_vals,  
                                           d=d_vals, 
                                           xt=c( 1,2,8, 240,240),
                                           minxt=c(0,0,0,240,240),
                                           maxxt=c(10,10,10,248,248),
                                           bUseGrouped_xt=1,#assume both groups have same sample times (useful for pcb to avoid unblinding)
                                           a=rbind(20,40),#two different doses, each group one
                                           notfixed_bpop=c(1,1,1,0),
                                           notfixed_d=c(1,1,1),
                                           notfixed_sigma=c(1,1),
                                           mina=0,   #can make vector here eg c(0,5)
                                           maxa=100)


plot_model_prediction(poped.db.4.odeSolve.MD,IPRED=T,DV=T,separate.groups=T)
plot_model_prediction(poped.db.4.CFS.MD,IPRED=T,DV=T,separate.groups=T)



#########################################
#I do not get MD in de deSolve setting
#So for time being simplify to SD
#########################################





#########################
#Question 4 Part II: SD
#########################

poped.db.4.odeSolve.SD <- create.poped.database(list(),
                                                ff_file="ff.ODE.SD",
                                                fError_file="feps.add.prop",
                                                groupsize=20,
                                                m=2,      #number of groups
                                                sigma=sigma_vals,
                                                bpop=bpop_vals,  
                                                d=d_vals, 
                                                xt=c( 1,2,8, 12,50),
                                                minxt=0,
                                                maxxt=50,
                                                bUseGrouped_xt=1,#assume both groups have same sample times (useful for pcb to avoid unblinding)
                                                a=rbind(20,40),#two different doses, each group one
                                                notfixed_bpop=c(1,1,1,0),
                                                notfixed_d=c(1,1,1),
                                                notfixed_sigma=c(1,1),
                                                mina=0,   #can make vector here eg c(0,5)
                                                maxa=100)


poped.db.4.CFS.SD <- create.poped.database(list(),
                                           ff_file="ff.1.comp.oral.SD",
                                           fError_file="feps.add.prop",
                                           groupsize=20,
                                           m=2,      #number of groups
                                           sigma=sigma_vals,
                                           bpop=bpop_vals,  
                                           d=d_vals, 
                                           xt=c( 1,2,8, 12,50),
                                           minxt=0,
                                           maxxt=50,
                                           bUseGrouped_xt=1,#assume both groups have same sample times (useful for pcb to avoid unblinding)
                                           a=rbind(20,40),#two different doses, each group one
                                           notfixed_bpop=c(1,1,1,0),
                                           notfixed_d=c(1,1,1),
                                           notfixed_sigma=c(1,1),
                                           mina=0,   #can make vector here eg c(0,5)
                                           maxa=100)


plot_model_prediction(poped.db.4.odeSolve.SD,IPRED=T,DV=T,separate.groups=T)
plot_model_prediction(poped.db.4.CFS.SD,IPRED=T,DV=T,separate.groups=T)



FIM.4.CSF.SD <- evaluate.fim(poped.db.4.CFS.SD)
FIM.4.odeSolve.SD <- evaluate.fim(poped.db.4.odeSolve.SD)

det(FIM.4.CSF.SD)
det(FIM.4.odeSolve.SD)

get_rse(FIM.4.CSF.SD, poped.db.4.CFS.SD)
get_rse(FIM.4.odeSolve.SD, poped.db.4.odeSolve.SD)

