library(PopED)

PK.one.cpt.oral.params.trans2 <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}

bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
d_vals <- c(V=0.09,KA=0.09,CL=0.25^2)

PK.one.cpt.oral.ff.model <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  y=xt
  with(as.list(parameters),{
    N = floor(xt/TAU)+1
    KE = CL/V
    y=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * TAU)) * (1 - exp(-N * KE * TAU))/(1 - exp(-KE * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU))) 
    
    return(list( y= y,globalStructure=globalStructure))
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
sigma_vals <- diag(c(0.04,5))

poped.db.1 <- create.poped.database(list(),
                                    fg_file="PK.one.cpt.oral.params.trans2",
                                    ff_file="PK.one.cpt.oral.ff.model",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    sigma=sigma_vals,
                                    bpop=bpop_vals,  
                                    d=d_vals, 
                                    xt=c( 1,2,8,16,24),
                                    maxxt=336, 
                                    minxt=0, 
                                    a=c(20,24))  

