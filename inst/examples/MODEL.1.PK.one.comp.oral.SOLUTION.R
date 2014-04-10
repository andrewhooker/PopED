library(PopED)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1])
  return( parameters ) 
}

bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9)
d_vals <- c(V=0.09,KA=0.09,CL=0.25^2)
notfixed_bpop <- c(1,1,1,0)


ff.1.comp.oral.md <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  y=xt
  with(as.list(parameters),{
    tau=24
    N = floor(xt/tau)+1
    KE = CL/V
    y=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * tau)) * (1 - exp(-N * KE * tau))/(1 - exp(-KE * tau)) - 
         exp(-KA * (xt - (N - 1) * tau)) * (1 - exp(-N * KA * tau))/(1 - exp(-KA * tau))) 
    
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
sigma_vals_correct <- diag(c(0.04,5e-6))
notfixed_sigma <- c(0,0)


# initial design
poped.db.1 <- create.poped.database(list(),
                                    ff_file="ff.1.comp.oral.md",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    sigma=sigma_vals,
                                    bpop=bpop_vals,  
                                    d=d_vals, 
                                    xt=c( 1,2,8,16,24),
                                    maxxt=336, 
                                    minxt=0, 
                                    a=20)  

# the design with corrected residual error
poped.db.2 <- create.poped.database(ff_file="ff.1.comp.oral.md",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    sigma=sigma_vals_correct,
                                    bpop=bpop_vals,  
                                    d=d_vals, 
                                    xt=c( 1,2,8,16,24),
                                    maxxt=336, 
                                    minxt=0, 
                                    a=20)  

# the design with fixed parameters
poped.db.3 <- create.poped.database(ff_file="ff.1.comp.oral.md",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    sigma=sigma_vals_correct,
                                    bpop=bpop_vals,  
                                    d=d_vals,
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=notfixed_sigma,
                                    xt=c( 1,2,8,16,24),
                                    maxxt=336, 
                                    minxt=0, 
                                    a=20)  


# the design with dose design space defined
poped.db.4 <- create.poped.database(ff_file="ff.1.comp.oral.md",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    sigma=sigma_vals_correct,
                                    bpop=bpop_vals,  
                                    d=d_vals,
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=notfixed_sigma,
                                    xt=c( 1,2,8,16,24),
                                    maxxt=336, 
                                    minxt=0, 
                                    a=20,
                                    maxa=200,
                                    mina=0)

poped.db.5 <- create.poped.database(list(),
                                    ff_file="ff.1.comp.oral.md",
                                    fError_file="feps.add.prop",
                                    groupsize=20,
                                    m=2,
                                    sigma=sigma_vals_correct,
                                    bpop=bpop_vals, 
                                    d=d_vals, 
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=notfixed_sigma,
                                    xt=c( 1,2,8,240,245),
                                    minxt=c(0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248),
                                    a=rbind(20,40),
                                    bUseGrouped_xt=1,
                                    maxa=200,
                                    mina=0)
