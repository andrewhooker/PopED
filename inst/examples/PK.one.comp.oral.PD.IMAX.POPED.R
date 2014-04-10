sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                E0=bpop[5]*exp(b[4]),
                IMAX=bpop[6],
                IC50=bpop[7]
  )
  return( parameters ) 
}

# must be in order of bpop numbers
bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,E0=1120,IMAX=0.807,IC50=0.0993)
# must be in order of b numbers
d_vals <- c(V=0.09,KA=0.09,CL=0.25^2,E0=0.09)

notfixed_bpop <- c(1,1,1,0,1,1,1)

ff <- function(model_switch,xt,parameters,globalStructure){
  ##-- Model: One comp first order absorption + inhibitory imax
  ## -- works for both mutiple and single dosing
  y=xt
  MS <- model_switch
  
  with(as.list(parameters),{
    
    # PK model
    tau=24
    N = floor(xt/tau)+1
    KE = CL/V
    CONC=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * tau)) * (1 - exp(-N * KE * tau))/(1 - exp(-KE * tau)) - 
         exp(-KA * (xt - (N - 1) * tau)) * (1 - exp(-N * KA * tau))/(1 - exp(-KA * tau))) 
    
    # PD model
    EFF = E0*(1 - CONC*IMAX/(IC50 + CONC))
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,globalStructure=globalStructure))
  })
}

feps <- function(model_switch,xt,parameters,epsi,globalStructure){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- feval(globalStructure$ff_pointer,model_switch,xt,parameters,globalStructure) 
  y <- returnArgs[[1]]
  globalStructure <- returnArgs[[2]]
  
  MS <- model_switch
  
  pk.dv <- y*(1+epsi[,1])+epsi[,2]
  pd.dv <-  y*(1+epsi[,3])+epsi[,4]
  
  y[MS==1] = pk.dv[MS==1]
  y[MS==2] = pd.dv[MS==2]
  
  return(list( y= y,globalStructure =globalStructure )) 
}

# -- Matrix defining the variances of the residual variability terms --
sigma_vals <- diag(c(0.04,5e-6,0.09,100))
notfixed_sigma <- c(0,0,0,0)


poped.db.1 <- create.poped.database(list(),
                                    ff_file="ff",
                                    fError_file="feps",
                                    groupsize=20,
                                    m=2,
                                    sigma=sigma_vals,
                                    bpop=bpop_vals,  # should be in model definition and be named values
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=notfixed_sigma,
                                    d=d_vals, # should be in model definition and be named values
                                    xt=c( 1,2,8,240,240,1,2,8,240,240),
                                    minxt=c(0,0,0,240,240,0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248,10,10,10,248,248),
                                    G_xt=c(1,2,3,4,5,1,2,3,4,5),
                                    model_switch=c(1,1,1,1,1,2,2,2,2,2),
                                    a=rbind(20,40),
                                    bUseGrouped_xt=1,
                                    maxa=200,# should associated with each xt more tightly
                                    mina=0,# should associated with each xt more tightly
                                    modtit='One comp, oral dose' 
)  

poped.db.2 <- create.poped.database(list(),
                                    ff_file="ff",
                                    fError_file="feps",
                                    groupsize=20,
                                    m=3,
                                    sigma=sigma_vals,
                                    bpop=bpop_vals,  # should be in model definition and be named values
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=notfixed_sigma,
                                    d=d_vals, # should be in model definition and be named values
                                    xt=c( 1,2,8,240,240,1,2,8,240,240),
                                    minxt=c(0,0,0,240,240,0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248,10,10,10,248,248),
                                    G_xt=c(1,2,3,4,5,1,2,3,4,5),
                                    model_switch=c(1,1,1,1,1,2,2,2,2,2),
                                    a=rbind(20,40,0),
                                    bUseGrouped_xt=1,
                                    ourzero=0,
                                    maxa=200,# should associated with each xt more tightly
                                    mina=0,# should associated with each xt more tightly
                                    modtit='One comp, oral dose' 
)  

bpop_vals_ed <- cbind(zeros(7,1),bpop_vals,zeros(7,1)) 
bpop_vals_ed["IC50",1] <- 1
bpop_vals_ed["IC50",3] <- (bpop_vals_ed["IC50",2]*0.1)^2
bpop_vals_ed

poped.db.3 <- create.poped.database(list(),
                                    ff_file="ff",
                                    fError_file="feps",
                                    groupsize=20,
                                    m=3,
                                    sigma=sigma_vals,
                                    bpop=bpop_vals_ed,  # should be in model definition and be named values
                                    notfixed_bpop=notfixed_bpop,
                                    notfixed_sigma=notfixed_sigma,
                                    d=d_vals, # should be in model definition and be named values
                                    xt=c( 1,2,8,240,240,1,2,8,240,240),
                                    minxt=c(0,0,0,240,240,0,0,0,240,240),
                                    maxxt=c(10,10,10,248,248,10,10,10,248,248),
                                    G_xt=c(1,2,3,4,5,1,2,3,4,5),
                                    model_switch=c(1,1,1,1,1,2,2,2,2,2),
                                    a=rbind(20,40,0),
                                    bUseGrouped_xt=1,
                                    ourzero=0,
                                    ED_samp_size=20,
                                    maxa=200,# should associated with each xt more tightly
                                    mina=0,# should associated with each xt more tightly
                                    modtit='One comp, oral dose' 
)  


