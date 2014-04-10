library(PopED)

sfg.1.comp.oral.md.imax <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU = a[2],
                E0=bpop[5]*exp(b[4]),
                IMAX=bpop[6],
                IC50=bpop[7])
  parameters["KE"]=parameters["CL"]/parameters["V"]
  return( parameters ) 
}

# must be in order of bpop numbers
bpop_vals <- c(V=72.8,KA=0.25,CL=3.75,Favail=0.9,E0=1120,IMAX=0.807,IC50=0.0993)
# must be in order of b numbers
d_vals <- c(V=0.09,KA=0.09,CL=0.25^2,E0=0.09)

notfixed_bpop <- c(1,1,1,0,1,1,1)

ff.1.comp.oral.md.imax <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption + inhibitory imax
  ## -- works for both mutiple and single dosing  
  with(as.list(parameters),{
    
    y=xt
    MS <- model_switch
    
    # PK model
    returnArgs=ff.1.comp.oral.md(model_switch,xt,parameters,poped.db)
    CONC=returnArgs$y
    
    # PD model
    EFF = E0*(1 - CONC*IMAX/(IC50 + CONC))
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

feps.1.comp.oral.md.imax <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  returnArgs <- feval(poped.db$ff_pointer,model_switch,xt,parameters,poped.db) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  MS <- model_switch
  
  pk.dv <- y*(1+epsi[,1])+epsi[,2]
  pd.dv <-  y*(1+epsi[,3])+epsi[,4]
  
  y[MS==1] = pk.dv[MS==1]
  y[MS==2] = pd.dv[MS==2]
  
  return(list( y= y,poped.db =poped.db )) 
}

# -- Matrix defining the variances of the residual variability terms --
sigma_vals <- diag(c(0.04,5e-6,0.09,100))
notfixed_sigma <- c(0,0,0,0)

poped.db <- create.poped.database(ff_file="ff.1.comp.oral.md.imax",
                                  fError_file="feps.1.comp.oral.md.imax",
                                  fg_file="sfg.1.comp.oral.md.imax",
                                  groupsize=20,
                                  m=3,
                                  sigma=sigma_vals,
                                  bpop=bpop_vals,  
                                  notfixed_bpop=notfixed_bpop,
                                  notfixed_sigma=notfixed_sigma,
                                  d=d_vals, 
                                  xt=c( 1,2,8,240,240,1,2,8,240,240),
                                  minxt=c(0,0,0,240,240,0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248,10,10,10,248,248),
                                  G_xt=c(1,2,3,4,5,1,2,3,4,5),
                                  model_switch=c(1,1,1,1,1,2,2,2,2,2),
                                  a=cbind(c(20,40,0),c(24,24,24)),
                                  bUseGrouped_xt=1,
                                  ourzero=0,
                                  maxa=c(200,40),
                                  mina=c(0,2))

print(plot_model_prediction(poped.db,facet_scales="free"))
print(plot_model_prediction(poped.db,IPRED=T,DV=T,facet_scales="free",separate.groups=T))

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
print(FIM)
print(det(FIM))
print(get_rse(FIM,poped.db))

