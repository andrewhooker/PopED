library(PopED)

PK.1.comp.maturation.fg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( CL=bpop[1]*exp(b[1]),
                V=bpop[2]*exp(b[2]),
                EMAX=bpop[3],
                EC50=bpop[4],
                HILL=bpop[5],
                WT=a[1])
  return( parameters ) 
}

bpop_vals <- c(CL=1.8,V=20,EMAX=2,EC50=25,HILL=5)
d_vals <- c(CL=0.05,V=0.05)

PK.1.comp.maturation.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    
    CL=CL+(EMAX*WT**HILL)/(EC50**HILL+WT**HILL)
    V=V*(WT/70)
    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*xt) 
       
    return(list( y= y,poped.db=poped.db))
  })
}

feps.add.prop <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- feval(poped.db$ff_pointer,model_switch,xt,parameters,poped.db) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list( y= y,poped.db =poped.db )) 
}

# -- Matrix defining the variances of the residual variability terms --
sigma_vals <- diag(c(0.015,0.0015))

# initial design
poped.db.1 <- create.poped.database(ff_file="PK.1.comp.maturation.ff",
                                    fError_file="feps.add.prop",
                                    fg_file="PK.1.comp.maturation.fg",
                                    groupsize=rbind(50,20,20,20),
                                    m=4,
                                    sigma=sigma_vals,
                                    bpop=bpop_vals, 
                                    d=d_vals, 
                                    xt=c( 1,2,4,6,8,24),
                                    minxt=0,
                                    maxxt=24,
                                    a=rbind(70,60,50,10),
                                    maxa=70,
                                    mina=1)


## Choose database
poped.db <- poped.db.1 # original model and design

##  create plot of model 
plot_model_prediction(poped.db)
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)


## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

# RS+SG+LS optimization of sample times
output <- poped_optimize(poped.db,opt_xt=T,opt_a=T)

# MFEA optimization with only integer times allowed
 mfea.output <- poped_optimize(poped.db,opt_a=1,opt_xt=0,
                               bUseExchangeAlgorithm=1,
                               EAStepSize=1)
plot_model_prediction(mfea.output$poped.db,IPRED=T,DV=T,separate.groups=F)

