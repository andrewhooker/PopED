rm(list=ls())
library(PopED)

## Introduction
#Lets assume that we have a model with a covariate included in the model description.  Here we define a one-compartment PK model that has weight on both clearance and volume of distribution. 


mod_1 <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    
    CL=CL*(WT/70)^(WT_CL)
    V=V*(WT/70)^(WT_V)
    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*xt) 
    
    return(list( y= y,poped.db=poped.db))
  })
}

par_1 <- function(x,a,bpop,b,bocc){
  parameters=c( CL=bpop[1]*exp(b[1]),
                V=bpop[2]*exp(b[2]),
                WT_CL=bpop[3],
                WT_V=bpop[4],
                WT=a[1])
  return( parameters ) 
}

#Now we define a design.  In this case one group of individuals, where we define the individuals' typical weight as 70 kg. 

poped_db <- create.poped.database(ff_fun=mod_1,
                                  fg_fun=par_1,
                                  fError_fun=feps.add.prop,
                                  groupsize=50,
                                  m=1,
                                  sigma=c(0.015,0.0015),
                                  notfixed_sigma = c(1,0),
                                  bpop=c(CL=3.8,V=20,WT_CL=0.75,WT_V=1), 
                                  d=c(CL=0.05,V=0.05), 
                                  xt=c( 1,2,4,6,8,24),
                                  minxt=0,
                                  maxxt=24,
                                  bUseGrouped_xt=1,
                                  a=c(WT=70)
)



#We can create a plot of the model typical predictions:

plot_model_prediction(poped_db)



#And evaluate the initial design

evaluate_design(poped_db)

#We see that the covariate parameters can not be estimated according to this design calculation (RSE of bpop[3]=0 and bpop[4]=0).  Why is that? Well, the calculation being done is assuming that every individual in the group has the same covariate (to speed up the calculation).  This is clearly a poor prediction in this case!

## distribution of covariates
#We can improve the computation by assuming a distribution of the covariate (WT) in the individuals in the study. We set `groupsize=1`, the number of groups to be 50 (`m=50`) and assume that WT is sampled from a normal distribution with mean=70 and sd=10 (`a=as.list(rnorm(50,mean = 70,sd=10)`).

poped_db_2 <- create.poped.database(ff_fun=mod_1,
fg_fun=par_1,
fError_fun=feps.add.prop,
groupsize=1,
m=50,
sigma=c(0.015,0.0015),
notfixed_sigma = c(1,0),
bpop=c(CL=3.8,V=20,WT_CL=0.75,WT_V=1), 
d=c(CL=0.05,V=0.05), 
xt=c( 1,2,4,6,8,24),
minxt=0,
maxxt=24,
bUseGrouped_xt=1,
a=as.list(rnorm(50,mean = 70,sd=10))
)


evaluate_design(poped_db_2)


#Here we see that, given this distribution of weights, the covariate effect parameters (bpop[3] and bpop[4]) would be well estimated.

#However, we are only looking at one sample of 50 individuals.  Maybe a better approach is to look at the distribution of RSEs over a number of experiments given the expected weight distribution. 


nsim <- 10
rse_list <- c()
for(i in 1:nsim){
poped_db_tmp <- create.poped.database(ff_fun=mod_1,
fg_fun=par_1,
fError_fun=feps.add.prop,
groupsize=1,
m=50,
sigma=c(0.015,0.0015),
notfixed_sigma = c(1,0),
bpop=c(CL=3.8,V=20,WT_CL=0.75,WT_V=1), 
d=c(CL=0.05,V=0.05), 
xt=c( 1,2,4,6,8,24),
minxt=0,
maxxt=24,
bUseGrouped_xt=1,
a=as.list(rnorm(50,mean = 70,sd=10)))
rse_tmp <- evaluate_design(poped_db_tmp)$rse
rse_list <- rbind(rse_list,rse_tmp)
}
apply(rse_list,2,quantile)



