## libraries
library(deSolve)
library(PopED)
library(mvtnorm)

source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/tools/source.poped.code.R")

#--- load model functions
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.POPED.R")

#--- load design functions
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.design.1.POPED.R") # 4-groups, add+prop, as pfim
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.design.2.POPED.R") # 1-group, add+prop, as pfim
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.design.3.POPED.R") # 1-group, prop, as pfim
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.design.4.POPED.R") # 1-group, add+prop, as pfim, ODE solution
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.design.5.POPED.R") # 4-group, add+prop, as pfim, ODE solution
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.design.1.red.POPED.R") # 4-group, add+prop, as pfim
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.design.all_in_one.POPED.R") # 4-group, add+prop, as pfim

# set working directory
setwd("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/")

## choose a model and design 
poped.db.2 <- create.poped.database(warfarin.design.2.input())
poped.db.3 <- create.poped.database(warfarin.design.3.input())
poped.db.1 <- create.poped.database(warfarin.design.1.input())
poped.db.4 <- create.poped.database(warfarin.design.4.input())
poped.db.5 <- create.poped.database(warfarin.design.5.input())

poped.db.6 <- create.poped.database(warfarin.design.1.red.input())
if(FALSE){
  poped.db.6 <- create.poped.database(warfarin.design.1.red.input(),bParallelMFEA=TRUE)
  poped.db.6$parallelSettings$bParallelMFEA
}
poped.db <- poped.db.6


poped.db <- poped.db.5
poped.db <- poped.db.4
poped.db <- poped.db.3
poped.db <- poped.db.2
poped.db <- poped.db.1

##  create plot of model with uniform grid of points
if(FALSE){
  poped.db$ga <- rbind(70,60,50,80)
}

plot_model_prediction(poped.db)
plot_model_prediction(poped.db,sample.times=F)
plot_model_prediction(poped.db,separate.groups=T)
plot_model_prediction(poped.db,IPRED=T,DV=T)
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)

## get pred from model
pred <- model_prediction(poped.db)
pred

## evaluate initial design
FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
det(FIM)

## check all FIM calculations
df <- c()
for(k in c(1,2,3,4,5,6)){
  tmp.db <- eval(parse(text=paste("poped.db",k,sep=".")))
  for(i in c(0,1,4,5,6,7)){
    for(j in c(0,1)){
      if(k==4 && j==0) next
      FIM <- evaluate.fim(tmp.db,fim.calc.type=i,deriv.type=j) # new name for function needed
      tmp <- data.frame("Model"= k,"calc.type"= i, "deriv.type"=j, "FIM"=det(FIM))
      if(size(df,2)==0) cat(paste("Model","calc.type","deriv.type","FIM\n",sep=" "))
      cat(paste(k,i,j,det(FIM),"\n",sep=" "))      
      df <- rbind(df,tmp)
    }
  }
}
print(df,digits=3,row.names=F,print.gap=3)

## difference between ODEs and not (more with full)
(subset(df,Model==2 & deriv.type==1)["FIM"] - subset(df,Model==4 & deriv.type==1)["FIM"])/subset(df,Model==2 & deriv.type==1)["FIM"]*100

# optimization -----------------------
rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=20)
rs.output$xt
names(rs.output)
result.db <- rs.output$poped.db
plot_model_prediction(result.db)

rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=20,fim.calc.type=0) 

if(FALSE){ ## long optimization
  rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=300,fim.calc.type=4) 
}

# using doptim with RS
rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                            bUseRandomSearch= 1,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0)
names(rs.output)
get_rse(rs.output$fmf,rs.output$poped.db)

# line search
ls.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 1,
                            ls_step_size=10)

# Stochastic gradient search
sg.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1, 
                            bUseRandomSearch= 0,bUseStochasticGradient = 1,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            sgit=20)

# BFGS search (need to add output to result file)
bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 1,bUseLineSearch = 0)
bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=0,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 1,bUseLineSearch = 0)

# MFEA search does not work for continuous optimization
mfea.output <- poped_optimize(poped.db,opt_xt=0,opt_a=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            bUseExchangeAlgorithm=1)

## discrete dose
poped.db <- poped.db.5
rs.output <- RS_opt(poped.db,opt_xt=0,opt_x=1,opt_a=0,rsit=20)

rs.output <- poped_optimize(poped.db,opt_xt=0,opt_a=0,opt_x=1,rsit=20,
                            bUseRandomSearch= 1,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            iter_tot=1)

ls.output <- poped_optimize(poped.db,opt_xt=0,opt_a=0, opt_x=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 1,
                            ls_step_size=10)

mfea.output <- poped_optimize(poped.db,opt_xt=0,opt_a=0,opt_x=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            bUseExchangeAlgorithm=1)


## ED evaluate. with MC and random sampling
poped.db <- poped.db.2
output <- evaluate.e.ofv.fim(poped.db)
output$E_ofv
output$E_fim

output <- evaluate.e.ofv.fim(poped.db,bLHS=1)
output$E_ofv
output$E_fim

## ED evaluate. with Laplace
################################################################################################
################################################################################################
################################################################################################
## not working! 
################################################################################################
################################################################################################
################################################################################################
poped.db <- poped.db.2
output <- evaluate.e.ofv.fim(poped.db,use.laplace=T)
output$E_ofv
output$E_fim

## ds evaluate
poped.db$ds_index <- cbind(0,0,0,1,1,1,1,1) # size is number_of_non_fixed_parameters 
## should be the size of num_of_parameters in a future version
poped.db$ofv_calc_type <- 6

FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
det(FIM)
ofv_fim(FIM,poped.db,ofv_calc_type=6,ds_index=cbind(0,0,0,1,1,1,1,1))

rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=20,fim.calc.type=4) 

# ## plot model BSV and RV
if(FALSE){
  tmp.db <- poped.db.1
  tmp.db$ga <- rbind(70,60,50,80)
  poped.db <- tmp.db
}

ipred <- model_prediction(poped.db,IPRED=T)
ipred

plot_model_prediction(poped.db,IPRED=F)
plot_model_prediction(poped.db,IPRED=T)

plot_model_prediction(poped.db,sample.times=F,IPRED=T)
plot_model_prediction(poped.db,separate.groups=T,IPRED=T)

dv <- model_prediction(poped.db,DV=T,model_num_points=100)
dv

plot_model_prediction(poped.db,separate.groups=T,IPRED=T,DV=T)
plot_model_prediction(poped.db,separate.groups=F,IPRED=F,DV=T)

plot_model_prediction(poped.db,separate.groups=T,IPRED=T,DV=T,sample.times.IPRED=T,sample.times=F,PRED=F)

## models in libary

## parameter estimation

## vignette

## combine all R functions into one location

## git repository

## make examples 

## plotting functions (during optimization?)

## make more simple input

## split the various optimization tools out and clean up output (both results files and return values)

## compiled code

## examoples with C++ models in differential equations?

## france stuff

## optimize over number of samples (Dpoptim)

## optimize over number of individuals in groups (DpopIndGrp) 

## ED optimize.
##[xt,x,a,fmf,dmf,globalStructure]=EDoptim(globalStructure,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf)

## ED evaluate. with BFGS

## ed optimize with laplace
## ED evaluate. user specified distribution

# FOCE



