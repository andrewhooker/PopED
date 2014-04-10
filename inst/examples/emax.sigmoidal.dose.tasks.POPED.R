## libraries
library(deSolve)
library(PopED)
library(mvtnorm)
#source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/tools/source.poped.code.R")

#--- load model functions
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/emax.sigmoidal.dose.model.POPED.R")

#--- load design functions
source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/emax.sigmoidal.dose.design.1.POPED.R") 

# set working directory
setwd("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/")

## choose a model and design 
poped.db <- create.poped.database(emax.sigmoidal.dose.design.1.input())

##  create plot of model 
plot_model_prediction(poped.db,IPRED=T,DV=T)

## get pred from model
pred <- model_prediction(poped.db)
pred

## evaluate initial design
FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
det(FIM)
get_rse(FIM,poped.db)

## check all FIM calculations
df <- c()
poped.db.1 <- poped.db
for(k in c(1)){
  tmp.db <- eval(parse(text=paste("poped.db",k,sep=".")))
  for(i in c(0,1,4,5,6,7)){
    for(j in c(0,1)){
      if(k==4 && j==0) next
      FIM <- evaluate.fim(tmp.db,fim.calc.type=i,deriv.type=j) # new name for function needed
      tmp <- data.frame("Model"= k,"fim.calc.type"= i, "deriv.type"=j, "FIM"=det(FIM))
      if(size(df,2)==0) cat(paste("Model","fim.calc.type","deriv.type","FIM\n",sep=" "))
      cat(paste(k,i,j,det(FIM),"\n",sep=" "))      
      df <- rbind(df,tmp)
    }
  }
}
print(df,digits=3,row.names=F,print.gap=3)

# optimization -----------------------
rs.output <- RS_opt(poped.db,opt_xt=1,rsit=20)
rs.output$xt
names(rs.output)
result.db <- rs.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)

rs.output <- RS_opt(poped.db,opt_xt=1,rsit=20,ofv.calc.type=0) 

if(FALSE){ ## long optimization
  rs.output <- RS_opt(poped.db,opt_xt=1,rsit=300,ofv.calc.type=4) 
}

# using doptim with RS
rs.output <- poped_optimize(poped.db,opt_xt=1,rsit=20,
                            bUseRandomSearch= 1,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0)
names(rs.output)
get_rse(rs.output$fmf,rs.output$poped.db)
result.db <- rs.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)

# line search
ls.output <- poped_optimize(poped.db,opt_xt=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 1,
                            ls_step_size=10)
get_rse(ls.output$fmf,ls.output$poped.db)
result.db <- ls.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)

# Stochastic gradient search
sg.output <- poped_optimize(poped.db,opt_xt=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 1,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            sgit=20)
get_rse(sg.output$fmf,sg.output$poped.db)
result.db <- sg.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)

# BFGS search (need to add output to result file)
bfgs.output <- poped_optimize(poped.db,opt_xt=1,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 1,bUseLineSearch = 0)

get_rse(bfgs.output$fmf,bfgs.output$poped.db)
result.db <- bfgs.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)


# RS+SG+LS
sg.output <- poped_optimize(poped.db,opt_xt=1,
                            bUseRandomSearch= 1,bUseStochasticGradient = 1,bUseBFGSMinimizer = 1,bUseLineSearch = 1,
                            sgit=10,rsit=10,ls_step_size=10,ofv_calc_type=1,fim.calc.type=0)
get_rse(sg.output$fmf,sg.output$poped.db)
result.db <- sg.output$poped.db
plot_model_prediction(result.db,IPRED=T,DV=T)



# MFEA search does not work for continuous optimization
mfea.output <- poped_optimize(poped.db,opt_xt=1,opt_a=0,
                            bUseRandomSearch= 0,bUseStochasticGradient = 0,bUseBFGSMinimizer = 0,bUseLineSearch = 0,
                            bUseExchangeAlgorithm=1)



## ED evaluate. with MC and random sampling
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
output <- evaluate.e.ofv.fim(poped.db,use.laplace=T)
output$E_ofv
output$E_fim

## ds evaluate
poped.db$ds_index <- cbind(0,0,0,1,1,1,1) # size is number_of_non_fixed_parameters 
## should be the size of num_of_parameters in a future version
poped.db$ofv_calc_type <- 6

FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
det(FIM)
ofv_fim(FIM,poped.db,ofv_calc_type=6,ds_index=cbind(0,0,0,1,1,1,1))

rs.output <- RS_opt(poped.db,opt_xt=1,rsit=20,fim.calc.type=4,ofv_calc_type=6) 

# ## plot model BSV and RV
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

## git

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



