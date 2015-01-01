popedInput <- function_input_fo_reduced_eval_pfim_way()
poped.db <- convert_popedInput(popedInput)
poped.db <- convert_variables(poped.db)
downsize.list <- downsizing_general_design(poped.db)
tmp.names <- c("ni", "xt", "model_switch", "x", "a", "bpop", "n", "d", "maxxt", "minxt","maxa","mina")
eval(parse(text=paste(tmp.names,"<-","downsize.list$",tmp.names,sep="")))
bpop_val=bpop[,2,drop=F]
d_val=getfulld(d[,2,drop=F],poped.db$parameters$covd)
docc_full = getfulld(poped.db$parameters$docc[,2],poped.db$parameters$covdocc)
groupsize=poped.db$design$groupsize
sigma=poped.db$parameters$sigma
output = mftot4(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db)
det(output$ret)

## about 0.015sec per computation
## in matlab it is 0.007sec per computation (2x faster)
## PFIM takes 0.25, but that is with graphics.  but anyway it is whithin not faster
tot.time <- 0
n <- 100
for(i in 1:n){
  times <- system.time(mftot4(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db))
  tot.time <- tot.time + times
}
avg.time <- tot.time/n
avg.time

parse_rprof(getwd())

Rprof()
mftot0(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db)
Rprof(NULL)
summaryRprof()

## tot.time <- 0
## n <- 100
## for(i in 1:n){
##     start.time <- proc.time()
##     mftot0(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db)
##     tot.time <- tot.time + proc.time() - start.time
## }
## avg.time <- tot.time/n
## avg.time

#     tot.time <- 0
#     n <- 100
#     for(i in 1:n){
#       times <- system.time(mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db))
#       tot.time <- tot.time + times
#     }
#     avg.time <- tot.time/n
#     avg.time
#     
#     system.time(for(i in 1:100) mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db))
#     system.time(mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db))
#     
#     system.time(for(i in 1:100){
#     evaluate.fim(poped.db,
#                  bpop.val=bpop,
#                  d_full=d,
#                  docc_full=docc_full,
#                  sigma_full=poped.db$parameters$sigma,
#                  model_switch=model_switch,
#                  ni=ni,
#                  xt=xt,
#                  x=x,
#                  a=a,
#                  groupsize=poped.db$design$groupsize,
#                  ...)

library(compiler)
enableJIT(1)
#enableJIT(0)
tot.time <- 0
n <- 100
for(i in 1:n){
  times <- system.time(mftot0(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db))
  tot.time <- tot.time + times
}
avg.time <- tot.time/n
avg.time

## #enableJIT(0)
tot.time <- 0
n <- 100
for(i in 1:n){
    start.time <- proc.time()
    mftot0(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db)
    tot.time <- tot.time + proc.time() - start.time
}
avg.time <- tot.time/n
avg.time

## enableJIT(0)
## library(compiler)
## mftot0.cmp <- cmpfun(mftot0)
## tot.time <- 0
## n <- 100
## for(i in 1:n){
##     start.time <- proc.time()
##     mftot0.cmp(model_switch,groupsize,ni,xt,x,a,bpop_val,d_val,sigma,docc_full,poped.db)
##     tot.time <- tot.time + proc.time() - start.time
## }
## avg.time <- tot.time/n
## avg.time
