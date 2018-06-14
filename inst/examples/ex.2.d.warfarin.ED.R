## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Evaluating with uncertainty around parameter values in the model

library(PopED)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list(y=y,poped.db=poped.db))
  })
}

# Adding 10% Uncertainty to all fixed effects (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                      bpop_vals,
                      ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed["Favail",] <- c(0,1,0)
bpop_vals_ed


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff,
                                  fg_fun=sfg,
                                  fError_fun=feps.add.prop,
                                  bpop=bpop_vals_ed, 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70,
                                  mina=0,
                                  maxa=100,
                                  ED_samp_size=20)

## ElnD: E(ln(det(FIM))) evaluate.
## result is inaccurate (run several times to see)
## increase ED_samp_size for a more accurate calculation
tic();evaluate_design(poped.db,d_switch=FALSE,ED_samp_size=20); toc()

## optimization with line search search 
output_ls <- poped_optim(poped.db, opt_xt=T, parallel=T, method = "LS", d_switch=F, ED_samp_size=20)

## ED: E(det(FIM)) using Laplace approximation 
## deterministic calculation, relatively fast
## can be more stable for optimization
tic(); evaluate_design(poped.db,d_switch=FALSE,use_laplace=TRUE); toc()

## optimization with Laplace
output_ls <- poped_optim(poped.db, opt_xt=T, parallel=T, method = "LS", d_switch=F, use_laplace=T)
