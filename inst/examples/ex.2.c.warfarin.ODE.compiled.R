## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error to 
##   avoid sample times at very low concentrations (time 0 or very late samoples).

## Model described with an ODE
library(PopED)
library(deSolve)

sfg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

# define the ODE
one.comp.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {    
    dA1  <- -KA*A1
    dA2  <- KA*A1 - CL/V*A2
    return(list(c(dA1, dA2)))
  })
}

ff.ODE <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE*Favail, A2 = 0)
    times <- drop(xt)##xt[,,drop=T] 
    times <- sort(times) 
    times <- c(0,times) ## add extra time for start of integration
    out   <- ode(A_ini, times, one.comp.ode, parameters)#,atol=1e-13,rtol=1e-13)
    y = out[,"A2"]/(V)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}


## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff.ODE,
                                  fg_fun = sfg,
                                  fError_fun = feps.add.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=c(DOSE=70),
                                  mina=c(DOSE=0),
                                  maxa=c(DOSE=100))

##  create plot of model without variability 
plot_model_prediction(poped.db)

##  create plot of model with variability 
plot_model_prediction(poped.db,IPRED=T,DV=T)

## evaluate initial design
# speed is quite slow with the ODE model if not using compiled code (see below)
tic(); design_ode <- evaluate_design(poped.db); toc()


##################################
# Compiled ODE
##################################

# to make speed reasonable for optimization we can use compiled C code

# compile and load the one_comp_oral_CL.c code.
# to set this up see the 
# "R Package deSolve, Writing Code in Compiled Languages" 
# vingette in the deSolve documentation

# make sure you are in the same directory as one_comp_oral_CL.c
system("R CMD SHLIB one_comp_oral_CL.c")
dyn.load(paste("one_comp_oral_CL", .Platform$dynlib.ext, sep = ""))

ff.ODE.compiled <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE*Favail, A2 = 0)
    times <- drop(xt)##xt[,,drop=T] 
    times <- sort(times) 
    times <- c(0,times) ## add extra time for start of integration
    out <- ode(A_ini, times, func = "derivs", parms = c(CL,V,KA),
               #jacfunc = "jac", # not really needed, speed up is minimal if this is defined or not.
               dllname = "one_comp_oral_CL",
               initfunc = "initmod", nout = 1, outnames = "Sum")    
    y = out[,"A2"]/(V)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}

## -- Update poped.db with compiled function
poped.db.compiled <- create.poped.database(poped.db, ff_fun=ff.ODE.compiled)

##  create plot of model without variability 
plot_model_prediction(poped.db.compiled)

##  create plot of model with variability 
plot_model_prediction(poped.db.compiled,IPRED=T,DV=T)

## evaluate initial design (much faster than pure R solution)
tic(); design_ode_compiled <- evaluate_design(poped.db.compiled); toc()

# no difference in computation
design_ode_compiled[["ofv"]] - design_ode[["ofv"]]

## making optimization times more resonable
output <- poped_optim(poped.db.compiled, opt_xt =TRUE, parallel=TRUE, method = c("LS"))


##################################
# Compiled ODE using Rcpp
################################## 

#You could also use somewhat slower, but more readable, inline code from Rcpp
library(Rcpp)

cppFunction('List one_comp_oral_ode(double Time, NumericVector A, NumericVector Pars) {
            int n = A.size();
            NumericVector dA(n);
            
            double CL = Pars[0];
            double V = Pars[1];
            double KA = Pars[2];
            
            dA[0] = -KA*A[0];
            dA[1] = KA*A[0] - (CL/V)*A[1];
            return List::create(dA);
            }')


ff.ode.rcpp <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp, linear abssorption, single dose
  ##-- Parameterized by CL, KA, F and V.
  with(as.list(parameters),{
    A_ini  <- c(A1 = DOSE*Favail, A2 = 0)
    times <- drop(xt)##xt[,,drop=T] 
    times <- sort(times) 
    times <- c(0,times) ## add extra time for start of integration
    out   <- ode(A_ini, times, one_comp_oral_ode, c(CL,V,KA))#,atol=1e-13,rtol=1e-13)
    y = out[,"A2"]/(V)
    y=y[-1] # remove initial time for start of integration
    y = cbind(y) # must be a column matrix 
    return(list( y= y,poped.db=poped.db)) 
  })
}

## -- Update poped.db with compiled function
poped.db.compiled.rcpp <- create.poped.database(poped.db,ff_fun=ff.ode.rcpp)

## evaluate initial design (much faster than pure R solution, slower than desilve method)
tic(); design_ode_compiled_rcpp <- evaluate_design(poped.db.compiled.rcpp); toc()

##################################
# comapre to the analytic solution
##################################
ff <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list(y=y,poped.db=poped.db))
  })
}

## -- Update poped.db with analytic function
poped.db.analytic <- create.poped.database(poped.db,ff_fun=ff)

## computation times are ~6x faster with the analytic solution
tic(); design_analytic <- evaluate_design(poped.db.analytic); toc()

## differences can be reduced by decreasing atol and rtol in the desolve "ode" function 
(design_analytic$ofv - design_ode_compiled$ofv)/design_ode_compiled$ofv*100
