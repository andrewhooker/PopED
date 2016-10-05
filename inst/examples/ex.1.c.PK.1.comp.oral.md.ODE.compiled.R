library(PopED)

library(deSolve)

## define the ODE
PK.1.comp.oral.ode <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {    
    dA1 <- -KA*A1
    dA2 <- KA*A1 - (CL/V)*A2
    return(list(c(dA1, dA2)))
  })
}

## define the initial conditions and the dosing
ff.ode <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt) #xt[,,drop=T] 
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE*Favail), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, PK.1.comp.oral.ode, parameters, events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
                Favail=bpop[4],
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}


poped.db <- create.poped.database(ff_fun=ff.ode,
                                  fError_fun=feps.add.prop,
                                  fg_fun=sfg,
                                  groupsize=20,
                                  m=2,      #number of groups
                                  sigma=c(0.04,5e-6),
                                  bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                  d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                  notfixed_bpop=c(1,1,1,0),
                                  notfixed_sigma=c(0,0),
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  discrete_xt = list(0:248),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                  maxa=c(DOSE=200,TAU=24),
                                  mina=c(DOSE=0,TAU=24))


plot_model_prediction(poped.db)

# calculations are noticeably slower than with the analytic solution
tic(); evaluate_design(poped.db); toc()



##################################
# Compiled ODE
##################################
# to make reasonable for optimization we can use compiled C code

# compile and load the one_comp_oral_CL.c code.
# to set this up see the 
# "R Package deSolve, Writing Code in Compiled Languages" 
# vingette in the deSolve documentation

# make sure you are in the same directory as one_comp_oral_CL.c
system("R CMD SHLIB one_comp_oral_CL.c")
dyn.load(paste("one_comp_oral_CL", .Platform$dynlib.ext, sep = ""))


#' define the initial conditions and the dosing
ff.ode.compiled <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt)
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE*Favail), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, func = "derivs", parms = c(CL,V,KA),
               #jacfunc = "jac", # not really needed, speed up is minimal if this is defined or not.
               dllname = "one_comp_oral_CL",
               initfunc = "initmod", nout = 1, outnames = "Sum",
               events = list(data = eventdat)) #atol=1e-13,rtol=1e-13))  
    y = out[, "A2"]/(V)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

## -- Update poped.db with compiled function
poped.db.compiled <- create.poped.database(poped.db,ff_fun=ff.ode.compiled)

plot_model_prediction(poped.db.compiled)

# calculations are faster than with the pure R solution
tic(); evaluate_design(poped.db.compiled); toc()

## optimize with line search
output <- poped_optim(poped.db.compiled, opt_xt =TRUE, parallel=TRUE, method = c("LS"))


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


ff.ode.rcpp <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt) #xt[,,drop=T] 
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE*Favail), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, one_comp_oral_ode, c(CL,V,KA), events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

## -- Update poped.db with compiled function
poped.db.compiled.rcpp <- create.poped.database(poped.db,ff_fun=ff.ode.rcpp)

plot_model_prediction(poped.db.compiled.rcpp)

# calculations are much faster than with the pure R solution but slower than .dll compiled solution in desolve
tic(); evaluate_design(poped.db.compiled.rcpp); toc()


  
  