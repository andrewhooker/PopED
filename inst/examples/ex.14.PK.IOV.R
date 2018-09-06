#library(PopED)

##-- Model: One comp first order absorption
library(Rcpp)

cppFunction(
  'List one_comp_oral_ode(double Time, NumericVector A, NumericVector Pars) {
   int n = A.size();
   NumericVector dA(n);
            
   double CL_OCC_1 = Pars[0];
   double CL_OCC_2 = Pars[1];
   double V = Pars[2];
   double KA = Pars[3];
   double TAU = Pars[4];
   double N,CL;
            
   N = floor(Time/TAU)+1;
   CL = CL_OCC_1;
   if(N>6) CL = CL_OCC_2;
   
   dA[0] = -KA*A[0];
   dA[1] = KA*A[0] - (CL/V)*A[1];
   return List::create(dA);
   }'
)

ff.ode.rcpp <- function(model_switch, xt, parameters, poped.db){
  with(as.list(parameters),{
    A_ini <- c(A1=0, A2=0)
    times_xt <- drop(xt) #xt[,,drop=T] 
    dose_times = seq(from=0,to=max(times_xt),by=TAU)
    eventdat <- data.frame(var = c("A1"), 
                           time = dose_times,
                           value = c(DOSE), method = c("add"))
    times <- sort(c(times_xt,dose_times))
    out <- ode(A_ini, times, one_comp_oral_ode, c(CL_OCC_1,CL_OCC_2,V,KA,TAU), 
               events = list(data = eventdat))#atol=1e-13,rtol=1e-13)
    y = out[, "A2"]/(V)
    y=y[match(times_xt,out[,"time"])]
    y=cbind(y)
    return(list(y=y,poped.db=poped.db))
  })
}

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c( CL_OCC_1=bpop[1]*exp(b[1]+bocc[1,1]),
                CL_OCC_2=bpop[1]*exp(b[1]+bocc[1,2]),
                V=bpop[2]*exp(b[2]),
                KA=bpop[3]*exp(b[3]),
                DOSE=a[1],
                TAU=a[2])
  return( parameters ) 
}


## -- Define design and design space
poped.db <- create.poped.database(ff_fun=ff.ode.rcpp,
                                  fError_fun=feps.add.prop,
                                  fg_fun=sfg,
                                  bpop=c(CL=3.75,V=72.8,KA=0.25), 
                                  d=c(CL=0.25^2,V=0.09,KA=0.09), 
                                  sigma=c(0.04,5e-6),
                                  notfixed_sigma=c(0,0),
                                  docc = matrix(c(0,0.09,0),nrow = 1),
                                  m=2,
                                  groupsize=20,
                                  xt=c( 1,2,8,240,245),
                                  minxt=c(0,0,0,240,240),
                                  maxxt=c(10,10,10,248,248),
                                  bUseGrouped_xt=1,
                                  a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
                                  maxa=c(DOSE=200,TAU=24),
                                  mina=c(DOSE=0,TAU=24))



##  create plot of model without variability 
plot_model_prediction(poped.db, model_num_points = 300)

##  Visulaize the IOV 
library(ggplot2)
set.seed(12345679)
plot_model_prediction(poped.db, PRED=F,IPRED=F, 
                      separate.groups=T, model_num_points = 500, 
                      groupsize_sim = 1,
                      IPRED.lines = T, alpha.IPRED.lines=0.6,
                      sample.times = F
) + geom_vline(xintercept = 24*6,color="red")

## evaluate initial design
evaluate_design(poped.db)
