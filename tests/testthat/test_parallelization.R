context("Parallelization")

test_that("plot_efficiency_of_windows works in parallel", {
  
  if(skip_parallel) skip("Parallel calculations skipped")
  
  sfg <- function(x,a,bpop,b,bocc){
    parameters=c(CL=bpop[1]*exp(b[1]),
                 V=bpop[2]*exp(b[2]),
                 KA=bpop[3]*exp(b[3]),
                 Favail=bpop[4],
                 DOSE=a[1])
    return(parameters) 
  }
  
  ## -- Define model, parameters, initial design
  poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                    fg_fun=sfg,
                                    fError_fun=feps.add.prop,
                                    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6), 
                                    sigma=c(prop=0.01,add=0.25),
                                    groupsize=32,
                                    xt=c( 0.5,1,2,6,24,36,72,120),
                                    minxt=0.01,
                                    maxxt=120,
                                    a=c(DOSE=70),
                                    mina=c(DOSE=0.01),
                                    maxa=c(DOSE=100))
  
  t1 <- system.time(p1 <- plot_efficiency_of_windows(poped.db,xt_windows=0.5,parallel = T,iNumSimulations = 100))
  t2 <- system.time(p1 <- plot_efficiency_of_windows(poped.db,xt_windows=0.5,parallel = F,iNumSimulations = 100))
  
  expect_lt(t1["elapsed"],t2["elapsed"])
  
  
})

