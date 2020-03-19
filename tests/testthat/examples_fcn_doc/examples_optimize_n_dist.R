poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=function(x,a,bpop,b,bocc){
                                    parameters=c(CL=bpop[1]*exp(b[1]),
                                                 V=bpop[2]*exp(b[2]),
                                                 KA=bpop[3]*exp(b[3]),
                                                 Favail=bpop[4],
                                                 DOSE=a[1])
                                    return(parameters) 
                                  },
                                  fError_fun=feps.add.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(PROP=0.01,ADD=0.25),
                                  xt=list(c(1,2,3),c(4,5,20,120)),
                                  groupsize=20,
                                  minxt=0.01,
                                  maxxt=120,
                                  a=70,
                                  mina=0.01,
                                  maxa=100)




plot_model_prediction(poped.db)

(n_opt <- optimize_n_dist(poped.db))
