% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize_n.R
\name{optimize_n_eff}
\alias{optimize_n_eff}
\title{Translate efficiency to number of subjects}
\usage{
optimize_n_eff(poped.db, ofv_ref, norm_group_fim = NULL, ...)
}
\arguments{
\item{poped.db}{A PopED database.}

\item{ofv_ref}{A reference OFV value to compare to.}

\item{norm_group_fim}{The FIM per individual in each design group. If \code{NULL}, then
these are computed.}

\item{...}{Arguments passed to \code{\link{evaluate.fim}} and \code{efficiency}.}
}
\value{
The number of individuals needed.
}
\description{
optimize HOW MANY n there should be to achieve efficiency=1 compared to a reference OFV
}
\examples{
# 2 design groups with either early or late samples
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
                                  sigma=c(0.01,0.25),
                                  xt=list(c(1,2,3),c(4,5,20,120)),
                                  groupsize=50,
                                  minxt=0.01,
                                  maxxt=120,
                                  a=70,
                                  mina=0.01,
                                  maxa=100)


plot_model_prediction(poped.db)

evaluate_design(poped.db)


# what are the optimal proportions of 
# individuals in the two groups in the study?
(n_opt <- optimize_groupsize(poped.db))

# How many individuals in the original design are needed to achieve an
# efficiency of 1 compared to the optimized design with n=100?
optimize_n_eff(poped.db,
               ofv_ref=n_opt$opt_ofv_with_n)



}
