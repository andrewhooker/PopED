\name{get_rse}
\alias{get_rse}
\title{Compute the expected parameter relative standard errors}
\usage{
get_rse(fmf, poped.db, bpop = poped.db$gbpop[, 2, drop = F],
  d = poped.db$gd[, 2, drop = F], docc = poped.db$docc,
  sigma = poped.db$sigma, use_percent = T)
}
\arguments{
  \item{use_percent}{Should RSE be reported as percent or
  not?}

  \item{poped.db}{A PopED database.}

  \item{fmf}{The initial value of the FIM. If set to zero
  then it is computed.}

  \item{bpop}{Matrix defining the fixed effects, per row
  (row number = parameter_number) we should have: \itemize{
  \item column 1 the type of the distribution for E-family
  designs (0 = Fixed, 1 = Normal, 2 = Uniform, 3 = User
  Defined Distribution, 4 = lognormal and 5 = truncated
  normal) \item column 2 defines the mean. \item column 3
  defines the variance of the distribution (or length of
  uniform distribution). } Can also just supply the
  parameter values as a vector \code{c()}}

  \item{d}{Matrix defining the diagnonals of the IIV (same
  logic as for the fixed efects). can also just supply the
  parameter values as a \code{c()}.}

  \item{docc}{Matrix defining the IOV, the IOV variances
  and the IOV distribution}

  \item{sigma}{Matrix defining the variances can
  covariances of the residual variability terms of the
  model.}
}
\value{
A named list of RSE values.
}
\description{
This function computes the expected relative standard
errors of a model given a design and a previously computed
FIM.
}
\seealso{
Other evaluate_design: \code{\link{evaluate.fim}};
\code{\link{model_prediction}};
\code{\link{plot_efficiency_of_windows}};
\code{\link{plot_model_prediction}}
}

