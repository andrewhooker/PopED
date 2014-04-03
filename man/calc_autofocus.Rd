\name{calc_autofocus}
\alias{calc_autofocus}
\title{Compute the autofocus portion of the stochastic gradient routine}
\usage{
calc_autofocus(m, ni_var, dmf, varopt, varopto, maxvar, minvar, gradvar,
  normgvar, avar, model_switch, groupsize, xtopt, xopt, aopt, ni, bpop, d,
  sigma, docc, poped.db)
}
\arguments{
  \item{ni_var}{The ni_var.}

  \item{varopt}{The varopt.}

  \item{varopto}{The varopto.}

  \item{maxvar}{The maxvar.}

  \item{minvar}{The minvar.}

  \item{gradvar}{The gradvar.}

  \item{normgvar}{The normgvar.}

  \item{avar}{The avar.}

  \item{xtopt}{The optimal sampling times matrix.}

  \item{xopt}{The optimal discrete design variables
  matrix.}

  \item{aopt}{The optimal continuous design variables
  matrix.}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{poped.db}{A PopED database.}

  \item{dmf}{The inital OFV. If set to zero then it is
  computed.}

  \item{m}{Number of groups/individuals}

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

  \item{sigma}{Matrix defining the variances can
  covariances of the residual variability terms of the
  model.}

  \item{docc}{Matrix defining the IOV, the IOV variances
  and the IOV distribution}
}
\value{
A list containing: \item{navar}{The autofocus parameter.}
\item{poped.db}{PopED database.}
}
\description{
Compute the autofocus portion of the stochastic gradient
routine
}
\seealso{
Other Optimization: \code{\link{Doptim}};
\code{\link{RS_opt}}; \code{\link{calc_ofv_and_grad}};
\code{\link{mfea}}; \code{\link{poped_optimize}}
}

