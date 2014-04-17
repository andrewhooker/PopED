\name{blockfinal_2}
\alias{blockfinal_2}
\title{Result function for optimization routines}
\usage{
blockfinal_2(fn, fmf, dmf, groupsize, ni, xt, x, a, model_switch, bpop, d, docc,
  sigma, m, poped.db, opt_xt = poped.db$optsw[2], opt_a = poped.db$optsw[4],
  opt_x = poped.db$optsw[4], fmf_init = NULL, dmf_init = NULL,
  param_cvs_init = NULL)
}
\arguments{
  \item{fmf_init}{Initial FIM.}

  \item{dmf_init}{Initial OFV.}

  \item{param_cvs_init}{The inital design parameter RSE
  values.}

  \item{fmf}{The initial value of the FIM. If set to zero
  then it is computed.}

  \item{dmf}{The inital OFV. If set to zero then it is
  computed.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{xt}{A matrix of sample times.  Each row is a vector
  of sample times for a group.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{poped.db}{A PopED database.}

  \item{opt_xt}{Should the sample times be optimized?}

  \item{opt_a}{Should the continuous design variables be
  optimized?}

  \item{opt_x}{Should the discrete design variables be
  optimized?}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

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
  model. can also just supply the diagnonal parameter
  values (variances) as a \code{c()}.}

  \item{m}{Number of groups/individuals}

  \item{fn}{The file handle to write to.}
}
\description{
Create some output to the screen and a text file that
summarizes the problem you solved.
}
\seealso{
Other Helper: \code{\link{blockexp}};
\code{\link{blockfinal}}; \code{\link{blockheader_2}};
\code{\link{blockheader}}; \code{\link{blockopt_2}};
\code{\link{blockopt}}
}

