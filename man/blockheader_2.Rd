\name{blockheader_2}
\alias{blockheader_2}
\title{Header function for optimization routines}
\usage{
blockheader_2(name, iter, poped.db, e_flag = FALSE,
  opt_xt = poped.db$optsw[2], opt_a = poped.db$optsw[4],
  opt_x = poped.db$optsw[4], opt_samps = poped.db$optsw[1],
  opt_inds = poped.db$optsw[5], fmf = 0, dmf = 0, bpop = NULL,
  d = NULL, docc = NULL, sigma = NULL, ...)
}
\arguments{
  \item{name}{The name used for the output file. Combined
  with \code{iter} and some other text.}

  \item{iter}{The last number in the name printed to the
  output file, combined with \code{name}.}

  \item{poped.db}{A PopED database.}

  \item{opt_xt}{Should the sample times be optimized?}

  \item{opt_a}{Should the continuous design variables be
  optimized?}

  \item{opt_x}{Should the discrete design variables be
  optimized?}

  \item{fmf}{The initial value of the FIM. If set to zero
  then it is computed.}

  \item{dmf}{The inital OFV. If set to zero then it is
  computed.}

  \item{...}{arguments passed to
  \code{\link{evaluate.fim}}.}

  \item{e_flag}{Shuould output be with uncertainty around
  parameters?}

  \item{opt_samps}{Are the nuber of sample times per group
  being optimized?}

  \item{opt_inds}{Are the nuber of individuals per group
  being optimized?}

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
\description{
Create some output to the screen and a text file that
summarizes the problem you are tying to solve.
}
\seealso{
Other Helper: \code{\link{blockexp}};
\code{\link{blockfinal_2}}; \code{\link{blockfinal}};
\code{\link{blockheader}}; \code{\link{blockopt_2}};
\code{\link{blockopt}}
}

