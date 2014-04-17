\name{ed_mftot}
\alias{ed_mftot}
\title{Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).}
\usage{
ed_mftot(model_switch, groupsize, ni, xtoptn, xoptn, aoptn, bpopdescr, ddescr,
  covd, sigma, docc, poped.db)
}
\arguments{
  \item{xtoptn}{The xtoptn value}

  \item{xoptn}{The xoptn}

  \item{aoptn}{The aoptn value}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{poped.db}{A PopED database.}

  \item{bpopdescr}{Matrix defining the fixed effects, per
  row (row number = parameter_number) we should have:
  \itemize{ \item column 1 the type of the distribution for
  E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform, 3 =
  User Defined Distribution, 4 = lognormal and 5 =
  truncated normal) \item column 2 defines the mean. \item
  column 3 defines the variance of the distribution (or
  length of uniform distribution). }}

  \item{ddescr}{Matrix defining the diagnonals of the IIV
  (same logic as for the \code{bpopdescr}).}

  \item{covd}{Matrix defining the covariances of the IIV
  variances. Set to zero if not defined.}

  \item{sigma}{Matrix defining the variances can
  covariances of the residual variability terms of the
  model. can also just supply the diagnonal parameter
  values (variances) as a \code{c()}.}

  \item{docc}{Matrix defining the IOV, the IOV variances
  and the IOV distribution}
}
\value{
A list containing the E(FIM) and E(OFV(FIM)) and the a
poped.db.
}
\description{
Compute the expectation of the FIM given the model,
parameters, distributions of parameter uncertainty, design
and methods defined in the PopED database.
}
\seealso{
Other E-family: \code{\link{ed_laplace_ofv}};
\code{\link{evaluate.e.ofv.fim}}

Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_laplace_ofv}};
\code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mftot}};
\code{\link{mf}}; \code{\link{ofv_criterion}};
\code{\link{ofv_fim}}
}

