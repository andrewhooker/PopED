\name{evaluate.e.ofv.fim}
\alias{evaluate.e.ofv.fim}
\title{Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).}
\usage{
evaluate.e.ofv.fim(poped.db, fim.calc.type = NULL, bpop = poped.db$gbpop,
  d = poped.db$gd, covd = poped.db$covd, docc = poped.db$docc,
  sigma = poped.db$sigma, model_switch = NULL, ni = NULL, xt = NULL,
  x = NULL, a = NULL, groupsize = poped.db$design$groupsize,
  deriv.type = NULL, bLHS = poped.db$bLHS,
  ofv_calc_type = poped.db$ofv_calc_type,
  ED_samp_size = poped.db$ED_samp_size, use.laplace = FALSE, ...)
}
\arguments{
  \item{use.laplace}{Should the Laplace method be used in
  calculating the expectation?  Currently experimental for
  the R-version of PopED.}

  \item{poped.db}{A PopED database.}

  \item{fim.calc.type}{The method used for calculating the
  FIM. Potential values: \itemize{ \item 0 = Full FIM.  No
  assumption that fixed and random effects are
  uncorrelated.  See \code{\link{mftot0}}. \item 1 =
  Reduced FIM. Assume that there is no correlation in the
  FIM between the fixed and random effects, and set these
  elements in the FIM to zero. See \code{\link{mftot1}}.
  \item 2 = weighted models (placeholder). \item 3 = Not
  currently used. \item 4 = Reduced FIM and computing all
  derivatives with respect to the standard deviation of the
  residual unexplained variation (sqrt(SIGMA) in NONMEM).
  This matches what is done in PFIM, and assumes that the
  standard deviation of the residual unexplained variation
  is the estimated parameter (NOTE: NONMEM estimates the
  variance of the resudual unexplained variation by
  default). See \code{\link{mftot4}}. \item 5 = Full FIM
  parameterized with A,B,C matrices & derivative of
  variance. See \code{\link{mftot5}}. \item 6 = Calculate
  one model switch at a time, good for large matrices. See
  \code{\link{mftot6}}. \item 7 = Reduced FIM parameterized
  with A,B,C matrices & derivative of variance See
  \code{\link{mftot7}}. }}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{xt}{A matrix of sample times.  Each row is a vector
  of sample times for a group.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

  \item{deriv.type}{A number indicating the type of
  derivative to use: \itemize{ \item 0=Complex difference
  \item 1=Central difference \item 20=Analytic derivative
  (placeholder) \item 30=Automatic differentiation
  (placeholder) }}

  \item{...}{Other arguments passed to the function.}

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

  \item{covd}{Matrix defining the covariances of the IIV
  variances. Set to zero if not defined.}

  \item{docc}{Matrix defining the IOV, the IOV variances
  and the IOV distribution}

  \item{sigma}{Matrix defining the variances can
  covariances of the residual variability terms of the
  model.}

  \item{bLHS}{How to sample from distributions in E-family
  calculations. 0=Random Sampling, 1=LatinHyperCube --}

  \item{ofv_calc_type}{OFV calculation type for FIM
  (1=Determinant of FIM,4=log determinant of
  FIM,6=determinant of interesting part of FIM (Ds))}

  \item{ED_samp_size}{Sample size for E-family sampling}
}
\value{
A list containing the E(FIM) and E(OFV(FIM)) and the a
poped.db updated according to the function arguments.
}
\description{
Compute the expectation of the FIM and OFV(FIM) given the
model, parameters, distributions of parameter uncertainty,
design and methods defined in the PopED database. Some of
the arguments coming from the PopED database can be
overwritten; by default these arguments are \code{NULL} in
the function, if they are supplied then they are used
instead of the arguments from the PopED database.
}
\examples{
\dontrun{
evaluate.ed.fim(poped.db)
}
}
\seealso{
Other E-family: \code{\link{ed_laplace_ofv}};
\code{\link{ed_mftot}}

Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_laplace_ofv}}; \code{\link{ed_mftot}};
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

