\name{mftot4}
\alias{mftot4}
\title{The reduced Fisher Information Matrix (FIM) using the SD of RUV as a parameter.}
\usage{
mftot4(model_switch, groupsize, ni, xt, x, a, bpop, d, sigma, docc,
  globalStructure)
}
\arguments{
  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{xt}{A matrix of sample times.  Each row is a vector
  of sample times for a group.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{bpop}{The fixed effects parameter values.  Supplied
  as a vector.}

  \item{d}{A between subject variability matrix (OMEGA in
  NONMEM).}

  \item{sigma}{A residual unexplained variability matrix
  (SIGMA in NONMEM).}

  \item{docc}{A between occasion variability matrix.}

  \item{globalStructure}{A PopED database.}
}
\value{
As a list: \item{ret}{The FIM} \item{globalStructure}{A
PopED database}
}
\description{
Compute the reduced FIM using the standard deviation of the
residual unexplained variability (RUV) terms as a
parameter, given specific model(s), parameters, design and
methods. This computation assumes that there is no
correlation in the FIM between the fixed and random
effects, and set these elements in the FIM to zero. In
addition all derivatives in the computation are made with
respect to the standard deviation of the RUV terms
(sqrt(SIGMA) in NONMEM). This matches what is done in PFIM,
and assumes that the standard deviation of the residual
unexplained variation is the estimated parameter (NOTE:
NONMEM estimates the variance of the resudual unexplained
variation by default).
}
\seealso{
For an easier function to use, please see
\code{\link{evaluate.fim}}.

Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_laplace_ofv}}; \code{\link{ed_mftot}};
\code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot5}};
\code{\link{mftot6}}; \code{\link{mftot7}};
\code{\link{mftot}}; \code{\link{mf}};
\code{\link{ofv_criterion}}; \code{\link{ofv_fim}}
}

