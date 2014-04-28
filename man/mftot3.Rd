\name{mftot3}
\alias{mftot3}
\title{The Fisher Information Matrix (FIM) some other method}
\usage{
mftot3(model_switch, groupsize, ni, xt, x, a, bpop, d, sigma, docc, poped.db)
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

  \item{poped.db}{A PopED database.}
}
\description{
Compute the FIM using some other method given specific
model(s), parameters, design and methods. This is a
placeholder.
}
\seealso{
For an easier function to use, please see
\code{\link{evaluate.fim}}.

Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_mftot}}; \code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot4}}; \code{\link{mftot5}};
\code{\link{mftot6}}; \code{\link{mftot7}};
\code{\link{mftot}}; \code{\link{mf}};
\code{\link{ofv_criterion}}; \code{\link{ofv_fim}}
}

