\name{mf7}
\alias{mf7}
\title{The full Fisher Information Matrix (FIM) for one individual Calculating one model switch at a time, good for large matrices.}
\usage{
mf7(model_switch, xt_ind, x, a, bpop, d, sigma, docc, globalStructure)
}
\arguments{
  \item{model_switch}{A vector that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{xt_ind}{A vector of sample times.}

  \item{x}{A vector for the discrete design variables.}

  \item{a}{A vector of covariates.}

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
As a list: \item{ret}{The FIM for one individual}
\item{globalStructure}{A PopED database}
}
\description{
Compute the full FIM for one individual given specific
model(s), parameters, design and methods. This computation
calculates the FIM for each model switch separately.
Correlations between the models parameters are assumed to
be zero.
}
\seealso{
Used by \code{\link{mftot6}}.

Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_laplace_ofv}}; \code{\link{ed_mftot}};
\code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mftot}};
\code{\link{mf}}; \code{\link{ofv_criterion}};
\code{\link{ofv_fim}}
}

