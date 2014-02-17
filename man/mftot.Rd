\name{mftot}
\alias{mftot}
\title{Evaluate the Fisher Information Matrix (FIM)}
\usage{
mftot(model_switch, groupsize, ni, xt, x, a, bpop, d, sigma, docc,
  globalStructure)
}
\arguments{
  \item{globalStructure}{A PopED database.}

  \item{bpop}{The fixed effects parameter values.  Supplied
  as a vector.}

  \item{d}{A between subject variability matrix (OMEGA in
  NONMEM).}

  \item{docc}{A between occasion variability matrix.}

  \item{sigma}{A residual unexplained variability matrix
  (SIGMA in NONMEM).}

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
}
\value{
As a list: \item{ret}{The FIM} \item{globalStructure}{A
PopED database}
}
\description{
Compute the FIM given specific model(s), parameters, design
and methods.
}
\seealso{
For an easier function to use, please see
\code{\link{evaluate.fim}}.

Other FIM: \code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{mf3}};
\code{\link{mf5}}; \code{\link{mf6}}; \code{\link{mf7}};
\code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mf}}
}

