\name{mftot5}
\alias{mftot5}
\title{The full Fisher Information Matrix (FIM) parameterized with A,B,C matrices & using the derivative of variance.}
\usage{
mftot5(model_switch, groupsize, ni, xt, x, a, bpop, d, sigma, docc,
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
Compute the full FIM given specific model(s), parameters,
design and methods. This computation parameterizes the FIM
calculation using A,B,C matrices (as in Retout \emph{et
al.}) but uses the derivative of variances. Should give the
same answer as \code{\link{mftot0}} but computation times
may be different.
}
\references{
S. Retout and F. Mentre, "Further developments of the
Fisher Information Matrix in nonlinear mixed effects models
with evaluation in population pharmacokinetics", J. of
Biopharm. Stats., 13(2), 2003.
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
\code{\link{mftot6}}; \code{\link{mftot7}};
\code{\link{mftot}}; \code{\link{mf}}
}

