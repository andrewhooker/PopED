\name{LinMatrixL}
\alias{LinMatrixL}
\title{The linearized matrix L}
\usage{
LinMatrixL(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, globalStructure)
}
\arguments{
  \item{bpop}{The fixed effects parameter values.  Supplied
  as a vector.}

  \item{b_ind}{The point at which to evaluate the
  derivative}

  \item{bocc_ind}{The point at which to evaluate the
  derivative}

  \item{globalStructure}{A PopED database.}

  \item{model_switch}{A vector that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{xt_ind}{A vector of sample times.}

  \item{x}{A vector for the discrete design variables.}

  \item{a}{A vector of covariates.}
}
\value{
As a list: \item{y}{A matrix of size (samples per
individual x number of random effects)}
\item{globalStructure}{A PopED database}
}
\description{
Function computes the derivative of the model with respect
to the between subject variability terms in the model (b's
and bocc's) evaluated at a defined point (b_ind and
bocc_ind).
}

