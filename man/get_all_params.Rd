\name{get_all_params}
\alias{get_all_params}
\title{Extract all model parameters from the PopED database.}
\usage{
get_all_params(poped.db)
}
\arguments{
  \item{poped.db}{A PopED database.}
}
\value{
A list containing: \item{bpop}{A vector of fixed effect
parameter values.} \item{d}{A vector of between subject
variability parameters} \item{covd}{A vector of the
covariances of the between subject variability parameters.
Row major format of the lower triangular portion of the D
(OMEGA) matrix} \item{docc}{A vector of the between
occasion variability (BOV) terms in the model}
\item{covdocc}{A vector of the covariances between the BOV
terms.  Row major of the lower triangular portion of the
BOV matrix. } \item{sigma}{A vector of the resudual
unexplained variances (RUV)} \item{covsigma}{A vector of
the covariances between the RUV terms} \item{all}{A vector
with all of the above, in the order of this list.}
}
\description{
Extract all model parameters from the PopED database.
}

