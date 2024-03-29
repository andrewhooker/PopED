% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_unfixed_params.R
\name{get_unfixed_params}
\alias{get_unfixed_params}
\title{Return all the unfixed parameters}
\usage{
get_unfixed_params(poped.db, params = NULL)
}
\arguments{
\item{poped.db}{a PopED database.}

\item{params}{If params is supplied then the parameters are taken from this vector.}
}
\value{
A list with the  parameters.  All unfixed parameters are also
  returned in the "\code{all} output with the specified order 
  (bpop,d,covd,docc,covdocc,sigma,covsigma). \code{var_derivative}  is a
  vector of 1's or 0's, 1 means derivative of parameter is taken with respect
  to the variance otherwise with respect to standard deviation.
}
\description{
all = vector of all unfixed params var derivative is a vector of 1 and 0, 1
means derivative of parameter is taken w.r.t. variance otherwise w.r.t. sd If
params is supplied then the parameter is taken from this vector instead of
poped.db
}
\keyword{internal}
