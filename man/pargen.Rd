\name{pargen}
\alias{pargen}
\title{Parameter simulation}
\usage{
pargen(par, user_dist_pointer, sample_size, bLHS, sample_number,
  globalStructure)
}
\arguments{
  \item{par}{A matrix describing the parameters. Each row
  is a parameter and the matrix has three columns:
  \enumerate{ \item First column - Type of distribution
  (0-fixed, 1-normal, 2-uniform, 3-user specified,
  4-lognormal, 5-Truncated normal).  \item Second column -
  Mean of distribution.  \item Third column - Variance or
  range of distribution.  }}

  \item{user_dist_pointer}{A text string of the name of a
  function that generates random samples from a user
  defined distribution.}

  \item{sample_size}{The number of random samples per
  parameter to generate}

  \item{bLHS}{Logical, indicating if Latin Hypercube
  Sampling should be used.}

  \item{sample_number}{The sample number to extract.}

  \item{globalStructure}{A PopED database.}
}
\value{
A matrix of random samples of size (sample_size x
number_of_parameters)
}
\description{
Function generates random samples for a list of parameters
}

