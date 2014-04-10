\name{a_line_search}
\alias{a_line_search}
\title{Optimize using line search}
\usage{
a_line_search(out_file, bED, diff, temp_mf, temp_detmf, globalStructure)
}
\arguments{
  \item{out_file}{The output file to write to.}

  \item{diff}{The stopping criteria for the algorithm.  If,
  after one iteration the difference in OFV is less than
  \code{diff} then the algorithm stops.}

  \item{temp_mf}{The initial value of the FIM.}

  \item{temp_detmf}{The initial value of the det(FIM).}

  \item{globalStructure}{A PopED database.}

  \item{bED}{If the algorithm should use E-family methods.
  Logical.}
}
\value{
A list containing: \item{fmf}{The FIM.} \item{dmf}{The
determinant of the FIM.} \item{best_changed}{If the
algorithm has found a better design than the starting
design.} \item{xt}{A matrix of sample times.  Each row is a
vector of sample times for a group.} \item{x}{A matrix for
the discrete design variables.  Each row is a group.}
\item{a}{A matrix of covariates.  Each row is a group.}
\item{globalStructure}{A PopED database.}
}
\description{
The function performs a grid search sequatially along
design variables.
}
\seealso{
Other Optimize: \code{\link{Doptim}}; \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{bfgsb_min}}; \code{\link{calc_autofocus}};
\code{\link{calc_ofv_and_grad}}; \code{\link{mfea}};
\code{\link{poped_optimize}}
}

