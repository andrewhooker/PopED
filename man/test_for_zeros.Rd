\name{test_for_zeros}
\alias{test_for_zeros}
\title{Test if any matrix element is zero.}
\usage{
test_for_zeros(mat, ourzero)
}
\arguments{
  \item{mat}{A matrix.}

  \item{ourzero}{A matrix the same size as mat with the
  value that zero should be reassigned to.}
}
\value{
A matrix
}
\description{
Function tests is any matrix element is zero.  For those
elements the function sets those values to the minimum
value allowed (not zero). This is used to avoid numerical
problems in the FIM calculation.
}
\seealso{
Other matrix_manipulation: \code{\link{diag_matlab}};
\code{\link{test_for_max}}; \code{\link{test_for_min}}
}

