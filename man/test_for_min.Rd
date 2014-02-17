\name{test_for_min}
\alias{test_for_min}
\title{Test if any matrix element is above a minimum value.}
\usage{
test_for_min(mat, min_mat)
}
\arguments{
  \item{mat}{A matrix.}

  \item{min_mat}{A matrix the same size as mat with the
  minimum allowed value of that element.}
}
\value{
A matrix
}
\description{
Function tests is any matrix element is above a minimum
value.  For those elements the function sets those values
to the minimum value.
}
\seealso{
Other matrix_manipulation: \code{\link{diag_matlab}};
\code{\link{test_for_max}}; \code{\link{test_for_zeros}}
}

