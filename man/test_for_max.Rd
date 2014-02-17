\name{test_for_max}
\alias{test_for_max}
\title{Test if any matrix element is above a max value.}
\usage{
test_for_max(mat, max_mat)
}
\arguments{
  \item{mat}{A matrix.}

  \item{max_mat}{A matrix the same size as mat with the
  maximum allowed value of that element.}
}
\value{
A matrix
}
\description{
Function tests is any matrix element is above a maximum
value.  For those elements the function sets those values
to the maximum value.
}
\seealso{
Other matrix_manipulation: \code{\link{diag_matlab}};
\code{\link{test_for_min}}; \code{\link{test_for_zeros}}
}

