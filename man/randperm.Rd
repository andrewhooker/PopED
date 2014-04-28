\name{randperm}
\alias{randperm}
\title{Function written to match MATLAB's randperm function}
\usage{
randperm(num)
}
\arguments{
  \item{num}{Either a vector of one or more elements from
  which to choose, or a positive integer.}
}
\value{
See \code{\link{sample}}
}
\description{
A wrapper for the \code{\link{sample}} function.
}
\examples{

randperm(c(2,3,4,5,6))

randperm(10)
}
\seealso{
Other MATLAB: \code{\link{cell}};
\code{\link{diag_matlab}}; \code{\link{feval}};
\code{\link{fileparts}}; \code{\link{isempty}};
\code{\link{isfield}}; \code{\link{ones}};
\code{\link{randn}}; \code{\link{rand}};
\code{\link{size}}; \code{\link{tic}}; \code{\link{toc}};
\code{\link{zeros}}
}

