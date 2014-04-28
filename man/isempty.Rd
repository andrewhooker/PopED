\name{isempty}
\alias{isempty}
\title{Function written to match MATLAB's isempty function}
\usage{
isempty(...)
}
\arguments{
  \item{...}{arguments to pass to the function. Typically a
  matrix.}
}
\value{
Logical. True if the passed object has any dimension that
is zero.
}
\description{
Function written to match MATLAB's isempty function
}
\examples{

isempty(zeros(2,3))

isempty(zeros(2,0))

isempty(c(1,2,3))

}
\seealso{
Other MATLAB: \code{\link{cell}};
\code{\link{diag_matlab}}; \code{\link{feval}};
\code{\link{fileparts}}; \code{\link{isfield}};
\code{\link{ones}}; \code{\link{randn}};
\code{\link{randperm}}; \code{\link{rand}};
\code{\link{size}}; \code{\link{tic}}; \code{\link{toc}};
\code{\link{zeros}}
}

