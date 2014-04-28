\name{feval}
\alias{feval}
\title{MATLAB feval function}
\usage{
feval(file.name, ...)
}
\arguments{
  \item{file.name}{A string that is either a function or a
  file name (with an .R extension).}

  \item{...}{Arguments for the function.  Multiple
  arguments separated by a comma.}
}
\value{
Output from the defined function.
}
\description{
This is just a wrapper for the \code{\link{do.call}}
function to behave like the feval function in MATLAB.
}
\examples{

feval("sin",pi/2)

}
\seealso{
Other MATLAB: \code{\link{cell}};
\code{\link{diag_matlab}}; \code{\link{fileparts}};
\code{\link{isempty}}; \code{\link{isfield}};
\code{\link{ones}}; \code{\link{randn}};
\code{\link{randperm}}; \code{\link{rand}};
\code{\link{size}}; \code{\link{tic}}; \code{\link{toc}};
\code{\link{zeros}}
}

