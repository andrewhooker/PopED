% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/size.R
\name{size}
\alias{size}
\title{Function written to match MATLAB's size function}
\usage{
size(obj, dimension.index = NULL)
}
\arguments{
\item{obj}{An object you want to know the various dimensions of.  Typically a matrix.}

\item{dimension.index}{Which dimension you are interested in.}
}
\value{
The dimensions of the object or specific dimension you are interested in.
}
\description{
Function written to match MATLAB's size function
}
\examples{

size(c(2,3,4,5,6))

size(10)

size(zeros(4,7))

}
\seealso{
Other MATLAB: 
\code{\link{cell}()},
\code{\link{diag_matlab}()},
\code{\link{feval}()},
\code{\link{fileparts}()},
\code{\link{isempty}()},
\code{\link{ones}()},
\code{\link{rand}()},
\code{\link{randn}()},
\code{\link{tic}()},
\code{\link{toc}()},
\code{\link{zeros}()}
}
\concept{MATLAB}
