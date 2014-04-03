\name{blockheader}
\alias{blockheader}
\title{Header function for optimization routines}
\usage{
blockheader(bDiscreteOpt, iter, poped.db)
}
\arguments{
  \item{bDiscreteOpt}{is discrete optimization being
  performed?}

  \item{iter}{The number of iterations entered into the
  \code{blockheader_2} function.}

  \item{poped.db}{A PopED database.}
}
\description{
Create some output to the screen and a text file that
summarizes the problem you are tying to solve.
}
\seealso{
Other Helper: \code{\link{blockexp}};
\code{\link{blockfinal_2}}; \code{\link{blockfinal}};
\code{\link{blockheader_2}}; \code{\link{blockopt_2}};
\code{\link{blockopt}}
}

