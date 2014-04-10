\name{blockopt_2}
\alias{blockopt_2}
\title{Summarize your optimization settings for optimization routines}
\usage{
blockopt_2(fn, poped.db, opt_method = "")
}
\arguments{
  \item{opt_method}{If "RS" (random search), "SG"
  (stochastic gradient) or "DO" (discrete optimization)
  then specifc output is produced.}

  \item{poped.db}{A PopED database.}

  \item{fn}{The file handle to write to.}
}
\description{
Create some output to the screen and a text file that
summarizes the optimization settings you will use to
optimize.
}
\seealso{
Other Helper: \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{blockexp}}; \code{\link{blockfinal_2}};
\code{\link{blockfinal}}; \code{\link{blockheader_2}};
\code{\link{blockheader}}; \code{\link{blockopt}}
}

