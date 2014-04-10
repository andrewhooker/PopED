\name{blockexp}
\alias{blockexp}
\title{Summarize your experiment for optimization routines}
\usage{
blockexp(fn, poped.db, e_flag = FALSE, opt_xt = poped.db$optsw[2],
  opt_a = poped.db$optsw[4], opt_x = poped.db$optsw[4],
  opt_samps = poped.db$optsw[1], opt_inds = poped.db$optsw[5])
}
\arguments{
  \item{fn}{The file handle to write to.}

  \item{e_flag}{Shuould output be with uncertainty around
  parameters?}

  \item{poped.db}{A PopED database.}

  \item{opt_xt}{Should the sample times be optimized?}

  \item{opt_a}{Should the continuous design variables be
  optimized?}

  \item{opt_x}{Should the discrete design variables be
  optimized?}

  \item{opt_samps}{Are the nuber of sample times per group
  being optimized?}

  \item{opt_inds}{Are the nuber of individuals per group
  being optimized?}
}
\description{
Create some output to the screen and a text file that
summarizes the initial design and the design space you will
use to optimize.
}
\seealso{
Other Helper: \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{blockfinal_2}}; \code{\link{blockfinal}};
\code{\link{blockheader_2}}; \code{\link{blockheader}};
\code{\link{blockopt_2}}; \code{\link{blockopt}}
}

