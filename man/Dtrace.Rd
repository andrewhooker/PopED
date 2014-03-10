\name{Dtrace}
\alias{Dtrace}
\title{Trace optimization routines}
\usage{
Dtrace(fn, it, ni, xtopt, xopt, aopt, gxt, ga, dmf, diff, ixt, ia, itvector,
  dmfvector, globalStructure, opt_xt = globalStructure$optsw[2],
  opt_a = globalStructure$optsw[4], opt_x = globalStructure$optsw[3],
  opt_samps = globalStructure$optsw[1], opt_inds = globalStructure$optsw[5],
  rsit = globalStructure$rsit)
}
\arguments{
  \item{fn}{A file to output information to. Can also be
  the screen if \code{''}.}

  \item{it}{the interation number.}

  \item{ni}{A vector of the number of samples in each
  group.}
}
\description{
A helper function for writing output to the screen and
files when optimizing.
}

