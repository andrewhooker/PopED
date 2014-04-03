\name{Dtrace}
\alias{Dtrace}
\title{Trace optimization routines}
\usage{
Dtrace(fn, it, ni, xtopt, xopt, aopt, gxt, ga, dmf, diff, ixt, ia, itvector,
  dmfvector, poped.db, opt_xt = poped.db$optsw[2],
  opt_a = poped.db$optsw[4], opt_x = poped.db$optsw[3],
  opt_samps = poped.db$optsw[1], opt_inds = poped.db$optsw[5],
  rsit = poped.db$rsit)
}
\arguments{
  \item{fn}{A file to output information to. Can also be
  the screen if \code{''}.}

  \item{it}{the interation number.}

  \item{xtopt}{The matrix defining optimal sampling
  schedule.}

  \item{xopt}{The cell structure defining the optimal
  discrete design variables.}

  \item{aopt}{The matrix defining the optimal continuous
  design variables.}

  \item{gxt}{The matrix defining the current sampling
  schedule.}

  \item{ga}{The matrix defining the current continuous
  design variables.}

  \item{dmf}{The current OFV.}

  \item{diff}{The difference from the previous iteration.}

  \item{ixt}{If xt Gradient Inversion Occured or not.}

  \item{ia}{If a Gradient Inversion Occured or not.}

  \item{itvector}{The iteration vector}

  \item{dmfvector}{The dmf vector.}

  \item{opt_samps}{Are the nuber of sample times per group
  being optimized?}

  \item{opt_inds}{Are the nuber of individuals per group
  being optimized?}

  \item{opt_xt}{Should the sample times be optimized?}

  \item{opt_a}{Should the continuous design variables be
  optimized?}

  \item{opt_x}{Should the discrete design variables be
  optimized?}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{poped.db}{A PopED database.}

  \item{rsit}{Number of Random search iterations}
}
\description{
A helper function for writing output to the screen and
files when optimizing.
}

