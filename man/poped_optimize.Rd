\name{poped_optimize}
\alias{poped_optimize}
\title{Optimization main module for PopED}
\usage{
poped_optimize(poped.db, ni = NULL, xt = NULL, model_switch = NULL,
  x = NULL, a = NULL, bpop = NULL, d = NULL, maxxt = NULL,
  minxt = NULL, maxa = NULL, mina = NULL, fmf = 0, dmf = 0,
  trflag = TRUE, opt_xt = poped.db$optsw[2], opt_a = poped.db$optsw[4],
  opt_x = poped.db$optsw[3], opt_samps = poped.db$optsw[1],
  opt_inds = poped.db$optsw[5], cfaxt = poped.db$cfaxt,
  cfaa = poped.db$cfaa, rsit = poped.db$rsit,
  rsit_output = poped.db$rsit_output,
  fim.calc.type = poped.db$iFIMCalculationType,
  ofv_calc_type = poped.db$ofv_calc_type,
  approx_type = poped.db$iApproximationMethod,
  bUseExchangeAlgorithm = poped.db$bUseExchangeAlgorithm, iter = 1, ...)
}
\arguments{
  \item{...}{arguments passed to other functions. See
  \link{\code{Doptim}}.}

  \item{poped.db}{A PopED database.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{xt}{A matrix of sample times.  Each row is a vector
  of sample times for a group.}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{cfaxt}{First step factor for sample times}
}
\description{
Optimize the objective function. The function works for
both discrete and continuous optimization variables. This
function takes information from the PopED database supplied
as an argument. The PopED database supplies information
about the the model, parameters, design and methods to use.
Some of the arguments coming from the PopED database can be
overwritten; if they are supplied then they are used
instead of the arguments from the PopED database.
}
\examples{
\dontrun{
}
}
\references{
\enumerate{ \item M. Foracchia, A.C. Hooker, P. Vicini and
A. Ruggeri, "PopED, a software fir optimal experimental
design in population kinetics", Computer Methods and
Programs in Biomedicine, 74, 2004. \item J. Nyberg, S.
Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C.
Hooker, "PopED: An extended, parallelized, nonlinear mixed
effects models optimal design tool", Computer Methods and
Programs in Biomedicine, 108, 2012. }
}
\seealso{
Other Optimization: \code{\link{Doptim}};
\code{\link{RS_opt}}
}

