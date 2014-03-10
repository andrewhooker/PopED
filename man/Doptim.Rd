\name{Doptim}
\alias{Doptim}
\title{D CONTINUOUS VARIABLE OPTIMIZATION FUNCTION}
\usage{
Doptim(poped.db, ni, xt, model_switch, x, a, bpopdescr, ddescr, maxxt, minxt,
  maxa, mina, fmf, dmf, trflag = TRUE,
  bUseRandomSearch = poped.db$bUseRandomSearch,
  bUseStochasticGradient = poped.db$bUseStochasticGradient,
  bUseBFGSMinimizer = poped.db$bUseBFGSMinimizer,
  bUseLineSearch = poped.db$bUseLineSearch, sgit = poped.db$sgit,
  ls_step_size = poped.db$ls_step_size,
  BFGSConvergenceCriteriaMinStep = poped.db$BFGSConvergenceCriteriaMinStep,
  BFGSProjectedGradientTol = poped.db$BFGSProjectedGradientTol,
  BFGSTolerancef = poped.db$BFGSTolerancef,
  BFGSToleranceg = poped.db$BFGSToleranceg,
  BFGSTolerancex = poped.db$BFGSTolerancex,
  iter_tot = poped.db$iNumSearchIterationsIfNotLineSearch)
}
\arguments{
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
}
\description{
Optimize the objective function. The function works for
continuous optimization variables. This function takes
information from the PopED database supplied as an
argument. The PopED database supplies information about the
the model, parameters, design and methods to use. Some of
the arguments coming from the PopED database can be
overwritten; if they are supplied then they are used
instead of the arguments from the PopED database.
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
Other Optimization: \code{\link{RS_opt}};
\code{\link{poped_optimize}}
}

