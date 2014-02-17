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
  ls_step_size = globalStructure$ls_step_size,
  BFGSConvergenceCriteriaMinStep = popedInput$BFGSConvergenceCriteriaMinStep,
  BFGSProjectedGradientTol = popedInput$BFGSProjectedGradientTol,
  BFGSTolerancef = popedInput$BFGSTolerancef,
  BFGSToleranceg = popedInput$BFGSToleranceg,
  BFGSTolerancex = popedInput$BFGSTolerancex,
  iter_tot = poped.db$iNumSearchIterationsIfNotLineSearch)
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
Other Optimization: \code{\link{RS_opt}};
\code{\link{poped_optimize}}
}

