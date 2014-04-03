\name{Doptim}
\alias{Doptim}
\title{D CONTINUOUS VARIABLE OPTIMIZATION FUNCTION}
\usage{
Doptim(poped.db, ni, xt, model_switch, x, a, bpopdescr, ddescr, maxxt, minxt,
  maxa, mina, fmf = 0, dmf = 0, trflag = TRUE,
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
  \item{bpopdescr}{Matrix defining the fixed effects, per
  row (row number = parameter_number) we should have:
  \itemize{ \item column 1 the type of the distribution for
  E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform, 3 =
  User Defined Distribution, 4 = lognormal and 5 =
  truncated normal) \item column 2 defines the mean. \item
  column 3 defines the variance of the distribution (or
  length of uniform distribution). }}

  \item{ddescr}{Matrix defining the diagnonals of the IIV
  (same logic as for the \code{bpopdescr}).}

  \item{fmf}{The initial value of the FIM. If set to zero
  then it is computed.}

  \item{dmf}{The inital OFV. If set to zero then it is
  computed.}

  \item{trflag}{Should the optimization be output to the
  screen and to a file?}

  \item{ls_step_size}{Number of grid points in the line
  search}

  \item{iter_tot}{Number of iterations of full Random
  search and full Stochastic Gradient if line search is not
  used.}

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

  \item{maxxt}{Matrix or single value defining the maximum
  value for each xt sample.  If a single value is supplied
  then all xt values are given the same maximum value.}

  \item{minxt}{Matrix or single value defining the minimum
  value for each xt sample.  If a single value is supplied
  then all xt values are given the same minimum value}

  \item{maxa}{Vector defining the max value for each
  covariate. IF a single value is supplied then all a
  values are given the same max value}

  \item{mina}{Vector defining the min value for each
  covariate. IF a single value is supplied then all a
  values are given the same max value}

  \item{bUseRandomSearch}{\itemize{ \item \bold{******START
  OF Optimization algorithm SPECIFICATION
  OPTIONS**********}} Use random search (1=TRUE, 0=FALSE)}

  \item{bUseStochasticGradient}{Use Stochastic Gradient
  search (1=TRUE, 0=FALSE)}

  \item{bUseBFGSMinimizer}{Use BFGS Minimizer (1=TRUE,
  0=FALSE)}

  \item{bUseLineSearch}{Use Line search (1=TRUE, 0=FALSE)}

  \item{sgit}{Number of stochastic gradient iterations}

  \item{BFGSConvergenceCriteriaMinStep}{BFGS Minimizer
  Convergence Criteria Minimum Step}

  \item{BFGSProjectedGradientTol}{BFGS Minimizer
  Convergence Criteria Normalized Projected Gradient
  Tolerance}

  \item{BFGSTolerancef}{BFGS Minimizer Line Search
  Tolerance f}

  \item{BFGSToleranceg}{BFGS Minimizer Line Search
  Tolerance g}

  \item{BFGSTolerancex}{BFGS Minimizer Line Search
  Tolerance x}
}
\description{
Optimize the objective function.  The function works for
continuous optimization variables. There are 4 different
optimization algorithms used in this function \enumerate{
\item Adaptive random search. See \code{\link{RS_opt}}.
\item Stochastic gradient. \item A Broyden Fletcher
Goldfarb Shanno (BFGS) method for nonlinear minimization
with box constraints. \item A line search. See
\code{\link{a_line_search}}. } The optimization algorithms
run in series, taking as input the output from the previous
method. The stopping rule used is to test if the line
search algorithm fids a better optimum then its inital
value. If so, then the chain of algorithms is run again.
If line search is not used then the argument
\code{iter_tot} defines the number of times the chain of
algorithms is run. This function takes information from the
PopED database supplied as an argument. The PopED database
supplies information about the the model, parameters,
design and methods to use. Some of the arguments coming
from the PopED database can be overwritten; if they are
supplied then they are used instead of the arguments from
the PopED database.
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
\code{\link{calc_autofocus}};
\code{\link{calc_ofv_and_grad}}; \code{\link{mfea}};
\code{\link{poped_optimize}}
}

