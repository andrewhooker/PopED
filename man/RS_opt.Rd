\name{RS_opt}
\alias{RS_opt}
\title{Optimize the objective function using an adaptive random search algorithm.}
\usage{
RS_opt(poped.db, ni = NULL, xt = NULL, model_switch = NULL, x = NULL,
  a = NULL, bpopdescr = NULL, ddescr = NULL, maxxt = NULL,
  minxt = NULL, maxa = NULL, mina = NULL, fmf = 0, dmf = 0,
  trflag = TRUE, opt_xt = poped.db$optsw[2], opt_a = poped.db$optsw[4],
  opt_x = poped.db$optsw[3], cfaxt = poped.db$cfaxt, cfaa = poped.db$cfaa,
  rsit = poped.db$rsit, rsit_output = poped.db$rsit_output,
  fim.calc.type = poped.db$iFIMCalculationType,
  ofv_calc_type = poped.db$ofv_calc_type,
  approx_type = poped.db$iApproximationMethod, iter = 1, ...)
}
\arguments{
  \item{cfaxt}{First step factor for sample times}

  \item{...}{arguments passed to
  \code{\link{evaluate.fim}}.}

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
Optimize the objective function using an adaptive random
search algorithm. The function works for both discrete and
continuous optimization variables. This function takes
information from the PopED database supplied as an
argument. The PopED database supplies information about the
the model, parameters, design and methods to use. Some of
the arguments coming from the PopED database can be
overwritten; by default these arguments are \code{NULL} in
the function, if they are supplied then they are used
instead of the arguments from the PopED database.
}
\examples{
\dontrun{
rs.output <- RS_opt(poped.db)
rs.output <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=20,fim.calc.type=0)
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
\code{\link{poped_optimize}}
}

