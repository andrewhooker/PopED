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
  bUseExchangeAlgorithm = poped.db$bUseExchangeAlgorithm, iter = 1,
  d_switch = poped.db$d_switch, ED_samp_size = poped.db$ED_samp_size,
  bLHS = poped.db$bLHS, ...)
}
\arguments{
  \item{...}{arguments passed to other functions. See
  \code{\link{Doptim}}.}

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

  \item{fmf}{The initial value of the FIM. If set to zero
  then it is computed.}

  \item{dmf}{The inital OFV. If set to zero then it is
  computed.}

  \item{trflag}{Should the optimization be output to the
  screen and to a file?}

  \item{opt_xt}{Should the sample times be optimized?}

  \item{opt_a}{Should the continuous design variables be
  optimized?}

  \item{opt_x}{Should the discrete design variables be
  optimized?}

  \item{cfaxt}{First step factor for sample times}

  \item{cfaa}{Stochastic Gradient search first step factor
  for covariates}

  \item{rsit}{Number of Random search iterations}

  \item{rsit_output}{Number of iterations in random search
  between screen output}

  \item{fim.calc.type}{The method used for calculating the
  FIM. Potential values: \itemize{ \item 0 = Full FIM.  No
  assumption that fixed and random effects are
  uncorrelated.  See \code{\link{mftot0}}. \item 1 =
  Reduced FIM. Assume that there is no correlation in the
  FIM between the fixed and random effects, and set these
  elements in the FIM to zero. See \code{\link{mftot1}}.
  \item 2 = weighted models (placeholder). \item 3 = Not
  currently used. \item 4 = Reduced FIM and computing all
  derivatives with respect to the standard deviation of the
  residual unexplained variation (sqrt(SIGMA) in NONMEM).
  This matches what is done in PFIM, and assumes that the
  standard deviation of the residual unexplained variation
  is the estimated parameter (NOTE: NONMEM estimates the
  variance of the resudual unexplained variation by
  default). See \code{\link{mftot4}}. \item 5 = Full FIM
  parameterized with A,B,C matrices & derivative of
  variance. See \code{\link{mftot5}}. \item 6 = Calculate
  one model switch at a time, good for large matrices. See
  \code{\link{mftot6}}. \item 7 = Reduced FIM parameterized
  with A,B,C matrices & derivative of variance See
  \code{\link{mftot7}}. }}

  \item{ofv_calc_type}{OFV calculation type for FIM
  (1=Determinant of FIM,4=log determinant of
  FIM,6=determinant of interesting part of FIM (Ds))}

  \item{approx_type}{Approximation method for model, 0=FO,
  1=FOCE, 2=FOCEI, 3=FOI.}

  \item{iter}{The number of iterations entered into the
  \code{blockheader_2} function.}

  \item{bpop}{Matrix defining the fixed effects, per row
  (row number = parameter_number) we should have: \itemize{
  \item column 1 the type of the distribution for E-family
  designs (0 = Fixed, 1 = Normal, 2 = Uniform, 3 = User
  Defined Distribution, 4 = lognormal and 5 = truncated
  normal) \item column 2 defines the mean. \item column 3
  defines the variance of the distribution (or length of
  uniform distribution). } Can also just supply the
  parameter values as a vector \code{c()}}

  \item{d}{Matrix defining the diagnonals of the IIV (same
  logic as for the fixed efects). can also just supply the
  parameter values as a \code{c()}.}

  \item{bUseExchangeAlgorithm}{Use Exchange algorithm
  (1=TRUE, 0=FALSE)}

  \item{d_switch}{\itemize{ \item \bold{******START OF
  CRITERION SPECIFICATION OPTIONS**********}} D-family
  design (1) or ED-familty design (0) (with or without
  parameter uncertainty)}

  \item{ED_samp_size}{Sample size for E-family sampling}

  \item{bLHS}{How to sample from distributions in E-family
  calculations. 0=Random Sampling, 1=LatinHyperCube --}

  \item{opt_samps}{Are the nuber of sample times per group
  being optimized?}

  \item{opt_inds}{Are the nuber of individuals per group
  being optimized?}
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
\code{\link{RS_opt}}; \code{\link{calc_autofocus}};
\code{\link{calc_ofv_and_grad}}; \code{\link{mfea}}
}

