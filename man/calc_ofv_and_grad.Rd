\name{calc_ofv_and_grad}
\alias{calc_ofv_and_grad}
\title{Compute an objective function and gradient}
\usage{
calc_ofv_and_grad(x, optxt, opta, model_switch, aa, axt, groupsize, ni, xtopto,
  xopto, aopto, bpop, d, sigma, docc_full, poped.db, only_fim = FALSE)
}
\arguments{
  \item{aa}{The aa value}

  \item{axt}{the axt value}

  \item{xtopto}{the xtopto value}

  \item{xopto}{the xopto value}

  \item{optxt}{If sampling times are optimized}

  \item{opta}{If continuous design variables are optimized}

  \item{aopto}{the aopto value}

  \item{only_fim}{Should the gradient be calculated?}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{docc_full}{A between occasion variability matrix.}

  \item{poped.db}{A PopED database.}

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

  \item{sigma}{Matrix defining the variances can
  covariances of the residual variability terms of the
  model.}
}
\value{
A list containing: \item{f}{The objective function.}
\item{g}{The gradient.}
}
\description{
Compute an objective function and gradient with respect to
the optimization parameters. This function can be passed to
the Broyden Fletcher Goldfarb Shanno (BFGS) method for
nonlinear minimization with box constraints implemented in
\code{\link{bfgsb_min}}.
}
\seealso{
Other Optimization: \code{\link{Doptim}};
\code{\link{RS_opt}}; \code{\link{calc_autofocus}};
\code{\link{mfea}}; \code{\link{poped_optimize}}

Other Optimize: \code{\link{a_line_search}};
\code{\link{bfgsb_min}}
}

