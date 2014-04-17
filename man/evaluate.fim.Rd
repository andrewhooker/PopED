\name{evaluate.fim}
\alias{evaluate.fim}
\title{Evaluate the Fisher Information Matrix (FIM)}
\usage{
evaluate.fim(poped.db, fim.calc.type = NULL, bpop.val = NULL,
  d_full = NULL, docc_full = NULL, sigma_full = NULL,
  model_switch = NULL, ni = NULL, xt = NULL, x = NULL, a = NULL,
  groupsize = NULL, deriv.type = NULL, ...)
}
\arguments{
  \item{poped.db}{A PopED database.}

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

  \item{bpop.val}{The fixed effects parameter values.
  Supplied as a vector.}

  \item{d_full}{A between subject variability matrix (OMEGA
  in NONMEM).}

  \item{docc_full}{A between occasion variability matrix.}

  \item{sigma_full}{A residual unexplained variability
  matrix (SIGMA in NONMEM).}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{ni}{A vector of the number of samples in each
  group.}

  \item{xt}{A matrix of sample times.  Each row is a vector
  of sample times for a group.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{groupsize}{A vector of the numer of individuals in
  each group.}

  \item{deriv.type}{A number indicating the type of
  derivative to use: \itemize{ \item 0=Complex difference
  \item 1=Central difference \item 20=Analytic derivative
  (placeholder) \item 30=Automatic differentiation
  (placeholder) }}

  \item{...}{Other arguments passed to the function.}
}
\value{
The FIM.
}
\description{
Compute the FIM given the model, parameters, design and
methods defined in the PopED database. Some of the
arguments coming from the PopED database can be
overwritten; by default these arguments are \code{NULL} in
the function, if they are supplied then they are used
instead of the arguments from the PopED database.
}
\examples{
\dontrun{
evaluate.fim(poped.db)
}
}
\seealso{
Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_laplace_ofv}}; \code{\link{ed_mftot}};
\code{\link{evaluate.e.ofv.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mftot}};
\code{\link{mf}}; \code{\link{ofv_criterion}};
\code{\link{ofv_fim}}

Other evaluate_FIM: \code{\link{evaluate.e.ofv.fim}};
\code{\link{ofv_fim}}

Other evaluate_design: \code{\link{get_rse}};
\code{\link{model_prediction}};
\code{\link{plot_efficiency_of_windows}};
\code{\link{plot_model_prediction}}
}

