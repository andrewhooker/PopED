\name{LinMatrixH}
\alias{LinMatrixH}
\title{Model linearization with respect to epsilon.}
\usage{
LinMatrixH(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, poped.db)
}
\arguments{
  \item{xt_ind}{A vector of the individual/group sample
  times}

  \item{b_ind}{vector of individual realization of the BSV
  terms b}

  \item{bocc_ind}{Vector of individual realizations of the
  BOV terms bocc}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{bpop}{The fixed effects parameter values.  Supplied
  as a vector.}

  \item{poped.db}{A PopED database.}
}
\value{
A matrix of size (samples per individual x number of
epsilons)
}
\description{
The function performs a linearization of the model with
respect to the residual variability. Derivative of model
w.r.t. eps evaluated at eps=0
}
\seealso{
Other FIM: \code{\link{LinMatrixLH}};
\code{\link{LinMatrixL_occ}}; \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{ed_laplace_ofv}}; \code{\link{ed_mftot}};
\code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mftot}};
\code{\link{mf}}; \code{\link{ofv_criterion}};
\code{\link{ofv_fim}}
}

