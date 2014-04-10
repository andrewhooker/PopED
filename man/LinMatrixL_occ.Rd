\name{LinMatrixL_occ}
\alias{LinMatrixL_occ}
\title{Model linearization with respect to occasion variablity parameters.}
\usage{
LinMatrixL_occ(model_switch, xt_ind, x, a, bpop, b_ind, bocc_ind, iCurrentOcc,
  globalStructure)
}
\arguments{
  \item{iCurrentOcc}{The current occasion.}

  \item{model_switch}{A matrix that is the same size as xt,
  specifying which model each sample belongs to.}

  \item{x}{A matrix for the discrete design variables.
  Each row is a group.}

  \item{a}{A matrix of covariates.  Each row is a group.}

  \item{bpop}{The fixed effects parameter values.  Supplied
  as a vector.}

  \item{globalStructure}{A PopED database.}

  \item{xt_ind}{A vector of the individual/group sample
  times}

  \item{b_ind}{vector of individual realization of the BSV
  terms b}

  \item{bocc_ind}{Vector of individual realizations of the
  BOV terms bocc}
}
\value{
A matrix of size (samples per individual x number of iovs)
}
\description{
The function performs a linearization of the model with
respect to the occation variability parameter.. Derivative
of model w.r.t. eta_occ, evaluated bocc_ind.
}
\seealso{
Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{PopED}},
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

