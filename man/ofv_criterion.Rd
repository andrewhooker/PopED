\name{ofv_criterion}
\alias{ofv_criterion}
\title{Normalize an objective function by the size of the FIM matrix}
\usage{
ofv_criterion(ofv_f, num_parameters, poped.db)
}
\arguments{
  \item{ofv_f}{An objective function}

  \item{num_parameters}{The number of parameters to use for
  normalization}

  \item{poped.db}{a poped database}
}
\value{
The specified criterion value.
}
\description{
Compute a criterion of the FIM given the non-normalized
FIM, model, parameters, design and methods defined in the
PopED database.
}
\examples{
\dontrun{
evaluate.fim(poped.db)
}
}
\seealso{
Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{PopED}}, \code{\link{PopED-package}},
\code{\link{poped}}; \code{\link{ed_laplace_ofv}};
\code{\link{ed_mftot}}; \code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mftot}};
\code{\link{mf}}; \code{\link{ofv_fim}}
}

