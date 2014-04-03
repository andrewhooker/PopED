\name{ofv_fim}
\alias{ofv_fim}
\title{Evaluate a criterion of the Fisher Information Matrix (FIM)}
\usage{
ofv_fim(fmf, poped.db, ofv_calc_type = poped.db$ofv_calc_type,
  ds_index = poped.db$ds_index)
}
\arguments{
  \item{fmf}{The FIM}

  \item{poped.db}{A poped database}

  \item{ofv_calc_type}{OFV calculation type for FIM
  (1=Determinant of FIM,4=log determinant of
  FIM,6=determinant of interesting part of FIM (Ds))}

  \item{ds_index}{Ds_index is a vector set to 1 if a
  parameter is uninteresting, otherwise 0. size=(1,num
  unfixed parameters). First unfixed bpop, then unfixed d,
  then unfixed docc and last unfixed sigma. Default is the
  fixed effects being important, everything else not
  important.  Must use in conjunction with
  \code{ofv_calc_type=6}.}
}
\value{
The specified criterion value.
}
\description{
Compute a criterion of the FIM given the model, parameters,
design and methods defined in the PopED database.
}
\examples{
\dontrun{
evaluate.fim(poped.db)
}
## Function translated using 'matlab.to.r()'
}
\seealso{
Other FIM: \code{\link{LinMatrixH}};
\code{\link{LinMatrixLH}}; \code{\link{LinMatrixL_occ}};
\code{\link{ed_laplace_ofv}}; \code{\link{ed_mftot}};
\code{\link{evaluate.e.ofv.fim}};
\code{\link{evaluate.fim}}; \code{\link{gradf_eps}};
\code{\link{mf3}}; \code{\link{mf5}}; \code{\link{mf6}};
\code{\link{mf7}}; \code{\link{mf8}}; \code{\link{mftot0}};
\code{\link{mftot1}}; \code{\link{mftot2}};
\code{\link{mftot3}}; \code{\link{mftot4}};
\code{\link{mftot5}}; \code{\link{mftot6}};
\code{\link{mftot7}}; \code{\link{mftot}};
\code{\link{mf}}; \code{\link{ofv_criterion}}
}

