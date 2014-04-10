\name{plot_efficiency_of_windows}
\alias{plot_efficiency_of_windows}
\title{Plot the efficience of windows}
\usage{
plot_efficiency_of_windows(poped.db, xt_windows, iNumSimulations = 100, ...)
}
\arguments{
  \item{poped.db}{A poped database}

  \item{iNumSimulations}{The number of efficiency
  calculations to make.}

  \item{...}{Extra arguments passed to \code{evaluate.fim}}

  \item{xt_windows}{The distance on one direction from the
  optimal sample times.  Can be a number or a matrix of the
  same size as the xt matrix found in \code{poped.db$gxt}.}
}
\value{
A \link[ggplot2]{ggplot2} object.
}
\description{
Function plots the efficiency of windows around the optimal
design points.
}
\seealso{
Other Graphics: \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{plot_model_prediction}}

Other Simulation: \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{model_prediction}};
\code{\link{plot_model_prediction}}

Other evaluate_design: \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{evaluate.fim}}; \code{\link{get_rse}};
\code{\link{model_prediction}};
\code{\link{plot_model_prediction}}
}

