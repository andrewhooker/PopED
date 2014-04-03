\name{plot_model_prediction}
\alias{plot_model_prediction}
\title{Plot model predictions}
\usage{
plot_model_prediction(poped.db, model_num_points = 100, separate.groups = F,
  sample.times = T, sample.times.IPRED = F, sample.times.DV = F,
  PRED = T, IPRED = F, IPRED.lines = F, alpha.IPRED.lines = 0.1,
  alpha.IPRED = 0.3, sample.times.size = 4, DV = F, alpha.DV = 0.3,
  DV.lines = F, DV.points = F, alpha.DV.lines = 0.3,
  alpha.DV.points = 0.3, sample.times.DV.points = F,
  sample.times.DV.lines = F, alpha.sample.times.DV.points = 0.3,
  alpha.sample.times.DV.lines = 0.3, y_lab = "Model Predictions",
  facet_scales = "fixed", facet_label_names = T, ...)
}
\arguments{
  \item{separate.groups}{Should there be separate plots for
  each group.}

  \item{sample.times}{Should sample times be shown on the
  plots.}

  \item{sample.times.IPRED}{Should sample times be shown
  based on the IPRED y-values.}

  \item{sample.times.DV}{Should sample times be shown based
  on the DV y-values.}

  \item{PRED}{Should a PRED line be drawn.}

  \item{IPRED.lines}{Should IPRED lines be drawn?}

  \item{alpha.IPRED.lines}{What should the transparency for
  the IPRED.lines be?}

  \item{alpha.IPRED}{What should the tranparency of the
  IPRED CI?}

  \item{sample.times.size}{What should the size of the
  sample.times be?}

  \item{alpha.DV}{What should the tranparency of the DV
  CI?}

  \item{DV.lines}{Should DV lines be drawn?}

  \item{DV.points}{Should DV points be drawn?}

  \item{alpha.DV.lines}{What should the transparency for
  the DV.lines be?}

  \item{alpha.DV.points}{What should the transparency for
  the DV.points be?}

  \item{sample.times.DV.points}{TRUE or FALSE.}

  \item{sample.times.DV.lines}{TRUE or FALSE.}

  \item{alpha.sample.times.DV.points}{What should the
  transparency for the sample.times.DV.points be?}

  \item{alpha.sample.times.DV.lines}{What should the
  transparency for the sample.times.DV.lines be?}

  \item{y_lab}{The label of the y-axis.}

  \item{facet_scales}{Can be "free", "fixed", "free_x" or
  "free_y"}

  \item{facet_label_names}{TRUE or FALSE}

  \item{poped.db}{A PopED database.}

  \item{...}{arguments passed to
  \code{\link{evaluate.fim}}.}

  \item{model_num_points}{How many points should be
  plotted.  If not a number then the design in poped.db is
  used.}

  \item{IPRED}{Should we simulate individual predictions?}

  \item{DV}{should we simulate observations?}
}
\value{
A \link[ggplot2]{ggplot2} object.
}
\description{
Function plots model predictions for the typical value in
the population, individual predictions and data
predictions.
}
\seealso{
Other Graphics: \code{\link{plot_efficiency_of_windows}}

Other Simulation: \code{\link{model_prediction}};
\code{\link{plot_efficiency_of_windows}}

Other evaluate_design: \code{\link{evaluate.fim}};
\code{\link{get_rse}}; \code{\link{model_prediction}};
\code{\link{plot_efficiency_of_windows}}
}

