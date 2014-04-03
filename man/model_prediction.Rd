\name{model_prediction}
\alias{model_prediction}
\title{Model predictions}
\usage{
model_prediction(poped.db, models_to_use = "all", model_num_points = NULL,
  model_minxt = NULL, model_maxxt = NULL, include_sample_times = T,
  groups_to_use = "all", IPRED = FALSE, DV = FALSE, num_ids = 100)
}
\arguments{
  \item{models_to_use}{Which model number should we use?}

  \item{model_num_points}{How many points should be
  plotted.  If not a number then the design in poped.db is
  used.}

  \item{model_minxt}{The minimum of the sample times for
  the predictions.}

  \item{model_maxxt}{The maximum of the sample times for
  the predictions.}

  \item{include_sample_times}{Should the sample times from
  poped.db be included in the predictions?}

  \item{IPRED}{Should we simulate individual predictions?}

  \item{DV}{should we simulate observations?}

  \item{num_ids}{The number of individuals to simulate if
  using IPRED or DV.}

  \item{groups_to_use}{Which groups should we use for
  predictions from the poped.db.}

  \item{poped.db}{A PopED database.}
}
\value{
A dataframe of simulated data, either with some dense grid
of samples or based on the design in the poped database.
}
\description{
Function generates model predictions for the typical value
in the population, individual predictions and data
predictions.
}
\seealso{
Other Simulation: \code{\link{plot_efficiency_of_windows}};
\code{\link{plot_model_prediction}}

Other evaluate_design: \code{\link{evaluate.fim}};
\code{\link{get_rse}};
\code{\link{plot_efficiency_of_windows}};
\code{\link{plot_model_prediction}}
}

