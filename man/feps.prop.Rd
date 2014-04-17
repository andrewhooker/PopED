\name{feps.prop}
\alias{feps.prop}
\title{RUV model:
Proportional.}
\usage{
feps.prop(model_switch, xt, parameters, epsi, poped.db)
}
\arguments{
  \item{epsi}{A matrix with the same number of rows as the
  \code{xt} vector, colums match the numbers defined in
  this function.}

  \item{model_switch}{a vector of values, the same size as
  \code{xt}, identifying which model response should be
  computed for the corresponding xt value.  Used for
  multiple response models.}

  \item{xt}{a vector of independent variable values (often
  time).}

  \item{parameters}{A named list of parameter values.}

  \item{poped.db}{a poped database.  This can be used to
  extract information that may be needed in the model
  file.}
}
\value{
A list consisting of: \enumerate{ \item y the values of the
model at the specified points. \item poped.db A
(potentially modified) poped database. }
}
\description{
This is a residual unexplained variability (RUV) model
function that encodes the model described above. The
function is suitable for input to the
\code{\link{create.poped.database}} function using the
\code{fError_file} argument.
}
\seealso{
Other RUV_models: \code{\link{feps.add.prop}};
\code{\link{feps.add}}

Other models: \code{\link{feps.add.prop}};
\code{\link{feps.add}};
\code{\link{ff.PK.1.comp.oral.md.CL}};
\code{\link{ff.PK.1.comp.oral.md.KE}};
\code{\link{ff.PK.1.comp.oral.sd.CL}};
\code{\link{ff.PK.1.comp.oral.sd.KE}};
\code{\link{ff.PKPD.1.comp.oral.md.CL.imax}};
\code{\link{ff.PKPD.1.comp.sd.CL.emax}}
}

