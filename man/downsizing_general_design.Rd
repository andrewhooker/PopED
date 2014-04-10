\name{downsizing_general_design}
\alias{downsizing_general_design}
\title{Downsize a general design to a specific design}
\usage{
downsizing_general_design(globalStructure)
}
\arguments{
  \item{globalStructure}{A PopED database}
}
\value{
A list containing: \item{ni}{A vector of the number of
samples in each group.} \item{xt}{A matrix of sample times.
Each row is a vector of sample times for a group.}
\item{model_switch}{A matrix that is the same size as xt,
specifying which model each sample belongs to.} \item{x}{A
matrix for the discrete design variables.  Each row is a
group.} \item{a}{A matrix of covariates.  Each row is a
group.} \item{bpop}{A matrix of fixed effect parameter
values.}
}
\description{
Function takes a design with potentially empty design
variables and resuces the design so that a FIM can be
calculated.
}
\seealso{
Other poped_input: \code{\link{PopED}},
\code{\link{PopED-package}}, \code{\link{poped}};
\code{\link{convert_variables}};
\code{\link{create.poped.database}};
\code{\link{poped.choose}}
}

