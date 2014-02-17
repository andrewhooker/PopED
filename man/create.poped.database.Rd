\name{create.poped.database}
\alias{create.poped.database}
\title{Create a PopED database from a poped input file}
\usage{
create.poped.database(popedInput, ...)
}
\arguments{
  \item{popedInput}{An input file to PopED}
}
\value{
A PopED database
}
\description{
This function takes the input file supplied by the user and
creates a database that can then be used to run all other
PopED functions.  The function supplies default values to
elements of the database that are not specified in the
input file.  The function also computes needed quatities.
}
\examples{
\dontrun{
source("warfarin.POPED.R")
poped.db <- create.poped.database(warfarin.input())
}
}
\seealso{
Other poped_input: \code{\link{convert_popedInput}};
\code{\link{convert_variables}};
\code{\link{downsizing_general_design}}
}

