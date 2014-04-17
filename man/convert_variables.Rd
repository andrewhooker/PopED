\name{convert_variables}
\alias{convert_variables}
\title{Create global variables in the PopED database}
\usage{
convert_variables(poped.db)
}
\arguments{
  \item{poped.db}{A PopED database}
}
\value{
A PopED database
}
\description{
Function takes design variables from input files and
converts them to the global variables needed in PopED.
Typically not used by the user.  Instead use the function
\code{\link{create.poped.database}}.
}
\seealso{
Other poped_input: \code{\link{create.poped.database}};
\code{\link{downsizing_general_design}};
\code{\link{poped.choose}}
}

