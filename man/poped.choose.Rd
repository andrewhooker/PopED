\name{poped.choose}
\alias{poped.choose}
\title{Choose between \code{arg1} and \code{arg2}}
\usage{
poped.choose(arg1, arg2)
}
\arguments{
  \item{arg1}{The first argument}

  \item{arg2}{The second argument}
}
\description{
Function chooses \code{arg1} unless it is \code{NULL} in
which case \code{arg2} is chosen.
}
\examples{

poped.choose(2,5)

poped.choose("foo",66)

poped.choose(NULL,"hello")
}
\seealso{
Other poped_input: \code{\link{convert_variables}};
\code{\link{create.poped.database}};
\code{\link{downsizing_general_design}}
}

