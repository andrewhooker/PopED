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
\examples{
## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.sd.CL

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                  fg_file="sfg",
                                  fError_file="feps.prop",
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=0.01,
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70)
}
\seealso{
Other poped_input: \code{\link{create.poped.database}};
\code{\link{downsizing_general_design}};
\code{\link{poped.choose}}
}

