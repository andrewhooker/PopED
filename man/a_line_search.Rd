% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a_line_search.R
\name{a_line_search}
\alias{a_line_search}
\title{Optimize using line search}
\usage{
a_line_search(
  poped.db,
  out_file = "",
  bED = FALSE,
  diff = 0,
  fmf_initial = 0,
  dmf_initial = 0,
  opt_xt = poped.db$settings$optsw[2],
  opt_a = poped.db$settings$optsw[4],
  opt_x = poped.db$settings$optsw[3],
  opt_samps = poped.db$settings$optsw[1],
  opt_inds = poped.db$settings$optsw[5],
  ls_step_size = poped.db$settings$ls_step_size
)
}
\arguments{
\item{poped.db}{A PopED database.}

\item{out_file}{The output file to write to.}

\item{bED}{If the algorithm should use E-family methods. Logical.}

\item{diff}{The OFV difference that is deemed significant for changing a design.  If,
by changing a design variable the difference between the new and old OFV is less than \code{diff} 
the change is not made.}

\item{fmf_initial}{The initial value of the FIM. If \code{0} then the FIM is calculated from poped.db.}

\item{dmf_initial}{The initial value of the objective function value (OFV). 
If \code{0} then the OFV is calculated from poped.db.}

\item{opt_xt}{Should the sample times be optimized?}

\item{opt_a}{Should the continuous design variables be optimized?}

\item{opt_x}{Should the discrete design variables be optimized?}

\item{opt_samps}{Are the number of sample times per group being optimized?}

\item{opt_inds}{Are the number of individuals per group being optimized?}

\item{ls_step_size}{Number of grid points in the line search.}
}
\value{
A list containing:
\item{fmf}{The FIM.}
\item{dmf}{The final value of the objective function value.}
\item{best_changed}{If the algorithm has found a better design than the starting design.}
\item{xt}{A matrix of sample times.  Each row is a vector of sample times for a group.}
\item{x}{A matrix for the discrete design variables.  Each row is a group.}
\item{a}{A matrix of covariates.  Each row is a group.}
\item{poped.db}{A PopED database.}
}
\description{
The function performs a grid search sequentially along 
design variables.  The grid is defined by ls_step_size.
}
\examples{
library(PopED)

############# START #################
## Create PopED database
## (warfarin model for optimization)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error  
## to avoid sample times at very low concentrations (time 0 or very late samples).

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
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.add.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(prop=0.01,add=0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0.01,
                                  maxxt=120,
                                  a=c(DOSE=70),
                                  mina=c(DOSE=0.01),
                                  maxa=c(DOSE=100))

############# END ###################
## Create PopED database
## (warfarin model for optimization)
#####################################


# very sparse grid to evaluate (4 points for each design valiable)
output <- a_line_search(poped.db, opt_xt=TRUE, opt_a=TRUE, ls_step_size=4)

\dontrun{  
  
  # longer run time
  output <- a_line_search(poped.db,opt_xt=TRUE)
  
  # output to a text file
  output <- a_line_search(poped.db,opt_xt=TRUE,out_file="tmp.txt")
  
}

}
\seealso{
Other Optimize: 
\code{\link{Doptim}()},
\code{\link{LEDoptim}()},
\code{\link{RS_opt}()},
\code{\link{bfgsb_min}()},
\code{\link{calc_autofocus}()},
\code{\link{calc_ofv_and_grad}()},
\code{\link{mfea}()},
\code{\link{optim_ARS}()},
\code{\link{optim_LS}()},
\code{\link{poped_optim}()},
\code{\link{poped_optim_1}()},
\code{\link{poped_optim_2}()},
\code{\link{poped_optim_3}()},
\code{\link{poped_optimize}()}
}
\concept{Optimize}
