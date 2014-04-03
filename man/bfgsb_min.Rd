\name{bfgsb_min}
\alias{bfgsb_min}
\title{Nonlinear minimization using BFGS with box constraints}
\usage{
bfgsb_min(f_name, f_options, x0, l, u, options = list())
}
\arguments{
  \item{f_name}{A function name (as a text string) that
  returns an objective function and the gradient of that
  objective function, in that order. See
  \code{\link{calc_ofv_and_grad}} as used in
  \code{\link{Doptim}}.}

  \item{f_options}{Options for the f_name argument.}

  \item{x0}{the intial values to optimize}

  \item{l}{the lower bounds}

  \item{u}{the upper bounds}

  \item{options}{a list of additional settings arguments}
}
\value{
A list containing: \item{x_k}{The objective function.}
\item{f_k}{The gradient.} \item{B_k}{The hessian.}
}
\description{
This is the implementation of a Broyden Fletcher Goldfarb
Shanno (BFGS) method for nonlinear minimization with box
constraints.
}
\seealso{
Other Optimize: \code{\link{a_line_search}};
\code{\link{calc_ofv_and_grad}}
}

