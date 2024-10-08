#' PopED - \bold{Pop}ulation (and individual) optimal \bold{E}xperimental \bold{D}esign. 
#'
#' PopED computes optimal experimental designs for both 
#' population  and individual studies based on nonlinear mixed-effect models.  
#' Often this is based on a computation of the Fisher Information Matrix (FIM). 
#' 
#'
#' To get started you need to define 
#' \enumerate{
#' \item A model.
#' \item An initial design (and design space if you want to optimize). 
#' \item The tasks to perform.  
#' }
#' There are a number of functions to help you with these tasks.  The user-level functions  defined below are 
#' meant to be run with a minimum of arguments (for beginners to advanced users).  Many of the other functions in the package
#' (and not listed here) are called by these user-level functions 
#' and are often not as user 
#' friendly (developer level or advanced user functions).
#' 
#' Define a structural model: 
#' \code{\link{ff.PK.1.comp.oral.md.CL}}, 
#'  \code{\link{ff.PK.1.comp.oral.md.KE}}, 
#'  \code{\link{ff.PK.1.comp.oral.sd.CL}}, 
#'  \code{\link{ff.PK.1.comp.oral.sd.KE}}, 
#'  \code{\link{ff.PKPD.1.comp.oral.md.CL.imax}}, 
#'  \code{\link{ff.PKPD.1.comp.sd.CL.emax}}.
#' 
#' Define a residual unexplained variability model (residual error model): 
#' \code{\link{feps.add.prop}},
#' \code{\link{feps.add}}, 
#' \code{\link{feps.prop}}.
#' 
#' Create an initial study design (and design space): 
#' \code{\link{create.poped.database}}.
#' 
#' Evaluate the model and/or design through simulation and graphics:
#' \code{\link{plot_model_prediction}}, 
#' \code{\link{model_prediction}}, 
#' \code{\link{plot_efficiency_of_windows}}.
#' 
#' Evaluate the design using the FIM:
#' \code{\link{evaluate_design}}, 
#' \code{\link{evaluate.fim}}, 
#' \code{\link{evaluate.e.ofv.fim}}, 
#' \code{\link{ofv_fim}},
#' \code{\link{get_rse}}.
#' 
#' Optimize the design (evaluate afterwards using the above functions): 
#' \code{\link{poped_optim}}, 
#' 
#' See the "Examples" section below for a short introduction to using the above functions. 
#' There are several other examples, as r-scripts, in the "examples" folder in the 
#' PopED installation directory located at (run at the R command line):
#' 
#' \code{system.file("examples", package="PopED")}.
#'
#' @references \enumerate{
#' \item J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. Hooker, "PopED: An extended, 
#' parallelized, nonlinear mixed effects models optimal design tool",  
#' Computer Methods and Programs in Biomedicine, 108, 2012.
#' \item M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a software for optimal 
#' experimental design in population kinetics", Computer Methods and Programs in Biomedicine, 74, 2004. 
#' \item \url{https://andrewhooker.github.io/PopED/}
#' }
#' 
#' 
#' @example tests/testthat/examples_fcn_doc/examples_poped-package.R
#' 
#' @importFrom stats dlnorm
#' @importFrom stats dnorm
#' @importFrom stats end
#' @importFrom stats optim
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom utils capture.output
#' @importFrom utils find
#' @importFrom utils packageVersion
#' @importFrom utils stack
#' @importFrom magrittr "%>%" 
# @importFrom rlang .data
# @importFrom devtools load_all
# @importFrom MASS write.matrix
# @importFrom mvtnorm rmvnorm
# @import ggplot2
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
