#' PopED - \bold{Pop}ulation (and individual) optimal \bold{E}xperimental \bold{D}esign. 
#'
#' PopED computes optimal experimental designs for both 
#' population  and individual studies based on nonlinear mixed-effect models.  
#' Often this is based on a computation of the Fisher Information Matrix (FIM). 
#'
#' 
#'
#' @references \enumerate{
#' \item J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. Hooker, "PopED: An extended, 
#' parallelized, nonlinear mixed effects models optimal design tool",  
#' Computer Methods and Programs in Biomedicine, 108, 2012.
#' \item M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a software for optimal 
#' experimental design in population kinetics", Computer Methods and Programs in Biomedicine, 74, 2004. 
#' \item \href{http://poped.sourceforge.net}{poped.sf.net}
#' \item \url{https://github.com/andrewhooker/PopED.git}
#' }
#' 
#' @family FIM 
#' @family Optimize 
#' @family Helper 
#' @family MATLAB
#' @family poped_input 
#' @family matrix_manipulation 
#' @family E-family 
#' @family evaluate_design 
#' @family Simulation 
#' @family Graphics 
#' 
#' @import ggplot2
#' @importFrom MASS write.matrix
#' @importFrom mvtnorm rmvnorm
#' @docType package
#' @name PopED-package
#' @aliases PopED poped
NULL