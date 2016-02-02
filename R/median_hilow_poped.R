#' Wrap summary functions from Hmisc and ggplot to work with stat_summary in ggplot
#' 
#' Created for back compatibility with older versions of ggplot, and so that 
#' PopED does not have to load ggplot when started.
#' 
#' @param x A numeric vector
#' @param ... Additional arguments passed to \code{\link[Hmisc]{smedian.hilow}} or 
#'   \code{\link[ggplot2]{median_hilow}}, depending on your version of ggplot.
#' 
#' @export


median_hilow_poped <- function (x, ...) 
{
  if(packageVersion("ggplot2") <= "1.0.1"){
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop("Hmisc package needed for this function to work. Please install it.",
           call. = FALSE)
    }
    result <- do.call(Hmisc::smedian.hilow, list(x = x, ...))
    return(dplyr::rename_(data.frame(t(result)),y="Median",ymin="Lower",ymax="Upper"))
  } else {
    return(do.call(ggplot2::median_hilow, list(x = x, ...)))
  }
}