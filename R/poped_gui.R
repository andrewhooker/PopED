#' Run the graphical interface for PopED
#'
#'
#' @export
#'

poped_gui <- function(){

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  shiny::runApp(system.file('shinyapp', package='PopED'))
  
}