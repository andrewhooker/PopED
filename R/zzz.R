.PopedNamespaceEnv <- new.env()

# Remove CRAN note on no visible binding for global variable
utils::globalVariables(c('.'))

.onAttach <- function(...) {
  
  if (!interactive()) return()
  
  text <- c('  For recent updates to PopED checkout:\n  https://andrewhooker.github.io/PopED/index.html',
            '  Submit suggestions and bug-reports at:\n  https://github.com/andrewhooker/PopED/issues',
            '  Learn more in this introduction to PopED:\n  https://andrewhooker.github.io/PopED/articles/intro-poped.html',
            "  Find out what's changed in PopED at:\n  https://github.com/andrewhooker/PopED/releases.")
  
  packageStartupMessage(sample(text, size = 1))
}
