# Check dependencies
library(devtools)
load_all('../PopED/')
library(dplyr)

# devtools::install_github("datastorm-open/DependenciesGraphs")
## Downloading GitHub repo datastorm-open/DependenciesGraphs@master
## from URL https://api.github.com/repos/datastorm-open/DependenciesGraphs/zipball/master
## Installing DependenciesGraphs
## Installing 16 packages: data.table, digest, htmltools, htmlwidgets, httpuv, jsonlite, mime, R6, Rcpp, rlist, shiny, shinydashboard, sourcetools, visNetwork, XML, yaml

library(DependenciesGraphs)

# Function to remove helper-functions/nodes that are called by too many different other functions
# to get a better overview
removeNode <- function(deps, nodeStr) {
  IDs = which(deps$Nomfun$label %in% nodeStr)
  rowRemove = deps$fromto$from %in% IDs | deps$fromto$to %in% IDs
  deps$fromto = deps$fromto[!rowRemove,]
  deps$Nomfun = deps$Nomfun[!1:nrow(deps$Nomfun) %in% IDs,]
  return(deps)
}

deps.orig <- envirDependencies("package:PopED")
N = nrow(deps.orig$Nomfun)

# Remove helper functions and display dependency graph
deps2 <- removeNode(deps.orig, nodeStr = c('zeros', 'size', 'isempty', 'diag_matlab', 'inv', 'trace_matrix', 'cell'))
plot(deps2)

# Display which functions are called by no other functions
deps.orig$Nomfun$label[seq(N) %in% setdiff(deps.orig$fromto$to, deps.orig$fromto$from)]
  
# Display which functions which are called but do not call other functions 
deps.orig$Nomfun$label[seq(N) %in% setdiff(deps.orig$fromto$from, deps.orig$fromto$to)]
