PopED: Population (and individual) Experimental Design in R
======

PopED computes optimal experimental designs for both population 
and individual studies based on nonlinear mixed-effect models. 
Often this is based on a computation of the Fisher Information Matrix (FIM). 

## Installation

You need to have R installed.  Download the latest version of R from www.r-project.org.
Install PopED in R using one of the following methods:

* latest stable release -- From CRAN.  Write at the R command line:

```
install.packages("PopED")
```

* Latest development version -- from Github. Note that the command below installs the "master" 
(development) branch; if you want the release branch from Github add `ref="release"` to the
`install_github()` call. The `install_github()` approach requires that you build from source, 
i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; 
you may also need to install dependencies manually.

```
devtools::install_github("PopED",username="andrewhooker")
```

## Getting started

To get started you need to define 

1. A model.
2. An initial design (and design space if you want to optimize). 
3. The tasks to perform.  

There are a number of functions to help you with these tasks.  See `?poped` for more information.  
 
There are several other examples, as r-scripts, in the "examples" folder in the 
PopED installation directory located at:

```
system.file("examples", package="PopED")
```

The same examples are located in the "inst/examples" directory of this repository.

