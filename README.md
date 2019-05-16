
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PopED <img src="man/figures/logo.png" align="right" alt="" width="240" />

[![Travis-CI Build
Status](https://travis-ci.org/andrewhooker/PopED.svg?branch=master)](https://travis-ci.org/andrewhooker/PopED)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/andrewhooker/PopED?branch=master&svg=true)](https://ci.appveyor.com/project/andrewhooker/PopED)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/PopED)](https://CRAN.R-project.org/package=PopED)
[![codecov.io](https://codecov.io/github/andrewhooker/PopED/coverage.svg?branch=master)](https://codecov.io/github/andrewhooker/PopED?branch=master)

PopED computes optimal experimental designs for both population and
individual studies based on nonlinear mixed-effect models. Often this is
based on a computation of the Fisher Information Matrix (FIM).

## Installation

You need to have R installed. Download the latest version of R from
www.r-project.org. You can install the released version of PopED from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PopED")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andrewhooker/PopED")
```

## Getting started

To get started you need to define

1.  A model.
2.  An initial design (and design space if you want to optimize).
3.  The tasks to perform.

Learn more in this [introduction to
PopED](https://andrewhooker.github.io/PopED/articles/intro-poped.html)

## Contact

You are welcome to:

  - submit suggestions and bug-reports at:
    <https://github.com/andrewhooker/PopED/issues>
  - send a pull request on: <https://github.com/andrewhooker/PopED>
  - compose a friendly e-mail to: <andrew.hooker@farmbio.uu.se>
