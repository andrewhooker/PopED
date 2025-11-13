# Wrap summary functions from Hmisc and ggplot to work with stat_summary in ggplot

Created for back compatibility with older versions of ggplot, and so
that PopED does not have to load ggplot when started.

## Usage

``` r
median_hilow_poped(x, ...)
```

## Arguments

- x:

  A numeric vector

- ...:

  Additional arguments passed to Hmisc's smedian.hilow function or
  ggplot2's median_hilow function, depending on your version of ggplot.
