# Start parallel computational processes

This tool chooses the type of parallelization process to use based on
the computer OS being used. For windows the default is "snow" and for
Linux-like systems the default is "multicore"

## Usage

``` r
start_parallel(
  parallel = TRUE,
  num_cores = NULL,
  parallel_type = NULL,
  seed = NULL,
  dlls = NULL,
  mrgsolve_model = NULL,
  ...
)
```

## Arguments

- parallel:

  Should the parallel functionality start up?

- num_cores:

  How many cores to use. Default is
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)-1
  . See [`detectCores`](https://rdrr.io/r/parallel/detectCores.html) for
  more information.

- parallel_type:

  Which type of parallelization should be used? Can be "snow" or
  "multicore". "snow" works on Linux-like systems & Windows. "multicore"
  works only on Linux-like systems. By default this is chosen for you
  depending on your operating system.

- seed:

  The random seed to use.

- dlls:

  If the computations require compiled code (DLL's) and you are using
  the "snow" method then you need to specify the name of the DLL's
  without the extension as a text vector `c("this_file","that_file")`.

- mrgsolve_model:

  If the computations require a mrgsolve model and you are using the
  "snow" method" then you need to specify the name of the model object
  created by `mread` or `mcode`

- ...:

  Arguments passed to
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)

## Value

An atomic vector (TRUE or FALSE) with two attributes: "type" and
"cores".
