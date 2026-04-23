# Optimize a function using a line search algorithm.

`optim_LS` performs sequential grid search optimization of an arbitrary
function with respect to each of the parameters to be optimized over.
The function works for both discrete and continuous optimization
parameters and allows for box-constraints (by using the upper and lower
function arguments) or sets of allowed values (by using the
allowed_values function argument) for all parameters, or on a parameter
per parameter basis.

## Usage

``` r
optim_LS(
  par,
  fn,
  lower = NULL,
  upper = NULL,
  allowed_values = NULL,
  line_length = 50,
  trace = TRUE,
  maximize = F,
  parallel = F,
  parallel_type = NULL,
  num_cores = NULL,
  mrgsolve_model = NULL,
  seed = round(runif(1, 0, 1e+07)),
  allow_replicates = TRUE,
  replicates_index = seq(1, length(par)),
  ofv_initial = NULL,
  closed_bounds = TRUE,
  ...
)
```

## Arguments

- par:

  A vector of initial values for the parameters to be optimized over.

- fn:

  A function to be minimized (or maximized), with first argument the
  vector of parameters over which minimization is to take place. It
  should return a scalar result.

- lower:

  Lower bounds on the parameters. A vector.

- upper:

  Upper bounds on the parameters. A vector.

- allowed_values:

  A list containing allowed values for each parameter
  `list(par1=c(2,3,4,5,6),par2=c(5,6,7,8))`. A vector containing allowed
  values for all parameters is also allowed `c(2,3,4,5,6)`.

- line_length:

  The number of different parameter values per parameter to evaluate.
  The values are selected as an evenly spaced grid between the upper and
  lower bounds.

- trace:

  Should the algorithm output results intermittently.

- maximize:

  Should the function be maximized? Default is to minimize.

- parallel:

  Should we use parallel computations?

- parallel_type:

  Which type of parallelization should be used? Can be "snow" or
  "multicore". "snow" works on Linux-like systems & Windows. "multicore"
  works only on Linux-like systems. By default this is chosen for you
  depending on your operating system. See
  [`start_parallel`](https://andrewhooker.github.io/PopED/dev/reference/start_parallel.md).

- num_cores:

  The number of cores to use in the parallelization. By default is set
  to the number output from
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).
  See
  [`start_parallel`](https://andrewhooker.github.io/PopED/dev/reference/start_parallel.md).

- mrgsolve_model:

  If the computations require a mrgsolve model and you are using the
  "snow" method then you need to specify the name of the model object
  created by `mread` or `mcode`.

- seed:

  The random seed to use in the algorithm,

- allow_replicates:

  Should the algorithm allow parameters to have the same value?

- replicates_index:

  A vector, the same length as the parameters. If two values are the
  same in this vector then the parameters may not assume the same value
  in the optimization.

- ofv_initial:

  An initial objective function value (OFV). If not NULL then the
  initial design is not evaluated and the OFV value is assumed to be
  this number.

- closed_bounds:

  Are the upper and lower limits open (boundaries not allowed) or closed
  (boundaries allowed) bounds?

- ...:

  Additional arguments passed to `fn` and `start_parallel`.

## References

1.  M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a
    software fir optimal experimental design in population kinetics",
    Computer Methods and Programs in Biomedicine, 74, 2004.

2.  J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and
    A.C. Hooker, "PopED: An extended, parallelized, nonlinear mixed
    effects models optimal design tool", Computer Methods and Programs
    in Biomedicine, 108, 2012.

## See also

Other Optimize:
[`Doptim()`](https://andrewhooker.github.io/PopED/dev/reference/Doptim.md),
[`LEDoptim()`](https://andrewhooker.github.io/PopED/dev/reference/LEDoptim.md),
[`RS_opt()`](https://andrewhooker.github.io/PopED/dev/reference/RS_opt.md),
[`a_line_search()`](https://andrewhooker.github.io/PopED/dev/reference/a_line_search.md),
[`bfgsb_min()`](https://andrewhooker.github.io/PopED/dev/reference/bfgsb_min.md),
[`calc_autofocus()`](https://andrewhooker.github.io/PopED/dev/reference/calc_autofocus.md),
[`calc_ofv_and_grad()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_grad.md),
[`mfea()`](https://andrewhooker.github.io/PopED/dev/reference/mfea.md),
[`optim_ARS()`](https://andrewhooker.github.io/PopED/dev/reference/optim_ARS.md),
[`poped_optim()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim.md),
[`poped_optim_1()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_1.md),
[`poped_optim_2()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_2.md),
[`poped_optim_3()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_3.md),
[`poped_optimize()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optimize.md)

## Examples

``` r
# "wild" function 
fw <- function(x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80

# Global minimum of 67.47 at about -15.81515
(fw_min <- fw(-15.81515))
#> [1] 67.46773

if (interactive()){  
  #plot of the function
  require(graphics)
  plot(fw, -50, 50, n = 10000, main = "Minimizing 'wild function'")
  
  # Known minimum
  points(-15.81515, fw_min, pch = 21, col = "red", cex = 1.5)
} 

# optimization with fewer function evaluations 
# compared to SANN: see examples in '?optim'
res1 <- optim_LS(50, fw,lower = -50, upper=50, line_length = 10000)
#> 
#>    Initial parameters: 50 
#>    Initial OFV: 159.0012 
#> 
#>    Searching parameter 1 
#>      Changed from 50 to -15.8166 ; OFV = 67.485 
#> 
#>    Elapsed time: 0.036 seconds.
#> 
#>    Final OFV =  67.48502 
#>    Parameters: -15.81658 
#> 

if (interactive()){ 
  require(graphics)
  plot(fw, -20, 0, n = 10000, main = "Minimizing 'wild function'")
  
  # Known minimum
  points(-15.81515, fw_min, pch = 21, col = "red", cex = 1.5)
  
  #plot of the optimization
  points(res1$par, res1$ofv, pch = 16, col = "green", cex = 1)
} 

# Upper and lower bounds and line_length should be considered carefully
res2 <- optim_LS(50, fw, lower=-Inf,upper=Inf,line_length = 10000)
#> 
#>    Initial parameters: 50 
#>    Initial OFV: 159.0012 
#> 
#>    Searching parameter 1 
#>      Changed from 50 to -5.0055 ; OFV = 69.8766 
#> 
#>    Elapsed time: 0.19 seconds.
#> 
#>    Final OFV =  69.87659 
#>    Parameters: -5.005501 
#> 

# Only integer values allowed
res_int <- optim_LS(50, fw, allowed_values = seq(-50,50,by=1))
#> 
#>    Initial parameters: 50 
#>    Initial OFV: 159.0012 
#> 
#>    Searching parameter 1 
#>      Changed from 50 to -17 ; OFV = 68.5368 
#> 
#>    Elapsed time: 0.001 seconds.
#> 
#>    Final OFV =  68.53679 
#>    Parameters: -17 
#> 


# Rosenbrock Banana function
# f(x, y) = (a-x)^2 + b*(y-x^2)^2
# global minimum at (x, y)=(a, a^2), where f(x, y)=0. 
# Usually a = 1 and b = 100 so x=1 and y=1
if (interactive()){ 
  fr <- function(x,a=1,b=100) {   
    x1 <- x[1]
    x2 <- x[2]
    b*(x2 - x1*x1)^2 + (a - x1)^2
  }
  
  res3 <- optim_LS(c(-1.2,1), fr,lower = -5, upper = 5, line_length = 1000)

  # plot the surface
  x <- seq(-50, 50, length= 30)
  y <- x
  f <- function(x,y){apply(cbind(x,y),1,fr)}
  z <- outer(x, y, f)
  persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue", ticktype="detailed") -> res
  points(trans3d(1, 1, 0, pmat = res), col = 2, pch = 16,cex=2)
  points(trans3d(res3$par[1], res3$par[1], res3$ofv, pmat = res), col = "green", pch = 16,cex=1.5)
}

# box constraints
flb <- function(x){
  p <- length(x) 
  sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) 
}

## 25-dimensional box constrained
if (interactive()){ 
  optim(rep(3, 25), flb,lower = rep(2, 25), upper = rep(4, 25),method = "L-BFGS-B") 
}
res_box <- optim_LS(rep(3, 25), flb,
                    lower = rep(2, 25), 
                    upper = rep(4, 25),
                    line_length = 1000) 
#> 
#>    Initial parameters: 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 
#>    Initial OFV: 3460 
#> 
#>    Searching parameter 7 
#>      Changed from 3 to 2.14414 ; OFV = 3370.22 
#>    Searching parameter 16 
#>      Changed from 3 to 2.14414 ; OFV = 3280.43 
#>    Searching parameter 2 
#>      Changed from 3 to 2.14414 ; OFV = 3190.65 
#>    Searching parameter 25 
#>      Changed from 3 to 4 ; OFV = 3146.65 
#>    Searching parameter 20 
#>      Changed from 3 to 2.14414 ; OFV = 3056.87 
#>    Searching parameter 11 
#>      Changed from 3 to 2.14414 ; OFV = 2967.09 
#>    Searching parameter 1 
#>      Changed from 3 to 2 ; OFV = 2789.85 
#>    Searching parameter 13 
#>      Changed from 3 to 2.14414 ; OFV = 2700.07 
#>    Searching parameter 18 
#>      Changed from 3 to 2.14414 ; OFV = 2610.29 
#>    Searching parameter 9 
#>      Changed from 3 to 2.14414 ; OFV = 2520.5 
#>    Searching parameter 21 
#>      Changed from 3 to 2 ; OFV = 2397.28 
#>    Searching parameter 4 
#>      Changed from 3 to 2.14414 ; OFV = 2307.5 
#>    Searching parameter 5 
#>      Changed from 3 to 2 ; OFV = 2184.28 
#>    Searching parameter 8 
#>      Changed from 3 to 2 ; OFV = 2026.82 
#>    Searching parameter 12 
#>      Changed from 3 to 2 ; OFV = 1869.37 
#>    Searching parameter 19 
#>      Changed from 3 to 2 ; OFV = 1711.91 
#>    Searching parameter 15 
#>      Changed from 3 to 2 ; OFV = 1589.68 
#>    Searching parameter 10 
#>      Changed from 3 to 2 ; OFV = 1432.22 
#>    Searching parameter 22 
#>      Changed from 3 to 2 ; OFV = 1304.22 
#>    Searching parameter 3 
#>      Changed from 3 to 2 ; OFV = 1146.77 
#>    Searching parameter 6 
#>      Changed from 3 to 2 ; OFV = 984.533 
#>    Searching parameter 17 
#>      Changed from 3 to 2 ; OFV = 827.077 
#>    Searching parameter 23 
#>      Changed from 3 to 2 ; OFV = 699.077 
#>    Searching parameter 24 
#>      Changed from 3 to 2.10811 ; OFV = 610.183 
#>    Searching parameter 14 
#>      Changed from 3 to 2 ; OFV = 446.962 
#> 
#>    Elapsed time: 0.135 seconds.
#> 
#>    Final OFV =  446.9622 
#>    Parameters: 2 2.144144 2 2.144144 2 2 2.144144 2 2.144144 2 2.144144 2 2.144144 2 2 2.144144 2 2.144144 2 2.144144 2 2 2 2.108108 4 
#> 

# one-dimensional function
if (interactive()){ 
  f <- function(x)  abs(x)+cos(x)
  res5 <- optim_LS(-20,f,lower=-20, upper=20, line_length = 500)
  
  curve(f, -20, 20)
  abline(v = res5$par, lty = 4,col="green")
}  

# one-dimensional function
f <- function(x)  (x^2+x)*cos(x) # -10 < x < 10
res_max <- optim_LS(0,f,lower=-10, upper=10,maximize=TRUE,line_length = 1000) 
#> 
#>    Initial parameters: 0 
#>    Initial OFV: 0 
#> 
#>    Searching parameter 1 
#>      Changed from 0 to 6.55656 ; OFV = 47.7052 
#> 
#>    Elapsed time: 0.009 seconds.
#> 
#>    Final OFV =  47.7052 
#>    Parameters: 6.556557 
#> 

if (interactive()){ 
  res_min <- optim_LS(0,f,lower=-10, upper=10, line_length = 1000) 
  
  curve(f, -10, 10)
  abline(v = res_min$par, lty = 4,col="green")
  abline(v = res_max$par, lty = 4,col="red")
}


# two-dimensional Rastrigin function
#It has a global minimum at f(x) = f(0) = 0.
if (interactive()){ 
  Rastrigin <- function(x1, x2){
    20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
  }
  
  
  x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
  z <- outer(x1, x2, Rastrigin)
  
  res6 <- optim_LS(c(-4,4),function(x) Rastrigin(x[1], x[2]),
                   lower=-5.12, upper=5.12, line_length = 1000)
  
  # color scale
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <-
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  # Generate the desired number of colors from this palette
  nbcol <- 100
  color <- jet.colors(nbcol)
  # Compute the z-value at the facet centres
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  # Recode facet z-values into color indices
  facetcol <- cut(zfacet, nbcol)
  persp(x1, x2, z, col = color[facetcol], phi = 30, theta = 30)
  filled.contour(x1, x2, z, color.palette = jet.colors)
}


## Parallel computation  
## works better when each evaluation takes longer
## here we have added extra time to the computations
## just to show that it works
if (interactive()){ 
  res7 <- optim_LS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                   lower=-5.12, upper=5.12, line_length = 200)
  res8 <- optim_LS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                   lower=-5.12, upper=5.12, line_length = 200, parallel = TRUE)
  res9 <- optim_LS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                   lower=-5.12, upper=5.12, line_length = 200, parallel = TRUE, 
                   parallel_type = "snow")
}
```
