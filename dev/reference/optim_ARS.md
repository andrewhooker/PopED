# Optimize a function using adaptive random search.

Optimize an objective function using an adaptive random search
algorithm. The function works for both discrete and continuous
optimization parameters and allows for box-constraints and sets of
allowed values.

## Usage

``` r
optim_ARS(
  par,
  fn,
  lower = NULL,
  upper = NULL,
  allowed_values = NULL,
  loc_fac = 4,
  no_bounds_sd = par,
  iter = 400,
  iter_adapt = 50,
  adapt_scale = 1,
  max_run = 200,
  trace = TRUE,
  trace_iter = 5,
  new_par_max_it = 200,
  maximize = F,
  parallel = F,
  parallel_type = NULL,
  num_cores = NULL,
  mrgsolve_model = NULL,
  seed = round(runif(1, 0, 1e+07)),
  allow_replicates = TRUE,
  replicates_index = seq(1, length(par)),
  generator = NULL,
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

- loc_fac:

  Locality factor for determining the standard deviation of the sampling
  distribution around the current position of the parameters. The
  initial standard deviation is normally calculated as
  `(upper - lower)/loc_fac` except in cases when there are no upper or
  lower limits (e.g. when `upper=Inf` or `lower=-Inf`).

- no_bounds_sd:

  The standard deviation of the sampling distribution around the current
  position of the parameters when there are no upper or lower limits
  (e.g. when `upper=Inf` or `lower=-Inf`).

- iter:

  The number of iterations for the algorithm to perform (this is a
  maximum number, it could be less).

- iter_adapt:

  The number of iterations before adapting (shrinking) the parameter
  search space.

- adapt_scale:

  The scale for adapting the size of the sampling distribution. The
  adaptation of the standard deviation of the sampling distribution
  around the current position of the parameters is done after
  `iter_adapt` iteration with no change in the best objective function.
  When adapting, the standard deviation of the sampling distribution is
  calculated as `(upper - lower)/(loc_fac*ff*adapt_scale)` where ff
  starts at 1 and increases by 1 for each adaptation.

- max_run:

  The maximum number of iterations to run without a change in the best
  parameter estimates.

- trace:

  Should the algorithm output results intermittently.

- trace_iter:

  How many iterations between each update to the screen about the result
  of the search.

- new_par_max_it:

  The algorithm randomly chooses samples based on the current best set
  of parameters. If when drawing these samples the new parameter set has
  already been tested then a new draw is performed. After
  `new_par_max_it` draws, with no new parameter sets, then the algorithm
  stops.

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

- generator:

  A user-defined function that generates new parameter sets to try in
  the algorithm. See examples below.

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
[`optim_LS()`](https://andrewhooker.github.io/PopED/dev/reference/optim_LS.md),
[`poped_optim()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim.md),
[`poped_optim_1()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_1.md),
[`poped_optim_2()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_2.md),
[`poped_optim_3()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optim_3.md),
[`poped_optimize()`](https://andrewhooker.github.io/PopED/dev/reference/poped_optimize.md)

## Examples

``` r
## "wild" function , global minimum at about -15.81515
fw <- function(x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80

# optimization with fewer function evaluations compared to SANN
res1 <- optim_ARS(50, fw,lower = -50, upper=100)
#> Initial OFV = 159.001
#> It.   5 | OFV = 83.4526
#> It.  10 | OFV = 76.4507
#> It.  15 | OFV = 76.4507
#> It.  20 | OFV = 76.4507
#> It.  25 | OFV = 76.4507
#> It.  30 | OFV = 70.9847
#> It.  35 | OFV = 70.9847
#> It.  40 | OFV = 70.9847
#> It.  45 | OFV = 70.9847
#> It.  50 | OFV = 70.9847
#> It.  55 | OFV = 70.9847
#> It.  60 | OFV = 70.9847
#> It.  65 | OFV = 70.9847
#> It.  70 | OFV = 70.9847
#> It.  75 | OFV = 70.9847
#> It.  80 | OFV = 70.9847
#> It.  85 | OFV = 70.9847
#> It.  90 | OFV = 70.9847
#> It.  95 | OFV = 70.9847
#> It. 100 | OFV = 70.9847
#> It. 105 | OFV = 68.8229
#> It. 110 | OFV = 68.8229
#> It. 115 | OFV = 68.8229
#> It. 120 | OFV = 68.8229
#> It. 125 | OFV = 68.8229
#> It. 130 | OFV = 68.8229
#> It. 135 | OFV = 68.8229
#> It. 140 | OFV = 68.8229
#> It. 145 | OFV = 68.8229
#> It. 150 | OFV = 68.5656
#> It. 155 | OFV = 68.5656
#> It. 160 | OFV = 68.5656
#> It. 165 | OFV = 68.5656
#> It. 170 | OFV = 68.5656
#> It. 175 | OFV = 68.5656
#> It. 180 | OFV = 68.5656
#> It. 185 | OFV = 68.5656
#> It. 190 | OFV = 68.5656
#> It. 195 | OFV = 68.5656
#> It. 200 | OFV = 68.5656
#> It. 205 | OFV = 68.4446
#> It. 210 | OFV = 68.4446
#> It. 215 | OFV = 68.4446
#> It. 220 | OFV = 68.4446
#> It. 225 | OFV = 68.4446
#> It. 230 | OFV = 68.4446
#> It. 235 | OFV = 68.4446
#> It. 240 | OFV = 68.4446
#> It. 245 | OFV = 68.4446
#> It. 250 | OFV = 68.4446
#> It. 255 | OFV = 68.4446
#> It. 260 | OFV = 68.4446
#> It. 265 | OFV = 68.4446
#> It. 270 | OFV = 68.4446
#> It. 275 | OFV = 68.4446
#> It. 280 | OFV = 68.4446
#> It. 285 | OFV = 68.4446
#> It. 290 | OFV = 68.4446
#> It. 295 | OFV = 68.4446
#> It. 300 | OFV = 68.4446
#> It. 305 | OFV = 68.4446
#> It. 310 | OFV = 68.4446
#> It. 315 | OFV = 68.4446
#> It. 320 | OFV = 68.4446
#> It. 325 | OFV = 68.4446
#> It. 330 | OFV = 68.4446
#> It. 335 | OFV = 68.4446
#> It. 340 | OFV = 68.4446
#> It. 345 | OFV = 68.4446
#> It. 350 | OFV = 68.4446
#> It. 355 | OFV = 68.4446
#> It. 360 | OFV = 68.4446
#> It. 365 | OFV = 68.4446
#> It. 370 | OFV = 68.4446
#> It. 375 | OFV = 68.4446
#> It. 380 | OFV = 68.4446
#> It. 385 | OFV = 68.4446
#> It. 390 | OFV = 68.4446
#> It. 395 | OFV = 67.5762
#> It. 400 | OFV = 67.5762
#> 
#> Total iterations: 400 
#> Elapsed time: 0.172 seconds.
#> 
#> Final OFV =  67.57618 
#> Parameters: -15.35198 
#> 

# often not as good performance when upper and lower bounds are poor
res2 <- optim_ARS(50, fw, lower=-Inf,upper=Inf)
#> Initial OFV = 159.001
#> It.   5 | OFV = 82.9206
#> It.  10 | OFV = 81.2685
#> It.  15 | OFV = 70.2732
#> It.  20 | OFV = 70.2732
#> It.  25 | OFV = 70.2732
#> It.  30 | OFV = 70.2732
#> It.  35 | OFV = 70.2732
#> It.  40 | OFV = 70.2732
#> It.  45 | OFV = 70.2732
#> It.  50 | OFV = 70.2732
#> It.  55 | OFV = 70.2732
#> It.  60 | OFV = 70.2732
#> It.  65 | OFV = 70.2641
#> It.  70 | OFV = 70.2641
#> It.  75 | OFV = 70.2641
#> It.  80 | OFV = 70.2641
#> It.  85 | OFV = 70.2641
#> It.  90 | OFV = 70.2641
#> It.  95 | OFV = 70.2641
#> It. 100 | OFV = 70.2641
#> It. 105 | OFV = 70.2641
#> It. 110 | OFV = 70.2641
#> It. 115 | OFV = 70.2641
#> It. 120 | OFV = 70.2641
#> It. 125 | OFV = 70.2641
#> It. 130 | OFV = 70.2641
#> It. 135 | OFV = 70.2641
#> It. 140 | OFV = 68.4052
#> It. 145 | OFV = 68.4052
#> It. 150 | OFV = 68.4052
#> It. 155 | OFV = 68.4052
#> It. 160 | OFV = 68.4052
#> It. 165 | OFV = 68.4052
#> It. 170 | OFV = 68.4052
#> It. 175 | OFV = 68.4052
#> It. 180 | OFV = 68.4052
#> It. 185 | OFV = 68.4052
#> It. 190 | OFV = 68.4052
#> It. 195 | OFV = 67.6696
#> It. 200 | OFV = 67.6696
#> It. 205 | OFV = 67.6696
#> It. 210 | OFV = 67.6696
#> It. 215 | OFV = 67.6696
#> It. 220 | OFV = 67.6696
#> It. 225 | OFV = 67.6696
#> It. 230 | OFV = 67.6696
#> It. 235 | OFV = 67.6696
#> It. 240 | OFV = 67.6696
#> It. 245 | OFV = 67.6696
#> It. 250 | OFV = 67.6696
#> It. 255 | OFV = 67.6696
#> It. 260 | OFV = 67.6696
#> It. 265 | OFV = 67.6696
#> It. 270 | OFV = 67.6696
#> It. 275 | OFV = 67.6696
#> It. 280 | OFV = 67.6696
#> It. 285 | OFV = 67.6696
#> It. 290 | OFV = 67.6696
#> It. 295 | OFV = 67.6696
#> It. 300 | OFV = 67.6696
#> It. 305 | OFV = 67.6696
#> It. 310 | OFV = 67.6696
#> It. 315 | OFV = 67.6696
#> It. 320 | OFV = 67.6696
#> It. 325 | OFV = 67.6696
#> It. 330 | OFV = 67.6696
#> It. 335 | OFV = 67.6696
#> It. 340 | OFV = 67.6696
#> It. 345 | OFV = 67.6696
#> It. 350 | OFV = 67.6696
#> It. 355 | OFV = 67.6696
#> It. 360 | OFV = 67.6696
#> It. 365 | OFV = 67.6696
#> It. 370 | OFV = 67.6696
#> It. 375 | OFV = 67.6696
#> It. 380 | OFV = 67.6696
#> It. 385 | OFV = 67.6696
#> It. 390 | OFV = 67.6696
#> Maximum number of identical optimal values reached (max_run=200), optimization stopped.
#> 
#> Total iterations: 392 
#> Elapsed time: 0.147 seconds.
#> 
#> Final OFV =  67.6696 
#> Parameters: -15.97183 
#> 

# Only integer values allowed
if (FALSE) { # \dontrun{ 
res_int <- optim_ARS(50, fw, allowed_values = seq(-50,100,by=1))
} # }

if (FALSE) { # \dontrun{ 
  #plot of the function and solutions
  require(graphics)
  plot(fw, -50, 50, n = 1000, main = "Minimizing 'wild function'")
  points(-15.81515, fw(-15.81515), pch = 16, col = "red", cex = 1)
  points(res1$par, res1$ofv, pch = 16, col = "green", cex = 1)
  points(res2$par, res2$ofv, pch = 16, col = "blue", cex = 1)
} # } 

# optim_ARS does not work great for hard to find minima on flat surface:
# Rosenbrock Banana function
# f(x, y) = (a-x)^2 + b(y-x^2)^2
# global minimum at (x, y)=(a, a^2), where f(x, y)=0. 
# Usually a = 1 and b = 100.
if (FALSE) { # \dontrun{ 
  fr <- function(x,a=1,b=100) {   
    x1 <- x[1]
    x2 <- x[2]
    b*(x2 - x1*x1)^2 + (a - x1)^2
  }
  
  res3 <- optim_ARS(c(-1.2,1), fr,lower = -5, upper = 5)
  
  # plot the surface
  x <- seq(-50, 50, length= 30)
  y <- x
  f <- function(x,y){apply(cbind(x,y),1,fr)}
  z <- outer(x, y, f)
  persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue", ticktype="detailed") -> res
  points(trans3d(1, 1, 0, pmat = res), col = 2, pch = 16,cex=2)
  points(trans3d(res3$par[1], res3$par[1], res3$ofv, pmat = res), col = "green", pch = 16,cex=2)
} # }

# box constraints
flb <- function(x){
  p <- length(x) 
  sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) 
}
## 25-dimensional box constrained
#optim(rep(3, 25), flb,lower = rep(2, 25), upper = rep(4, 25),method = "L-BFGS-B") 
res_box <- optim_ARS(rep(3, 25), flb,lower = rep(2, 25), upper = rep(4, 25)) 
#> Initial OFV = 3460
#> It.   5 | OFV = 3460
#> It.  10 | OFV = 3460
#> It.  15 | OFV = 3460
#> It.  20 | OFV = 3110.35
#> It.  25 | OFV = 3110.35
#> It.  30 | OFV = 2769.9
#> It.  35 | OFV = 2508.65
#> It.  40 | OFV = 1914.21
#> It.  45 | OFV = 1886.45
#> It.  50 | OFV = 1886.45
#> It.  55 | OFV = 1886.45
#> It.  60 | OFV = 1886.45
#> It.  65 | OFV = 1763.46
#> It.  70 | OFV = 1763.46
#> It.  75 | OFV = 1763.46
#> It.  80 | OFV = 902.615
#> It.  85 | OFV = 902.615
#> It.  90 | OFV = 902.615
#> It.  95 | OFV = 902.615
#> It. 100 | OFV = 902.615
#> It. 105 | OFV = 902.615
#> It. 110 | OFV = 902.615
#> It. 115 | OFV = 902.615
#> It. 120 | OFV = 902.615
#> It. 125 | OFV = 902.615
#> It. 130 | OFV = 902.615
#> It. 135 | OFV = 902.615
#> It. 140 | OFV = 865.598
#> It. 145 | OFV = 865.598
#> It. 150 | OFV = 772.584
#> It. 155 | OFV = 772.584
#> It. 160 | OFV = 772.584
#> It. 165 | OFV = 772.584
#> It. 170 | OFV = 772.584
#> It. 175 | OFV = 715.613
#> It. 180 | OFV = 699.241
#> It. 185 | OFV = 699.241
#> It. 190 | OFV = 699.241
#> It. 195 | OFV = 699.241
#> It. 200 | OFV = 699.241
#> It. 205 | OFV = 621.414
#> It. 210 | OFV = 621.414
#> It. 215 | OFV = 621.414
#> It. 220 | OFV = 621.414
#> It. 225 | OFV = 621.414
#> It. 230 | OFV = 621.414
#> It. 235 | OFV = 621.414
#> It. 240 | OFV = 606.201
#> It. 245 | OFV = 606.201
#> It. 250 | OFV = 606.201
#> It. 255 | OFV = 606.201
#> It. 260 | OFV = 606.201
#> It. 265 | OFV = 606.201
#> It. 270 | OFV = 606.201
#> It. 275 | OFV = 606.201
#> It. 280 | OFV = 606.201
#> It. 285 | OFV = 606.201
#> It. 290 | OFV = 536.585
#> It. 295 | OFV = 536.585
#> It. 300 | OFV = 536.585
#> It. 305 | OFV = 536.585
#> It. 310 | OFV = 536.585
#> It. 315 | OFV = 536.585
#> It. 320 | OFV = 536.585
#> It. 325 | OFV = 536.585
#> It. 330 | OFV = 536.585
#> It. 335 | OFV = 536.585
#> It. 340 | OFV = 536.585
#> It. 345 | OFV = 536.585
#> It. 350 | OFV = 536.585
#> It. 355 | OFV = 536.585
#> It. 360 | OFV = 536.585
#> It. 365 | OFV = 512.881
#> It. 370 | OFV = 512.881
#> It. 375 | OFV = 506.252
#> It. 380 | OFV = 506.252
#> It. 385 | OFV = 506.252
#> It. 390 | OFV = 506.252
#> It. 395 | OFV = 506.252
#> It. 400 | OFV = 503.399
#> 
#> Total iterations: 400 
#> Elapsed time: 0.172 seconds.
#> 
#> Final OFV =  503.399 
#> Parameters: 2.130566 2 2 2.059427 2 2 2.167358 2 2 2.142739 2 2 2.030765 2.254329 2 2 2.56842 2.045564 2.276299 2 2 2 2 2 3.690364 
#> 


## Combinatorial optimization: Traveling salesman problem
eurodistmat <- as.matrix(eurodist)

distance <- function(sq) {  # Target function
  sq2 <- embed(sq, 2)
  sum(eurodistmat[cbind(sq2[,2], sq2[,1])])
}

genseq <- function(sq) {  # Generate new candidate sequence
  idx <- seq(2, NROW(eurodistmat)-1)
  changepoints <- sample(idx, size = 2, replace = FALSE)
  tmp <- sq[changepoints[1]]
  sq[changepoints[1]] <- sq[changepoints[2]]
  sq[changepoints[2]] <- tmp
  sq
}

sq <- c(1:nrow(eurodistmat), 1)  # Initial sequence: alphabetic
res3 <- optim_ARS(sq,distance,generator=genseq) # Near optimum distance around 12842
#> Initial OFV = 29625
#> It.   5 | OFV = 28885
#> It.  10 | OFV = 28498
#> It.  15 | OFV = 27136
#> It.  20 | OFV = 26242
#> It.  25 | OFV = 25509
#> It.  30 | OFV = 25509
#> It.  35 | OFV = 25509
#> It.  40 | OFV = 25089
#> It.  45 | OFV = 23414
#> It.  50 | OFV = 23414
#> It.  55 | OFV = 23414
#> It.  60 | OFV = 20070
#> It.  65 | OFV = 19725
#> It.  70 | OFV = 19102
#> It.  75 | OFV = 18976
#> It.  80 | OFV = 18976
#> It.  85 | OFV = 18976
#> It.  90 | OFV = 18976
#> It.  95 | OFV = 18410
#> It. 100 | OFV = 18410
#> It. 105 | OFV = 18410
#> It. 110 | OFV = 18410
#> It. 115 | OFV = 18410
#> It. 120 | OFV = 17812
#> It. 125 | OFV = 16702
#> It. 130 | OFV = 16702
#> It. 135 | OFV = 16702
#> It. 140 | OFV = 16665
#> It. 145 | OFV = 16665
#> It. 150 | OFV = 16665
#> It. 155 | OFV = 16665
#> It. 160 | OFV = 16665
#> It. 165 | OFV = 16665
#> It. 170 | OFV = 16665
#> It. 175 | OFV = 16665
#> It. 180 | OFV = 16665
#> It. 185 | OFV = 16665
#> It. 190 | OFV = 16665
#> It. 195 | OFV = 16665
#> It. 200 | OFV = 16665
#> It. 205 | OFV = 16665
#> It. 210 | OFV = 16665
#> It. 215 | OFV = 16665
#> It. 220 | OFV = 16665
#> It. 225 | OFV = 16665
#> It. 230 | OFV = 16665
#> It. 235 | OFV = 16665
#> It. 240 | OFV = 16665
#> It. 245 | OFV = 16665
#> It. 250 | OFV = 16665
#> It. 255 | OFV = 16665
#> It. 260 | OFV = 16665
#> It. 265 | OFV = 16665
#> It. 270 | OFV = 16665
#> It. 275 | OFV = 16665
#> It. 280 | OFV = 16665
#> It. 285 | OFV = 16665
#> It. 290 | OFV = 16665
#> It. 295 | OFV = 16665
#> It. 300 | OFV = 16665
#> It. 305 | OFV = 16665
#> It. 310 | OFV = 16063
#> It. 315 | OFV = 16063
#> It. 320 | OFV = 16063
#> It. 325 | OFV = 16063
#> It. 330 | OFV = 16063
#> It. 335 | OFV = 15999
#> It. 340 | OFV = 15999
#> It. 345 | OFV = 15999
#> It. 350 | OFV = 15999
#> It. 355 | OFV = 15999
#> It. 360 | OFV = 15999
#> It. 365 | OFV = 15999
#> It. 370 | OFV = 15999
#> It. 375 | OFV = 15999
#> It. 380 | OFV = 15999
#> It. 385 | OFV = 15264
#> It. 390 | OFV = 15264
#> It. 395 | OFV = 15264
#> It. 400 | OFV = 15264
#> 
#> Total iterations: 400 
#> Elapsed time: 0.218 seconds.
#> 
#> Final OFV =  15264 
#> Parameters: 1 19 8 5 4 3 18 2 14 9 12 15 13 16 17 6 11 7 20 10 21 1 
#> 

if (FALSE) { # \dontrun{ 
  # plot of initial sequence
  # rotate for conventional orientation
  loc <- -cmdscale(eurodist, add = TRUE)$points
  x <- loc[,1]; y <- loc[,2]
  s <- seq_len(nrow(eurodistmat))
  tspinit <- loc[sq,]
  
  plot(x, y, type = "n", asp = 1, xlab = "", ylab = "",
       main = paste("Initial sequence of traveling salesman problem\n",
                    "Distance =",distance(sq)), axes = FALSE)
  arrows(tspinit[s,1], tspinit[s,2], tspinit[s+1,1], tspinit[s+1,2],
         angle = 10, col = "green")
  text(x, y, labels(eurodist), cex = 0.8)
  
  # plot of final sequence from optim_ARS
  tspres <- loc[res3$par,]
  plot(x, y, type = "n", asp = 1, xlab = "", ylab = "",
       main = paste("optim_ARS() 'solving' traveling salesman problem\n",
                    "Distance =",distance(c(1,res3$par,1))),axes = FALSE)
  arrows(tspres[s,1], tspres[s,2], tspres[s+1,1], tspres[s+1,2],
         angle = 10, col = "red")
  text(x, y, labels(eurodist), cex = 0.8)
  
  # using optim
  set.seed(123) # chosen to get a good soln relatively quickly
  (res4 <- optim(sq, distance, genseq, method = "SANN",
                 control = list(maxit = 30000, temp = 2000, trace = TRUE,
                                REPORT = 500))) 
  
  tspres <- loc[res4$par,]
  plot(x, y, type = "n", asp = 1, xlab = "", ylab = "",
       main = paste("optim() 'solving' traveling salesman problem\n",
                    "Distance =",distance(res4$par)),axes = FALSE)
  arrows(tspres[s,1], tspres[s,2], tspres[s+1,1], tspres[s+1,2],
         angle = 10, col = "red")
  text(x, y, labels(eurodist), cex = 0.8)
} # }  

# one-dimensional function
if (FALSE) { # \dontrun{ 
  f <- function(x)  abs(x)+cos(x)
  res5 <- optim_ARS(-20,f,lower=-20, upper=20)
  
  curve(f, -20, 20)
  abline(v = res5$par, lty = 4,col="green")
} # }  

# one-dimensional function
f <- function(x)  (x^2+x)*cos(x) # -10 < x < 10
res_max <- optim_ARS(0,f,lower=-10, upper=10,maximize=TRUE) # sometimes to local maxima
#> Initial OFV = 0
#> It.   5 | OFV = 23.0891
#> It.  10 | OFV = 23.0891
#> It.  15 | OFV = 34.1068
#> It.  20 | OFV = 34.1068
#> It.  25 | OFV = 34.6902
#> It.  30 | OFV = 34.6902
#> It.  35 | OFV = 34.6902
#> It.  40 | OFV = 34.6902
#> It.  45 | OFV = 34.6902
#> It.  50 | OFV = 34.6902
#> It.  55 | OFV = 34.6902
#> It.  60 | OFV = 34.6902
#> It.  65 | OFV = 34.6902
#> It.  70 | OFV = 34.6902
#> It.  75 | OFV = 35.4986
#> It.  80 | OFV = 46.4501
#> It.  85 | OFV = 46.4501
#> It.  90 | OFV = 46.4501
#> It.  95 | OFV = 46.4501
#> It. 100 | OFV = 46.6179
#> It. 105 | OFV = 46.6179
#> It. 110 | OFV = 47.6976
#> It. 115 | OFV = 47.6976
#> It. 120 | OFV = 47.6976
#> It. 125 | OFV = 47.6976
#> It. 130 | OFV = 47.6976
#> It. 135 | OFV = 47.6976
#> It. 140 | OFV = 47.6976
#> It. 145 | OFV = 47.6976
#> It. 150 | OFV = 47.6976
#> It. 155 | OFV = 47.6976
#> It. 160 | OFV = 47.6976
#> It. 165 | OFV = 47.6976
#> It. 170 | OFV = 47.6976
#> It. 175 | OFV = 47.6976
#> It. 180 | OFV = 47.6976
#> It. 185 | OFV = 47.6976
#> It. 190 | OFV = 47.6976
#> It. 195 | OFV = 47.6976
#> It. 200 | OFV = 47.6976
#> It. 205 | OFV = 47.7032
#> It. 210 | OFV = 47.7055
#> It. 215 | OFV = 47.7055
#> It. 220 | OFV = 47.7055
#> It. 225 | OFV = 47.7055
#> It. 230 | OFV = 47.7055
#> It. 235 | OFV = 47.7055
#> It. 240 | OFV = 47.7055
#> It. 245 | OFV = 47.7055
#> It. 250 | OFV = 47.7055
#> It. 255 | OFV = 47.7055
#> It. 260 | OFV = 47.7055
#> It. 265 | OFV = 47.7055
#> It. 270 | OFV = 47.7055
#> It. 275 | OFV = 47.7055
#> It. 280 | OFV = 47.7055
#> It. 285 | OFV = 47.7055
#> It. 290 | OFV = 47.7055
#> It. 295 | OFV = 47.7055
#> It. 300 | OFV = 47.7055
#> It. 305 | OFV = 47.7055
#> It. 310 | OFV = 47.7055
#> It. 315 | OFV = 47.7055
#> It. 320 | OFV = 47.7055
#> It. 325 | OFV = 47.7055
#> It. 330 | OFV = 47.7055
#> It. 335 | OFV = 47.7055
#> It. 340 | OFV = 47.7055
#> It. 345 | OFV = 47.7055
#> It. 350 | OFV = 47.7055
#> It. 355 | OFV = 47.7055
#> It. 360 | OFV = 47.7055
#> It. 365 | OFV = 47.7055
#> It. 370 | OFV = 47.7055
#> It. 375 | OFV = 47.7055
#> It. 380 | OFV = 47.7055
#> It. 385 | OFV = 47.7055
#> It. 390 | OFV = 47.7055
#> It. 395 | OFV = 47.7055
#> It. 400 | OFV = 47.7055
#> 
#> Total iterations: 400 
#> Elapsed time: 0.165 seconds.
#> 
#> Final OFV =  47.7055 
#> Parameters: 6.562641 
#> 

if (FALSE) { # \dontrun{ 
  res_min <- optim_ARS(0,f,lower=-10, upper=10) # sometimes to local minima
  
  curve(f, -10, 10)
  abline(v = res_min$par, lty = 4,col="green")
  abline(v = res_max$par, lty = 4,col="red")
} # }


# two-dimensional Rastrigin function
#It has a global minimum at f(x) = f(0) = 0.
if (FALSE) { # \dontrun{ 
  Rastrigin <- function(x1, x2){
    20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
  }
  
  
  x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
  z <- outer(x1, x2, Rastrigin)
  
  res6 <- optim_ARS(c(-4,4),function(x) Rastrigin(x[1], x[2]),lower=-5.12, upper=5.12)
  
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
} # }


## Parallel computation  
## works better when each evaluation takes longer
## here we have added extra time to the computations
## just to show that it works
if (FALSE) { # \dontrun{ 
  res7 <- optim_ARS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                    lower=-5.12, upper=5.12)
  res8 <- optim_ARS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                    lower=-5.12, upper=5.12,parallel = T)
  res9 <- optim_ARS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                    lower=-5.12, upper=5.12,parallel = T,parallel_type = "snow")
} # }
```
