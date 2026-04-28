# Model predictions

Function generates a data frame of model predictions for the typical
value in the population, individual predictions and data predictions.
The function can also be used to generate datasets without predictions
using the design specified in the arguments.

## Usage

``` r
model_prediction(
  poped.db = NULL,
  design = list(xt = poped.db$design[["xt"]], groupsize = poped.db$design$groupsize, m =
    poped.db$design[["m"]], x = poped.db$design[["x"]], a = poped.db$design[["a"]], ni =
    poped.db$design$ni, model_switch = poped.db$design$model_switch),
  model = list(fg_pointer = poped.db$model$fg_pointer, ff_pointer =
    poped.db$model$ff_pointer, ferror_pointer = poped.db$model$ferror_pointer),
  parameters = list(docc = poped.db$parameters$docc, d = poped.db$parameters$d, bpop =
    poped.db$parameters$bpop, covd = poped.db$parameters$covd, covdocc =
    poped.db$parameters$covdocc, sigma = poped.db$parameters$sigma),
  IPRED = FALSE,
  DV = FALSE,
  dosing = NULL,
  predictions = NULL,
  filename = NULL,
  models_to_use = "all",
  model_num_points = NULL,
  model_minxt = NULL,
  model_maxxt = NULL,
  include_sample_times = T,
  groups_to_use = "all",
  include_a = TRUE,
  include_x = TRUE,
  manipulation = NULL,
  PI = FALSE,
  PI_conf_level = 0.95,
  PI_ln_dist = TRUE
)
```

## Arguments

- poped.db:

  A PopED database created by
  [`create.poped.database`](https://andrewhooker.github.io/PopED/dev/reference/create.poped.database.md).

- design:

  A list that is passed as arguments to the function
  [`create_design`](https://andrewhooker.github.io/PopED/dev/reference/create_design.md)
  to create a design object.

- model:

  A list containing the model elements to use for the predictions

- parameters:

  A list of parameters to use in the model predictions.

- IPRED:

  Should we simulate individual predictions?

- DV:

  should we simulate observations?

- dosing:

  A list of lists that adds dosing records to the data frame (Each inner
  list corresponding to a group in the design).

- predictions:

  Should the resulting data frame have predictions? Either `TRUE` or
  `FALSE` or `NULL` in which case the function decides based on other
  arguments.

- filename:

  A filename that the data frame should be written to in comma separate
  value (csv) format.

- models_to_use:

  Which model numbers should we use? Model numbers are defined in
  `design` below using `model_switch`. For an explanation see
  [`create_design`](https://andrewhooker.github.io/PopED/dev/reference/create_design.md).

- model_num_points:

  How many extra observation rows should be created in the data frame
  for each group or individual per model. If used then the points are
  placed evenly between `model_minxt` and `model_maxxt`. This option is
  used by
  [`plot_model_prediction`](https://andrewhooker.github.io/PopED/dev/reference/plot_model_prediction.md)
  to simulate the response of the model on a finer grid then the defined
  design. If `NULL` then only the input design is used. Can be a single
  value or a vector the same length as the number of models.

- model_minxt:

  The minimum time value for extra observation rows indicated by
  `model_num_points`. A vector the same length as the number of models

- model_maxxt:

  The minimum time value for extra observation rows indicated by
  `model_num_points`. A vector the same length as the number of models

- include_sample_times:

  Should observations rows in the output data frame include the times
  indicated in the input design?

- groups_to_use:

  Which groups should we include in the output data frame?Allowed values
  are `"all"` or a vector of numbers indicating the groups to include,
  e.g. `c(1,3,6)`.

- include_a:

  Should we include the continuous design variables in the output?

- include_x:

  Should we include the discrete design variables in the output?

- manipulation:

  A list of one or more
  [`expression`](https://rdrr.io/r/base/expression.html) arguments. Each
  expression is evaluated using the code
  `for(i in 1:length(manipulation)){df <- within(df,{eval(manipulation[[i]])})}`.
  Can be used to transform or create new columns in the resulting data
  frame. Note that these transformations are created after any model
  predictions occur, so transformations in columns having to do with
  input to model predictions will not affect the predictions.

- PI:

  Compute prediction intervals for the data given the model. Predictions
  are based on first-order approximations to the model variance and a
  log-normality assumption of that variance (by default), if all
  predictions are positive, otherwise a normal distribution is assumed.

- PI_conf_level:

  The confidence level for the prediction interval computed.

- PI_ln_dist:

  Should the PI calculation be done assuming log-normal or a normal
  distribution. TRUE is the default and indicates a log-normal
  distribution. If any of the PRED values from the model are negative
  then a normal distribution is assumed.

## Value

A dataframe containing a design and (potentially) simulated data with
some dense grid of samples and/or based on the input design.

## See also

Other evaluate_design:
[`evaluate.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md),
[`evaluate_design()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate_design.md),
[`evaluate_power()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate_power.md),
[`get_rse()`](https://andrewhooker.github.io/PopED/dev/reference/get_rse.md),
[`plot_efficiency_of_windows()`](https://andrewhooker.github.io/PopED/dev/reference/plot_efficiency_of_windows.md),
[`plot_model_prediction()`](https://andrewhooker.github.io/PopED/dev/reference/plot_model_prediction.md)

Other Simulation:
[`plot_efficiency_of_windows()`](https://andrewhooker.github.io/PopED/dev/reference/plot_efficiency_of_windows.md),
[`plot_model_prediction()`](https://andrewhooker.github.io/PopED/dev/reference/plot_model_prediction.md)

## Examples

``` r
## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.md.CL
#> function (model_switch, xt, parameters, poped.db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         N = floor(xt/TAU) + 1
#>         y = (DOSE * Favail/V) * (KA/(KA - CL/V)) * (exp(-CL/V * 
#>             (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - 
#>             exp(-CL/V * TAU)) - exp(-KA * (xt - (N - 1) * TAU)) * 
#>             (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
#>         return(list(y = y, poped.db = poped.db))
#>     })
#> }
#> <bytecode: 0x55bedb304e30>
#> <environment: namespace:PopED>

## -- parameter definition function 
## -- names match parameters in function ff
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
    return(parameters) 
}

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.prop,
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=0.01,
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70)

## data frame with model predictions
df_1 <- model_prediction(poped.db)
head(df_1,n=20)
#>    Time      PRED Group Model a_i
#> 1   0.5 3.4254357     1     1  70
#> 2   1.0 5.4711041     1     1  70
#> 3   2.0 7.3821834     1     1  70
#> 4   6.0 7.9462805     1     1  70
#> 5  24.0 5.6858561     1     1  70
#> 6  36.0 4.5402483     1     1  70
#> 7  72.0 2.3116966     1     1  70
#> 8 120.0 0.9398657     1     1  70

##  data frame with variability 
df_2 <- model_prediction(poped.db,DV=TRUE)
head(df_2,n=20)
#>    ID  Time         DV    IPRED      PRED Group Model a_i
#> 1   1   0.5  2.3800236 2.352885 3.4254357     1     1  70
#> 2   1   1.0  4.2029869 4.167358 5.4711041     1     1  70
#> 3   1   2.0  6.0250484 6.628460 7.3821834     1     1  70
#> 4   1   6.0 10.6003247 9.428073 7.9462805     1     1  70
#> 5   1  24.0  7.3829039 7.111946 5.6858561     1     1  70
#> 6   1  36.0  5.1180355 5.665993 4.5402483     1     1  70
#> 7   1  72.0  2.6919875 2.865000 2.3116966     1     1  70
#> 8   1 120.0  1.1947993 1.154134 0.9398657     1     1  70
#> 9   2   0.5  2.5705709 2.416020 3.4254357     1     1  70
#> 10  2   1.0  4.7509404 4.166896 5.4711041     1     1  70
#> 11  2   2.0  5.7171194 6.334689 7.3821834     1     1  70
#> 12  2   6.0  7.6625930 8.183819 7.9462805     1     1  70
#> 13  2  24.0  6.5161958 6.013418 5.6858561     1     1  70
#> 14  2  36.0  4.9060740 4.806975 4.5402483     1     1  70
#> 15  2  72.0  2.4332146 2.455401 2.3116966     1     1  70
#> 16  2 120.0  0.9933456 1.002590 0.9398657     1     1  70
#> 17  3   0.5  5.4419156 5.831451 3.4254357     1     1  70
#> 18  3   1.0  8.0770953 8.064518 5.4711041     1     1  70
#> 19  3   2.0 10.1800510 9.125600 7.3821834     1     1  70
#> 20  3   6.0  8.6966384 8.472325 7.9462805     1     1  70

## data frame with variability (only IPRED, no DV)
df_3 <- model_prediction(poped.db,IPRED=TRUE)
head(df_3,n=20)
#>    ID  Time    IPRED      PRED Group Model a_i
#> 1   1   0.5 1.331747 3.4254357     1     1  70
#> 2   1   1.0 2.425106 5.4711041     1     1  70
#> 3   1   2.0 4.053278 7.3821834     1     1  70
#> 4   1   6.0 6.564106 7.9462805     1     1  70
#> 5   1  24.0 5.606070 5.6858561     1     1  70
#> 6   1  36.0 4.652204 4.5402483     1     1  70
#> 7   1  72.0 2.657258 2.3116966     1     1  70
#> 8   1 120.0 1.259309 0.9398657     1     1  70
#> 9   2   0.5 5.041980 3.4254357     1     1  70
#> 10  2   1.0 7.300742 5.4711041     1     1  70
#> 11  2   2.0 8.723086 7.3821834     1     1  70
#> 12  2   6.0 8.731333 7.9462805     1     1  70
#> 13  2  24.0 7.192552 5.6858561     1     1  70
#> 14  2  36.0 6.320159 4.5402483     1     1  70
#> 15  2  72.0 4.288082 2.3116966     1     1  70
#> 16  2 120.0 2.556484 0.9398657     1     1  70
#> 17  3   0.5 1.373968 3.4254357     1     1  70
#> 18  3   1.0 2.497014 5.4711041     1     1  70
#> 19  3   2.0 4.156027 7.3821834     1     1  70
#> 20  3   6.0 6.602481 7.9462805     1     1  70

## data frame with model predictions, no continuous design variables in data frame
df_4 <- model_prediction(poped.db,include_a = FALSE)
head(df_4,n=20)
#>    Time      PRED Group Model
#> 1   0.5 3.4254357     1     1
#> 2   1.0 5.4711041     1     1
#> 3   2.0 7.3821834     1     1
#> 4   6.0 7.9462805     1     1
#> 5  24.0 5.6858561     1     1
#> 6  36.0 4.5402483     1     1
#> 7  72.0 2.3116966     1     1
#> 8 120.0 0.9398657     1     1

## -- 2 groups
poped.db.2 <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                    fg_fun=sfg,
                                    fError_fun=feps.prop,
                                    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6), 
                                    sigma=0.01,
                                    groupsize=rbind(3,3),
                                    m=2,
                                    xt=c( 0.5,1,2,6,24,36,72,120),
                                    minxt=0,
                                    maxxt=120,
                                    a=rbind(70,50))

df_5 <- model_prediction(poped.db.2,DV=TRUE)
head(df_5,n=20)
#>    ID  Time        DV     IPRED      PRED Group Model a_i
#> 1   1   0.5 1.8824573 2.6986425 3.4254357     1     1  70
#> 2   1   1.0 4.4606407 4.5224462 5.4711041     1     1  70
#> 3   1   2.0 6.1021726 6.5605583 7.3821834     1     1  70
#> 4   1   6.0 8.0079320 7.7839864 7.9462805     1     1  70
#> 5   1  24.0 5.9242971 5.5574735 5.6858561     1     1  70
#> 6   1  36.0 4.1273096 4.4035807 4.5402483     1     1  70
#> 7   1  72.0 2.2569618 2.1907440 2.3116966     1     1  70
#> 8   1 120.0 0.7775830 0.8635865 0.9398657     1     1  70
#> 9   2   0.5 1.3147314 1.2950745 3.4254357     1     1  70
#> 10  2   1.0 2.3422395 2.3370827 5.4711041     1     1  70
#> 11  2   2.0 3.9529858 3.8424651 7.3821834     1     1  70
#> 12  2   6.0 6.6757835 5.9423279 7.9462805     1     1  70
#> 13  2  24.0 6.2962206 4.7862925 5.6858561     1     1  70
#> 14  2  36.0 3.9901656 3.8830024 4.5402483     1     1  70
#> 15  2  72.0 1.9881554 2.0728924 2.3116966     1     1  70
#> 16  2 120.0 0.8062582 0.8976806 0.9398657     1     1  70
#> 17  3   0.5 2.6556066 2.6104406 3.4254357     1     1  70
#> 18  3   1.0 4.9511436 4.5376033 5.4711041     1     1  70
#> 19  3   2.0 5.9386580 6.9882512 7.3821834     1     1  70
#> 20  3   6.0 9.0706361 9.2475260 7.9462805     1     1  70

## without a poped database, just describing the design
## Useful for creating datasets for use in other software (like NONMEM)
design_1 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  m=2,
  groupsize=3)

design_2 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  m=2,
  groupsize=3,
  a=c(WT=70,AGE=50))

design_3 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  m=2,
  groupsize=3,
  a=list(c(WT=70,AGE=50),c(AGE=45,WT=60)))

(df_6 <- model_prediction(design=design_1))
#>     Time PRED Group Model
#> 1    0.5   NA     1     1
#> 2    1.0   NA     1     1
#> 3    2.0   NA     1     1
#> 4    6.0   NA     1     1
#> 5   24.0   NA     1     1
#> 6   36.0   NA     1     1
#> 7   72.0   NA     1     1
#> 8  120.0   NA     1     1
#> 9    0.5   NA     2     1
#> 10   1.0   NA     2     1
#> 11   2.0   NA     2     1
#> 12   6.0   NA     2     1
#> 13  24.0   NA     2     1
#> 14  36.0   NA     2     1
#> 15  72.0   NA     2     1
#> 16 120.0   NA     2     1
(df_7 <- model_prediction(design=design_2))
#>     Time PRED Group Model WT AGE
#> 1    0.5   NA     1     1 70  50
#> 2    1.0   NA     1     1 70  50
#> 3    2.0   NA     1     1 70  50
#> 4    6.0   NA     1     1 70  50
#> 5   24.0   NA     1     1 70  50
#> 6   36.0   NA     1     1 70  50
#> 7   72.0   NA     1     1 70  50
#> 8  120.0   NA     1     1 70  50
#> 9    0.5   NA     2     1 70  50
#> 10   1.0   NA     2     1 70  50
#> 11   2.0   NA     2     1 70  50
#> 12   6.0   NA     2     1 70  50
#> 13  24.0   NA     2     1 70  50
#> 14  36.0   NA     2     1 70  50
#> 15  72.0   NA     2     1 70  50
#> 16 120.0   NA     2     1 70  50
(df_8 <- model_prediction(design=design_3))
#>     Time PRED Group Model WT AGE
#> 1    0.5   NA     1     1 70  50
#> 2    1.0   NA     1     1 70  50
#> 3    2.0   NA     1     1 70  50
#> 4    6.0   NA     1     1 70  50
#> 5   24.0   NA     1     1 70  50
#> 6   36.0   NA     1     1 70  50
#> 7   72.0   NA     1     1 70  50
#> 8  120.0   NA     1     1 70  50
#> 9    0.5   NA     2     1 60  45
#> 10   1.0   NA     2     1 60  45
#> 11   2.0   NA     2     1 60  45
#> 12   6.0   NA     2     1 60  45
#> 13  24.0   NA     2     1 60  45
#> 14  36.0   NA     2     1 60  45
#> 15  72.0   NA     2     1 60  45
#> 16 120.0   NA     2     1 60  45
(df_9 <- model_prediction(design=design_3,DV=TRUE))
#>    ID  Time DV IPRED PRED Group Model WT AGE
#> 1   1   0.5 NA    NA   NA     1     1 70  50
#> 2   1   1.0 NA    NA   NA     1     1 70  50
#> 3   1   2.0 NA    NA   NA     1     1 70  50
#> 4   1   6.0 NA    NA   NA     1     1 70  50
#> 5   1  24.0 NA    NA   NA     1     1 70  50
#> 6   1  36.0 NA    NA   NA     1     1 70  50
#> 7   1  72.0 NA    NA   NA     1     1 70  50
#> 8   1 120.0 NA    NA   NA     1     1 70  50
#> 9   2   0.5 NA    NA   NA     1     1 70  50
#> 10  2   1.0 NA    NA   NA     1     1 70  50
#> 11  2   2.0 NA    NA   NA     1     1 70  50
#> 12  2   6.0 NA    NA   NA     1     1 70  50
#> 13  2  24.0 NA    NA   NA     1     1 70  50
#> 14  2  36.0 NA    NA   NA     1     1 70  50
#> 15  2  72.0 NA    NA   NA     1     1 70  50
#> 16  2 120.0 NA    NA   NA     1     1 70  50
#> 17  3   0.5 NA    NA   NA     1     1 70  50
#> 18  3   1.0 NA    NA   NA     1     1 70  50
#> 19  3   2.0 NA    NA   NA     1     1 70  50
#> 20  3   6.0 NA    NA   NA     1     1 70  50
#> 21  3  24.0 NA    NA   NA     1     1 70  50
#> 22  3  36.0 NA    NA   NA     1     1 70  50
#> 23  3  72.0 NA    NA   NA     1     1 70  50
#> 24  3 120.0 NA    NA   NA     1     1 70  50
#> 25  4   0.5 NA    NA   NA     2     1 60  45
#> 26  4   1.0 NA    NA   NA     2     1 60  45
#> 27  4   2.0 NA    NA   NA     2     1 60  45
#> 28  4   6.0 NA    NA   NA     2     1 60  45
#> 29  4  24.0 NA    NA   NA     2     1 60  45
#> 30  4  36.0 NA    NA   NA     2     1 60  45
#> 31  4  72.0 NA    NA   NA     2     1 60  45
#> 32  4 120.0 NA    NA   NA     2     1 60  45
#> 33  5   0.5 NA    NA   NA     2     1 60  45
#> 34  5   1.0 NA    NA   NA     2     1 60  45
#> 35  5   2.0 NA    NA   NA     2     1 60  45
#> 36  5   6.0 NA    NA   NA     2     1 60  45
#> 37  5  24.0 NA    NA   NA     2     1 60  45
#> 38  5  36.0 NA    NA   NA     2     1 60  45
#> 39  5  72.0 NA    NA   NA     2     1 60  45
#> 40  5 120.0 NA    NA   NA     2     1 60  45
#> 41  6   0.5 NA    NA   NA     2     1 60  45
#> 42  6   1.0 NA    NA   NA     2     1 60  45
#> 43  6   2.0 NA    NA   NA     2     1 60  45
#> 44  6   6.0 NA    NA   NA     2     1 60  45
#> 45  6  24.0 NA    NA   NA     2     1 60  45
#> 46  6  36.0 NA    NA   NA     2     1 60  45
#> 47  6  72.0 NA    NA   NA     2     1 60  45
#> 48  6 120.0 NA    NA   NA     2     1 60  45

# generate random deviations in WT for each individual
df_10 <- model_prediction(design=design_3,DV=TRUE,
                          manipulation=expression({for(id in unique(ID)) 
                              WT[ID==id] = rnorm(1,WT[ID==id],WT[ID==id]*0.1);id <- NULL}))
head(df_10,n=20)
#>    ID  Time DV IPRED PRED Group Model       WT AGE
#> 1   1   0.5 NA    NA   NA     1     1 69.30103  50
#> 2   1   1.0 NA    NA   NA     1     1 69.30103  50
#> 3   1   2.0 NA    NA   NA     1     1 69.30103  50
#> 4   1   6.0 NA    NA   NA     1     1 69.30103  50
#> 5   1  24.0 NA    NA   NA     1     1 69.30103  50
#> 6   1  36.0 NA    NA   NA     1     1 69.30103  50
#> 7   1  72.0 NA    NA   NA     1     1 69.30103  50
#> 8   1 120.0 NA    NA   NA     1     1 69.30103  50
#> 9   2   0.5 NA    NA   NA     1     1 70.32465  50
#> 10  2   1.0 NA    NA   NA     1     1 70.32465  50
#> 11  2   2.0 NA    NA   NA     1     1 70.32465  50
#> 12  2   6.0 NA    NA   NA     1     1 70.32465  50
#> 13  2  24.0 NA    NA   NA     1     1 70.32465  50
#> 14  2  36.0 NA    NA   NA     1     1 70.32465  50
#> 15  2  72.0 NA    NA   NA     1     1 70.32465  50
#> 16  2 120.0 NA    NA   NA     1     1 70.32465  50
#> 17  3   0.5 NA    NA   NA     1     1 77.19992  50
#> 18  3   1.0 NA    NA   NA     1     1 77.19992  50
#> 19  3   2.0 NA    NA   NA     1     1 77.19992  50
#> 20  3   6.0 NA    NA   NA     1     1 77.19992  50

# generate random deviations in WT and AGE for each individual
df_11 <- model_prediction(design=design_3,DV=TRUE,
                          manipulation=list(
                            expression(for(id in unique(ID)) 
                              WT[ID==id] = rnorm(1,WT[ID==id],WT[ID==id]*0.1)),
                            expression(for(id in unique(ID)) 
                              AGE[ID==id] = rnorm(1,AGE[ID==id],AGE[ID==id]*0.2)),
                            expression(id <- NULL)
                          ))
head(df_10,n=20)
#>    ID  Time DV IPRED PRED Group Model       WT AGE
#> 1   1   0.5 NA    NA   NA     1     1 69.30103  50
#> 2   1   1.0 NA    NA   NA     1     1 69.30103  50
#> 3   1   2.0 NA    NA   NA     1     1 69.30103  50
#> 4   1   6.0 NA    NA   NA     1     1 69.30103  50
#> 5   1  24.0 NA    NA   NA     1     1 69.30103  50
#> 6   1  36.0 NA    NA   NA     1     1 69.30103  50
#> 7   1  72.0 NA    NA   NA     1     1 69.30103  50
#> 8   1 120.0 NA    NA   NA     1     1 69.30103  50
#> 9   2   0.5 NA    NA   NA     1     1 70.32465  50
#> 10  2   1.0 NA    NA   NA     1     1 70.32465  50
#> 11  2   2.0 NA    NA   NA     1     1 70.32465  50
#> 12  2   6.0 NA    NA   NA     1     1 70.32465  50
#> 13  2  24.0 NA    NA   NA     1     1 70.32465  50
#> 14  2  36.0 NA    NA   NA     1     1 70.32465  50
#> 15  2  72.0 NA    NA   NA     1     1 70.32465  50
#> 16  2 120.0 NA    NA   NA     1     1 70.32465  50
#> 17  3   0.5 NA    NA   NA     1     1 77.19992  50
#> 18  3   1.0 NA    NA   NA     1     1 77.19992  50
#> 19  3   2.0 NA    NA   NA     1     1 77.19992  50
#> 20  3   6.0 NA    NA   NA     1     1 77.19992  50

## create dosing rows 
dosing_1 <- list(list(AMT=1000,RATE=NA,Time=0.5),list(AMT=3000,RATE=NA,Time=0.5))
dosing_2 <- list(list(AMT=1000,RATE=NA,Time=0.5))
dosing_3 <- list(list(AMT=1000,Time=0.5))
dosing_4 <- list(list(AMT=c(1000,20),Time=c(0.5,10))) # multiple dosing


(df_12 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_1))
#>    ID  Time DV IPRED PRED  AMT RATE Group Model WT AGE
#> 1   1   0.5 NA    NA   NA 1000   NA     1     1 70  50
#> 2   1   0.5 NA    NA   NA   NA   NA     1     1 70  50
#> 3   1   1.0 NA    NA   NA   NA   NA     1     1 70  50
#> 4   1   2.0 NA    NA   NA   NA   NA     1     1 70  50
#> 5   1   6.0 NA    NA   NA   NA   NA     1     1 70  50
#> 6   1  24.0 NA    NA   NA   NA   NA     1     1 70  50
#> 7   1  36.0 NA    NA   NA   NA   NA     1     1 70  50
#> 8   1  72.0 NA    NA   NA   NA   NA     1     1 70  50
#> 9   1 120.0 NA    NA   NA   NA   NA     1     1 70  50
#> 10  2   0.5 NA    NA   NA 1000   NA     1     1 70  50
#> 11  2   0.5 NA    NA   NA   NA   NA     1     1 70  50
#> 12  2   1.0 NA    NA   NA   NA   NA     1     1 70  50
#> 13  2   2.0 NA    NA   NA   NA   NA     1     1 70  50
#> 14  2   6.0 NA    NA   NA   NA   NA     1     1 70  50
#> 15  2  24.0 NA    NA   NA   NA   NA     1     1 70  50
#> 16  2  36.0 NA    NA   NA   NA   NA     1     1 70  50
#> 17  2  72.0 NA    NA   NA   NA   NA     1     1 70  50
#> 18  2 120.0 NA    NA   NA   NA   NA     1     1 70  50
#> 19  3   0.5 NA    NA   NA 1000   NA     1     1 70  50
#> 20  3   0.5 NA    NA   NA   NA   NA     1     1 70  50
#> 21  3   1.0 NA    NA   NA   NA   NA     1     1 70  50
#> 22  3   2.0 NA    NA   NA   NA   NA     1     1 70  50
#> 23  3   6.0 NA    NA   NA   NA   NA     1     1 70  50
#> 24  3  24.0 NA    NA   NA   NA   NA     1     1 70  50
#> 25  3  36.0 NA    NA   NA   NA   NA     1     1 70  50
#> 26  3  72.0 NA    NA   NA   NA   NA     1     1 70  50
#> 27  3 120.0 NA    NA   NA   NA   NA     1     1 70  50
#> 28  4   0.5 NA    NA   NA 3000   NA     2     1 60  45
#> 29  4   0.5 NA    NA   NA   NA   NA     2     1 60  45
#> 30  4   1.0 NA    NA   NA   NA   NA     2     1 60  45
#> 31  4   2.0 NA    NA   NA   NA   NA     2     1 60  45
#> 32  4   6.0 NA    NA   NA   NA   NA     2     1 60  45
#> 33  4  24.0 NA    NA   NA   NA   NA     2     1 60  45
#> 34  4  36.0 NA    NA   NA   NA   NA     2     1 60  45
#> 35  4  72.0 NA    NA   NA   NA   NA     2     1 60  45
#> 36  4 120.0 NA    NA   NA   NA   NA     2     1 60  45
#> 37  5   0.5 NA    NA   NA 3000   NA     2     1 60  45
#> 38  5   0.5 NA    NA   NA   NA   NA     2     1 60  45
#> 39  5   1.0 NA    NA   NA   NA   NA     2     1 60  45
#> 40  5   2.0 NA    NA   NA   NA   NA     2     1 60  45
#> 41  5   6.0 NA    NA   NA   NA   NA     2     1 60  45
#> 42  5  24.0 NA    NA   NA   NA   NA     2     1 60  45
#> 43  5  36.0 NA    NA   NA   NA   NA     2     1 60  45
#> 44  5  72.0 NA    NA   NA   NA   NA     2     1 60  45
#> 45  5 120.0 NA    NA   NA   NA   NA     2     1 60  45
#> 46  6   0.5 NA    NA   NA 3000   NA     2     1 60  45
#> 47  6   0.5 NA    NA   NA   NA   NA     2     1 60  45
#> 48  6   1.0 NA    NA   NA   NA   NA     2     1 60  45
#> 49  6   2.0 NA    NA   NA   NA   NA     2     1 60  45
#> 50  6   6.0 NA    NA   NA   NA   NA     2     1 60  45
#> 51  6  24.0 NA    NA   NA   NA   NA     2     1 60  45
#> 52  6  36.0 NA    NA   NA   NA   NA     2     1 60  45
#> 53  6  72.0 NA    NA   NA   NA   NA     2     1 60  45
#> 54  6 120.0 NA    NA   NA   NA   NA     2     1 60  45
(df_13 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_2))
#>    ID  Time DV IPRED PRED  AMT RATE Group Model WT AGE
#> 1   1   0.5 NA    NA   NA 1000   NA     1     1 70  50
#> 2   1   0.5 NA    NA   NA   NA   NA     1     1 70  50
#> 3   1   1.0 NA    NA   NA   NA   NA     1     1 70  50
#> 4   1   2.0 NA    NA   NA   NA   NA     1     1 70  50
#> 5   1   6.0 NA    NA   NA   NA   NA     1     1 70  50
#> 6   1  24.0 NA    NA   NA   NA   NA     1     1 70  50
#> 7   1  36.0 NA    NA   NA   NA   NA     1     1 70  50
#> 8   1  72.0 NA    NA   NA   NA   NA     1     1 70  50
#> 9   1 120.0 NA    NA   NA   NA   NA     1     1 70  50
#> 10  2   0.5 NA    NA   NA 1000   NA     1     1 70  50
#> 11  2   0.5 NA    NA   NA   NA   NA     1     1 70  50
#> 12  2   1.0 NA    NA   NA   NA   NA     1     1 70  50
#> 13  2   2.0 NA    NA   NA   NA   NA     1     1 70  50
#> 14  2   6.0 NA    NA   NA   NA   NA     1     1 70  50
#> 15  2  24.0 NA    NA   NA   NA   NA     1     1 70  50
#> 16  2  36.0 NA    NA   NA   NA   NA     1     1 70  50
#> 17  2  72.0 NA    NA   NA   NA   NA     1     1 70  50
#> 18  2 120.0 NA    NA   NA   NA   NA     1     1 70  50
#> 19  3   0.5 NA    NA   NA 1000   NA     1     1 70  50
#> 20  3   0.5 NA    NA   NA   NA   NA     1     1 70  50
#> 21  3   1.0 NA    NA   NA   NA   NA     1     1 70  50
#> 22  3   2.0 NA    NA   NA   NA   NA     1     1 70  50
#> 23  3   6.0 NA    NA   NA   NA   NA     1     1 70  50
#> 24  3  24.0 NA    NA   NA   NA   NA     1     1 70  50
#> 25  3  36.0 NA    NA   NA   NA   NA     1     1 70  50
#> 26  3  72.0 NA    NA   NA   NA   NA     1     1 70  50
#> 27  3 120.0 NA    NA   NA   NA   NA     1     1 70  50
#> 28  4   0.5 NA    NA   NA 1000   NA     2     1 60  45
#> 29  4   0.5 NA    NA   NA   NA   NA     2     1 60  45
#> 30  4   1.0 NA    NA   NA   NA   NA     2     1 60  45
#> 31  4   2.0 NA    NA   NA   NA   NA     2     1 60  45
#> 32  4   6.0 NA    NA   NA   NA   NA     2     1 60  45
#> 33  4  24.0 NA    NA   NA   NA   NA     2     1 60  45
#> 34  4  36.0 NA    NA   NA   NA   NA     2     1 60  45
#> 35  4  72.0 NA    NA   NA   NA   NA     2     1 60  45
#> 36  4 120.0 NA    NA   NA   NA   NA     2     1 60  45
#> 37  5   0.5 NA    NA   NA 1000   NA     2     1 60  45
#> 38  5   0.5 NA    NA   NA   NA   NA     2     1 60  45
#> 39  5   1.0 NA    NA   NA   NA   NA     2     1 60  45
#> 40  5   2.0 NA    NA   NA   NA   NA     2     1 60  45
#> 41  5   6.0 NA    NA   NA   NA   NA     2     1 60  45
#> 42  5  24.0 NA    NA   NA   NA   NA     2     1 60  45
#> 43  5  36.0 NA    NA   NA   NA   NA     2     1 60  45
#> 44  5  72.0 NA    NA   NA   NA   NA     2     1 60  45
#> 45  5 120.0 NA    NA   NA   NA   NA     2     1 60  45
#> 46  6   0.5 NA    NA   NA 1000   NA     2     1 60  45
#> 47  6   0.5 NA    NA   NA   NA   NA     2     1 60  45
#> 48  6   1.0 NA    NA   NA   NA   NA     2     1 60  45
#> 49  6   2.0 NA    NA   NA   NA   NA     2     1 60  45
#> 50  6   6.0 NA    NA   NA   NA   NA     2     1 60  45
#> 51  6  24.0 NA    NA   NA   NA   NA     2     1 60  45
#> 52  6  36.0 NA    NA   NA   NA   NA     2     1 60  45
#> 53  6  72.0 NA    NA   NA   NA   NA     2     1 60  45
#> 54  6 120.0 NA    NA   NA   NA   NA     2     1 60  45
(df_14 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_3))
#>    ID  Time DV IPRED PRED  AMT Group Model WT AGE
#> 1   1   0.5 NA    NA   NA 1000     1     1 70  50
#> 2   1   0.5 NA    NA   NA   NA     1     1 70  50
#> 3   1   1.0 NA    NA   NA   NA     1     1 70  50
#> 4   1   2.0 NA    NA   NA   NA     1     1 70  50
#> 5   1   6.0 NA    NA   NA   NA     1     1 70  50
#> 6   1  24.0 NA    NA   NA   NA     1     1 70  50
#> 7   1  36.0 NA    NA   NA   NA     1     1 70  50
#> 8   1  72.0 NA    NA   NA   NA     1     1 70  50
#> 9   1 120.0 NA    NA   NA   NA     1     1 70  50
#> 10  2   0.5 NA    NA   NA 1000     1     1 70  50
#> 11  2   0.5 NA    NA   NA   NA     1     1 70  50
#> 12  2   1.0 NA    NA   NA   NA     1     1 70  50
#> 13  2   2.0 NA    NA   NA   NA     1     1 70  50
#> 14  2   6.0 NA    NA   NA   NA     1     1 70  50
#> 15  2  24.0 NA    NA   NA   NA     1     1 70  50
#> 16  2  36.0 NA    NA   NA   NA     1     1 70  50
#> 17  2  72.0 NA    NA   NA   NA     1     1 70  50
#> 18  2 120.0 NA    NA   NA   NA     1     1 70  50
#> 19  3   0.5 NA    NA   NA 1000     1     1 70  50
#> 20  3   0.5 NA    NA   NA   NA     1     1 70  50
#> 21  3   1.0 NA    NA   NA   NA     1     1 70  50
#> 22  3   2.0 NA    NA   NA   NA     1     1 70  50
#> 23  3   6.0 NA    NA   NA   NA     1     1 70  50
#> 24  3  24.0 NA    NA   NA   NA     1     1 70  50
#> 25  3  36.0 NA    NA   NA   NA     1     1 70  50
#> 26  3  72.0 NA    NA   NA   NA     1     1 70  50
#> 27  3 120.0 NA    NA   NA   NA     1     1 70  50
#> 28  4   0.5 NA    NA   NA 1000     2     1 60  45
#> 29  4   0.5 NA    NA   NA   NA     2     1 60  45
#> 30  4   1.0 NA    NA   NA   NA     2     1 60  45
#> 31  4   2.0 NA    NA   NA   NA     2     1 60  45
#> 32  4   6.0 NA    NA   NA   NA     2     1 60  45
#> 33  4  24.0 NA    NA   NA   NA     2     1 60  45
#> 34  4  36.0 NA    NA   NA   NA     2     1 60  45
#> 35  4  72.0 NA    NA   NA   NA     2     1 60  45
#> 36  4 120.0 NA    NA   NA   NA     2     1 60  45
#> 37  5   0.5 NA    NA   NA 1000     2     1 60  45
#> 38  5   0.5 NA    NA   NA   NA     2     1 60  45
#> 39  5   1.0 NA    NA   NA   NA     2     1 60  45
#> 40  5   2.0 NA    NA   NA   NA     2     1 60  45
#> 41  5   6.0 NA    NA   NA   NA     2     1 60  45
#> 42  5  24.0 NA    NA   NA   NA     2     1 60  45
#> 43  5  36.0 NA    NA   NA   NA     2     1 60  45
#> 44  5  72.0 NA    NA   NA   NA     2     1 60  45
#> 45  5 120.0 NA    NA   NA   NA     2     1 60  45
#> 46  6   0.5 NA    NA   NA 1000     2     1 60  45
#> 47  6   0.5 NA    NA   NA   NA     2     1 60  45
#> 48  6   1.0 NA    NA   NA   NA     2     1 60  45
#> 49  6   2.0 NA    NA   NA   NA     2     1 60  45
#> 50  6   6.0 NA    NA   NA   NA     2     1 60  45
#> 51  6  24.0 NA    NA   NA   NA     2     1 60  45
#> 52  6  36.0 NA    NA   NA   NA     2     1 60  45
#> 53  6  72.0 NA    NA   NA   NA     2     1 60  45
#> 54  6 120.0 NA    NA   NA   NA     2     1 60  45
(df_15 <- model_prediction(design=design_3,DV=TRUE,dosing=dosing_4))
#>    ID  Time DV IPRED PRED  AMT Group Model WT AGE
#> 1   1   0.5 NA    NA   NA 1000     1     1 70  50
#> 2   1   0.5 NA    NA   NA   NA     1     1 70  50
#> 3   1   1.0 NA    NA   NA   NA     1     1 70  50
#> 4   1   2.0 NA    NA   NA   NA     1     1 70  50
#> 5   1   6.0 NA    NA   NA   NA     1     1 70  50
#> 6   1  10.0 NA    NA   NA   20     1     1 70  50
#> 7   1  24.0 NA    NA   NA   NA     1     1 70  50
#> 8   1  36.0 NA    NA   NA   NA     1     1 70  50
#> 9   1  72.0 NA    NA   NA   NA     1     1 70  50
#> 10  1 120.0 NA    NA   NA   NA     1     1 70  50
#> 11  2   0.5 NA    NA   NA 1000     1     1 70  50
#> 12  2   0.5 NA    NA   NA   NA     1     1 70  50
#> 13  2   1.0 NA    NA   NA   NA     1     1 70  50
#> 14  2   2.0 NA    NA   NA   NA     1     1 70  50
#> 15  2   6.0 NA    NA   NA   NA     1     1 70  50
#> 16  2  10.0 NA    NA   NA   20     1     1 70  50
#> 17  2  24.0 NA    NA   NA   NA     1     1 70  50
#> 18  2  36.0 NA    NA   NA   NA     1     1 70  50
#> 19  2  72.0 NA    NA   NA   NA     1     1 70  50
#> 20  2 120.0 NA    NA   NA   NA     1     1 70  50
#> 21  3   0.5 NA    NA   NA 1000     1     1 70  50
#> 22  3   0.5 NA    NA   NA   NA     1     1 70  50
#> 23  3   1.0 NA    NA   NA   NA     1     1 70  50
#> 24  3   2.0 NA    NA   NA   NA     1     1 70  50
#> 25  3   6.0 NA    NA   NA   NA     1     1 70  50
#> 26  3  10.0 NA    NA   NA   20     1     1 70  50
#> 27  3  24.0 NA    NA   NA   NA     1     1 70  50
#> 28  3  36.0 NA    NA   NA   NA     1     1 70  50
#> 29  3  72.0 NA    NA   NA   NA     1     1 70  50
#> 30  3 120.0 NA    NA   NA   NA     1     1 70  50
#> 31  4   0.5 NA    NA   NA 1000     2     1 60  45
#> 32  4   0.5 NA    NA   NA   NA     2     1 60  45
#> 33  4   1.0 NA    NA   NA   NA     2     1 60  45
#> 34  4   2.0 NA    NA   NA   NA     2     1 60  45
#> 35  4   6.0 NA    NA   NA   NA     2     1 60  45
#> 36  4  10.0 NA    NA   NA   20     2     1 60  45
#> 37  4  24.0 NA    NA   NA   NA     2     1 60  45
#> 38  4  36.0 NA    NA   NA   NA     2     1 60  45
#> 39  4  72.0 NA    NA   NA   NA     2     1 60  45
#> 40  4 120.0 NA    NA   NA   NA     2     1 60  45
#> 41  5   0.5 NA    NA   NA 1000     2     1 60  45
#> 42  5   0.5 NA    NA   NA   NA     2     1 60  45
#> 43  5   1.0 NA    NA   NA   NA     2     1 60  45
#> 44  5   2.0 NA    NA   NA   NA     2     1 60  45
#> 45  5   6.0 NA    NA   NA   NA     2     1 60  45
#> 46  5  10.0 NA    NA   NA   20     2     1 60  45
#> 47  5  24.0 NA    NA   NA   NA     2     1 60  45
#> 48  5  36.0 NA    NA   NA   NA     2     1 60  45
#> 49  5  72.0 NA    NA   NA   NA     2     1 60  45
#> 50  5 120.0 NA    NA   NA   NA     2     1 60  45
#> 51  6   0.5 NA    NA   NA 1000     2     1 60  45
#> 52  6   0.5 NA    NA   NA   NA     2     1 60  45
#> 53  6   1.0 NA    NA   NA   NA     2     1 60  45
#> 54  6   2.0 NA    NA   NA   NA     2     1 60  45
#> 55  6   6.0 NA    NA   NA   NA     2     1 60  45
#> 56  6  10.0 NA    NA   NA   20     2     1 60  45
#> 57  6  24.0 NA    NA   NA   NA     2     1 60  45
#> 58  6  36.0 NA    NA   NA   NA     2     1 60  45
#> 59  6  72.0 NA    NA   NA   NA     2     1 60  45
#> 60  6 120.0 NA    NA   NA   NA     2     1 60  45


model_prediction(design=design_3,DV=TRUE,dosing=dosing_4,model_num_points = 10)
#>     ID      Time DV IPRED PRED  AMT Group Model WT AGE
#> 1    1   0.50000 NA    NA   NA 1000     1     1 70  50
#> 2    1   0.50000 NA    NA   NA   NA     1     1 70  50
#> 3    1   1.00000 NA    NA   NA   NA     1     1 70  50
#> 4    1   2.00000 NA    NA   NA   NA     1     1 70  50
#> 5    1   6.00000 NA    NA   NA   NA     1     1 70  50
#> 6    1  10.00000 NA    NA   NA   20     1     1 70  50
#> 7    1  13.77778 NA    NA   NA   NA     1     1 70  50
#> 8    1  24.00000 NA    NA   NA   NA     1     1 70  50
#> 9    1  27.05556 NA    NA   NA   NA     1     1 70  50
#> 10   1  36.00000 NA    NA   NA   NA     1     1 70  50
#> 11   1  40.33333 NA    NA   NA   NA     1     1 70  50
#> 12   1  53.61111 NA    NA   NA   NA     1     1 70  50
#> 13   1  66.88889 NA    NA   NA   NA     1     1 70  50
#> 14   1  72.00000 NA    NA   NA   NA     1     1 70  50
#> 15   1  80.16667 NA    NA   NA   NA     1     1 70  50
#> 16   1  93.44444 NA    NA   NA   NA     1     1 70  50
#> 17   1 106.72222 NA    NA   NA   NA     1     1 70  50
#> 18   1 120.00000 NA    NA   NA   NA     1     1 70  50
#> 19   2   0.50000 NA    NA   NA 1000     1     1 70  50
#> 20   2   0.50000 NA    NA   NA   NA     1     1 70  50
#> 21   2   1.00000 NA    NA   NA   NA     1     1 70  50
#> 22   2   2.00000 NA    NA   NA   NA     1     1 70  50
#> 23   2   6.00000 NA    NA   NA   NA     1     1 70  50
#> 24   2  10.00000 NA    NA   NA   20     1     1 70  50
#> 25   2  13.77778 NA    NA   NA   NA     1     1 70  50
#> 26   2  24.00000 NA    NA   NA   NA     1     1 70  50
#> 27   2  27.05556 NA    NA   NA   NA     1     1 70  50
#> 28   2  36.00000 NA    NA   NA   NA     1     1 70  50
#> 29   2  40.33333 NA    NA   NA   NA     1     1 70  50
#> 30   2  53.61111 NA    NA   NA   NA     1     1 70  50
#> 31   2  66.88889 NA    NA   NA   NA     1     1 70  50
#> 32   2  72.00000 NA    NA   NA   NA     1     1 70  50
#> 33   2  80.16667 NA    NA   NA   NA     1     1 70  50
#> 34   2  93.44444 NA    NA   NA   NA     1     1 70  50
#> 35   2 106.72222 NA    NA   NA   NA     1     1 70  50
#> 36   2 120.00000 NA    NA   NA   NA     1     1 70  50
#> 37   3   0.50000 NA    NA   NA 1000     1     1 70  50
#> 38   3   0.50000 NA    NA   NA   NA     1     1 70  50
#> 39   3   1.00000 NA    NA   NA   NA     1     1 70  50
#> 40   3   2.00000 NA    NA   NA   NA     1     1 70  50
#> 41   3   6.00000 NA    NA   NA   NA     1     1 70  50
#> 42   3  10.00000 NA    NA   NA   20     1     1 70  50
#> 43   3  13.77778 NA    NA   NA   NA     1     1 70  50
#> 44   3  24.00000 NA    NA   NA   NA     1     1 70  50
#> 45   3  27.05556 NA    NA   NA   NA     1     1 70  50
#> 46   3  36.00000 NA    NA   NA   NA     1     1 70  50
#> 47   3  40.33333 NA    NA   NA   NA     1     1 70  50
#> 48   3  53.61111 NA    NA   NA   NA     1     1 70  50
#> 49   3  66.88889 NA    NA   NA   NA     1     1 70  50
#> 50   3  72.00000 NA    NA   NA   NA     1     1 70  50
#> 51   3  80.16667 NA    NA   NA   NA     1     1 70  50
#> 52   3  93.44444 NA    NA   NA   NA     1     1 70  50
#> 53   3 106.72222 NA    NA   NA   NA     1     1 70  50
#> 54   3 120.00000 NA    NA   NA   NA     1     1 70  50
#> 55   4   0.50000 NA    NA   NA 1000     2     1 60  45
#> 56   4   0.50000 NA    NA   NA   NA     2     1 60  45
#> 57   4   1.00000 NA    NA   NA   NA     2     1 60  45
#> 58   4   2.00000 NA    NA   NA   NA     2     1 60  45
#> 59   4   6.00000 NA    NA   NA   NA     2     1 60  45
#> 60   4  10.00000 NA    NA   NA   20     2     1 60  45
#> 61   4  13.77778 NA    NA   NA   NA     2     1 60  45
#> 62   4  24.00000 NA    NA   NA   NA     2     1 60  45
#> 63   4  27.05556 NA    NA   NA   NA     2     1 60  45
#> 64   4  36.00000 NA    NA   NA   NA     2     1 60  45
#> 65   4  40.33333 NA    NA   NA   NA     2     1 60  45
#> 66   4  53.61111 NA    NA   NA   NA     2     1 60  45
#> 67   4  66.88889 NA    NA   NA   NA     2     1 60  45
#> 68   4  72.00000 NA    NA   NA   NA     2     1 60  45
#> 69   4  80.16667 NA    NA   NA   NA     2     1 60  45
#> 70   4  93.44444 NA    NA   NA   NA     2     1 60  45
#> 71   4 106.72222 NA    NA   NA   NA     2     1 60  45
#> 72   4 120.00000 NA    NA   NA   NA     2     1 60  45
#> 73   5   0.50000 NA    NA   NA 1000     2     1 60  45
#> 74   5   0.50000 NA    NA   NA   NA     2     1 60  45
#> 75   5   1.00000 NA    NA   NA   NA     2     1 60  45
#> 76   5   2.00000 NA    NA   NA   NA     2     1 60  45
#> 77   5   6.00000 NA    NA   NA   NA     2     1 60  45
#> 78   5  10.00000 NA    NA   NA   20     2     1 60  45
#> 79   5  13.77778 NA    NA   NA   NA     2     1 60  45
#> 80   5  24.00000 NA    NA   NA   NA     2     1 60  45
#> 81   5  27.05556 NA    NA   NA   NA     2     1 60  45
#> 82   5  36.00000 NA    NA   NA   NA     2     1 60  45
#> 83   5  40.33333 NA    NA   NA   NA     2     1 60  45
#> 84   5  53.61111 NA    NA   NA   NA     2     1 60  45
#> 85   5  66.88889 NA    NA   NA   NA     2     1 60  45
#> 86   5  72.00000 NA    NA   NA   NA     2     1 60  45
#> 87   5  80.16667 NA    NA   NA   NA     2     1 60  45
#> 88   5  93.44444 NA    NA   NA   NA     2     1 60  45
#> 89   5 106.72222 NA    NA   NA   NA     2     1 60  45
#> 90   5 120.00000 NA    NA   NA   NA     2     1 60  45
#> 91   6   0.50000 NA    NA   NA 1000     2     1 60  45
#> 92   6   0.50000 NA    NA   NA   NA     2     1 60  45
#> 93   6   1.00000 NA    NA   NA   NA     2     1 60  45
#> 94   6   2.00000 NA    NA   NA   NA     2     1 60  45
#> 95   6   6.00000 NA    NA   NA   NA     2     1 60  45
#> 96   6  10.00000 NA    NA   NA   20     2     1 60  45
#> 97   6  13.77778 NA    NA   NA   NA     2     1 60  45
#> 98   6  24.00000 NA    NA   NA   NA     2     1 60  45
#> 99   6  27.05556 NA    NA   NA   NA     2     1 60  45
#> 100  6  36.00000 NA    NA   NA   NA     2     1 60  45
#> 101  6  40.33333 NA    NA   NA   NA     2     1 60  45
#> 102  6  53.61111 NA    NA   NA   NA     2     1 60  45
#> 103  6  66.88889 NA    NA   NA   NA     2     1 60  45
#> 104  6  72.00000 NA    NA   NA   NA     2     1 60  45
#> 105  6  80.16667 NA    NA   NA   NA     2     1 60  45
#> 106  6  93.44444 NA    NA   NA   NA     2     1 60  45
#> 107  6 106.72222 NA    NA   NA   NA     2     1 60  45
#> 108  6 120.00000 NA    NA   NA   NA     2     1 60  45
model_prediction(design=design_3,DV=TRUE,dosing=dosing_4,model_num_points = 10,model_minxt=20)
#>     ID      Time DV IPRED PRED  AMT Group Model WT AGE
#> 1    1   0.50000 NA    NA   NA 1000     1     1 70  50
#> 2    1   0.50000 NA    NA   NA   NA     1     1 70  50
#> 3    1   1.00000 NA    NA   NA   NA     1     1 70  50
#> 4    1   2.00000 NA    NA   NA   NA     1     1 70  50
#> 5    1   6.00000 NA    NA   NA   NA     1     1 70  50
#> 6    1  10.00000 NA    NA   NA   20     1     1 70  50
#> 7    1  20.00000 NA    NA   NA   NA     1     1 70  50
#> 8    1  24.00000 NA    NA   NA   NA     1     1 70  50
#> 9    1  31.11111 NA    NA   NA   NA     1     1 70  50
#> 10   1  36.00000 NA    NA   NA   NA     1     1 70  50
#> 11   1  42.22222 NA    NA   NA   NA     1     1 70  50
#> 12   1  53.33333 NA    NA   NA   NA     1     1 70  50
#> 13   1  64.44444 NA    NA   NA   NA     1     1 70  50
#> 14   1  72.00000 NA    NA   NA   NA     1     1 70  50
#> 15   1  75.55556 NA    NA   NA   NA     1     1 70  50
#> 16   1  86.66667 NA    NA   NA   NA     1     1 70  50
#> 17   1  97.77778 NA    NA   NA   NA     1     1 70  50
#> 18   1 108.88889 NA    NA   NA   NA     1     1 70  50
#> 19   1 120.00000 NA    NA   NA   NA     1     1 70  50
#> 20   2   0.50000 NA    NA   NA 1000     1     1 70  50
#> 21   2   0.50000 NA    NA   NA   NA     1     1 70  50
#> 22   2   1.00000 NA    NA   NA   NA     1     1 70  50
#> 23   2   2.00000 NA    NA   NA   NA     1     1 70  50
#> 24   2   6.00000 NA    NA   NA   NA     1     1 70  50
#> 25   2  10.00000 NA    NA   NA   20     1     1 70  50
#> 26   2  20.00000 NA    NA   NA   NA     1     1 70  50
#> 27   2  24.00000 NA    NA   NA   NA     1     1 70  50
#> 28   2  31.11111 NA    NA   NA   NA     1     1 70  50
#> 29   2  36.00000 NA    NA   NA   NA     1     1 70  50
#> 30   2  42.22222 NA    NA   NA   NA     1     1 70  50
#> 31   2  53.33333 NA    NA   NA   NA     1     1 70  50
#> 32   2  64.44444 NA    NA   NA   NA     1     1 70  50
#> 33   2  72.00000 NA    NA   NA   NA     1     1 70  50
#> 34   2  75.55556 NA    NA   NA   NA     1     1 70  50
#> 35   2  86.66667 NA    NA   NA   NA     1     1 70  50
#> 36   2  97.77778 NA    NA   NA   NA     1     1 70  50
#> 37   2 108.88889 NA    NA   NA   NA     1     1 70  50
#> 38   2 120.00000 NA    NA   NA   NA     1     1 70  50
#> 39   3   0.50000 NA    NA   NA 1000     1     1 70  50
#> 40   3   0.50000 NA    NA   NA   NA     1     1 70  50
#> 41   3   1.00000 NA    NA   NA   NA     1     1 70  50
#> 42   3   2.00000 NA    NA   NA   NA     1     1 70  50
#> 43   3   6.00000 NA    NA   NA   NA     1     1 70  50
#> 44   3  10.00000 NA    NA   NA   20     1     1 70  50
#> 45   3  20.00000 NA    NA   NA   NA     1     1 70  50
#> 46   3  24.00000 NA    NA   NA   NA     1     1 70  50
#> 47   3  31.11111 NA    NA   NA   NA     1     1 70  50
#> 48   3  36.00000 NA    NA   NA   NA     1     1 70  50
#> 49   3  42.22222 NA    NA   NA   NA     1     1 70  50
#> 50   3  53.33333 NA    NA   NA   NA     1     1 70  50
#> 51   3  64.44444 NA    NA   NA   NA     1     1 70  50
#> 52   3  72.00000 NA    NA   NA   NA     1     1 70  50
#> 53   3  75.55556 NA    NA   NA   NA     1     1 70  50
#> 54   3  86.66667 NA    NA   NA   NA     1     1 70  50
#> 55   3  97.77778 NA    NA   NA   NA     1     1 70  50
#> 56   3 108.88889 NA    NA   NA   NA     1     1 70  50
#> 57   3 120.00000 NA    NA   NA   NA     1     1 70  50
#> 58   4   0.50000 NA    NA   NA 1000     2     1 60  45
#> 59   4   0.50000 NA    NA   NA   NA     2     1 60  45
#> 60   4   1.00000 NA    NA   NA   NA     2     1 60  45
#> 61   4   2.00000 NA    NA   NA   NA     2     1 60  45
#> 62   4   6.00000 NA    NA   NA   NA     2     1 60  45
#> 63   4  10.00000 NA    NA   NA   20     2     1 60  45
#> 64   4  20.00000 NA    NA   NA   NA     2     1 60  45
#> 65   4  24.00000 NA    NA   NA   NA     2     1 60  45
#> 66   4  31.11111 NA    NA   NA   NA     2     1 60  45
#> 67   4  36.00000 NA    NA   NA   NA     2     1 60  45
#> 68   4  42.22222 NA    NA   NA   NA     2     1 60  45
#> 69   4  53.33333 NA    NA   NA   NA     2     1 60  45
#> 70   4  64.44444 NA    NA   NA   NA     2     1 60  45
#> 71   4  72.00000 NA    NA   NA   NA     2     1 60  45
#> 72   4  75.55556 NA    NA   NA   NA     2     1 60  45
#> 73   4  86.66667 NA    NA   NA   NA     2     1 60  45
#> 74   4  97.77778 NA    NA   NA   NA     2     1 60  45
#> 75   4 108.88889 NA    NA   NA   NA     2     1 60  45
#> 76   4 120.00000 NA    NA   NA   NA     2     1 60  45
#> 77   5   0.50000 NA    NA   NA 1000     2     1 60  45
#> 78   5   0.50000 NA    NA   NA   NA     2     1 60  45
#> 79   5   1.00000 NA    NA   NA   NA     2     1 60  45
#> 80   5   2.00000 NA    NA   NA   NA     2     1 60  45
#> 81   5   6.00000 NA    NA   NA   NA     2     1 60  45
#> 82   5  10.00000 NA    NA   NA   20     2     1 60  45
#> 83   5  20.00000 NA    NA   NA   NA     2     1 60  45
#> 84   5  24.00000 NA    NA   NA   NA     2     1 60  45
#> 85   5  31.11111 NA    NA   NA   NA     2     1 60  45
#> 86   5  36.00000 NA    NA   NA   NA     2     1 60  45
#> 87   5  42.22222 NA    NA   NA   NA     2     1 60  45
#> 88   5  53.33333 NA    NA   NA   NA     2     1 60  45
#> 89   5  64.44444 NA    NA   NA   NA     2     1 60  45
#> 90   5  72.00000 NA    NA   NA   NA     2     1 60  45
#> 91   5  75.55556 NA    NA   NA   NA     2     1 60  45
#> 92   5  86.66667 NA    NA   NA   NA     2     1 60  45
#> 93   5  97.77778 NA    NA   NA   NA     2     1 60  45
#> 94   5 108.88889 NA    NA   NA   NA     2     1 60  45
#> 95   5 120.00000 NA    NA   NA   NA     2     1 60  45
#> 96   6   0.50000 NA    NA   NA 1000     2     1 60  45
#> 97   6   0.50000 NA    NA   NA   NA     2     1 60  45
#> 98   6   1.00000 NA    NA   NA   NA     2     1 60  45
#> 99   6   2.00000 NA    NA   NA   NA     2     1 60  45
#> 100  6   6.00000 NA    NA   NA   NA     2     1 60  45
#> 101  6  10.00000 NA    NA   NA   20     2     1 60  45
#> 102  6  20.00000 NA    NA   NA   NA     2     1 60  45
#> 103  6  24.00000 NA    NA   NA   NA     2     1 60  45
#> 104  6  31.11111 NA    NA   NA   NA     2     1 60  45
#> 105  6  36.00000 NA    NA   NA   NA     2     1 60  45
#> 106  6  42.22222 NA    NA   NA   NA     2     1 60  45
#> 107  6  53.33333 NA    NA   NA   NA     2     1 60  45
#> 108  6  64.44444 NA    NA   NA   NA     2     1 60  45
#> 109  6  72.00000 NA    NA   NA   NA     2     1 60  45
#> 110  6  75.55556 NA    NA   NA   NA     2     1 60  45
#> 111  6  86.66667 NA    NA   NA   NA     2     1 60  45
#> 112  6  97.77778 NA    NA   NA   NA     2     1 60  45
#> 113  6 108.88889 NA    NA   NA   NA     2     1 60  45
#> 114  6 120.00000 NA    NA   NA   NA     2     1 60  45

design_4 <- list(
  xt=c( 0.5,1,2,6,24,36,72,120),
  model_switch=c(1,1,1,1,2,2,2,2),
  m=2,
  groupsize=3,
  a=list(c(WT=70,AGE=50),c(AGE=45,WT=60)))

model_prediction(design=design_4,DV=TRUE,dosing=dosing_4)
#>    ID  Time DV IPRED PRED  AMT Group Model WT AGE
#> 1   1   0.5 NA    NA   NA 1000     1     1 70  50
#> 2   1   0.5 NA    NA   NA   NA     1     1 70  50
#> 3   1   1.0 NA    NA   NA   NA     1     1 70  50
#> 4   1   2.0 NA    NA   NA   NA     1     1 70  50
#> 5   1   6.0 NA    NA   NA   NA     1     1 70  50
#> 6   1  10.0 NA    NA   NA   20     1     1 70  50
#> 7   1  24.0 NA    NA   NA   NA     1     2 70  50
#> 8   1  36.0 NA    NA   NA   NA     1     2 70  50
#> 9   1  72.0 NA    NA   NA   NA     1     2 70  50
#> 10  1 120.0 NA    NA   NA   NA     1     2 70  50
#> 11  2   0.5 NA    NA   NA 1000     1     1 70  50
#> 12  2   0.5 NA    NA   NA   NA     1     1 70  50
#> 13  2   1.0 NA    NA   NA   NA     1     1 70  50
#> 14  2   2.0 NA    NA   NA   NA     1     1 70  50
#> 15  2   6.0 NA    NA   NA   NA     1     1 70  50
#> 16  2  10.0 NA    NA   NA   20     1     1 70  50
#> 17  2  24.0 NA    NA   NA   NA     1     2 70  50
#> 18  2  36.0 NA    NA   NA   NA     1     2 70  50
#> 19  2  72.0 NA    NA   NA   NA     1     2 70  50
#> 20  2 120.0 NA    NA   NA   NA     1     2 70  50
#> 21  3   0.5 NA    NA   NA 1000     1     1 70  50
#> 22  3   0.5 NA    NA   NA   NA     1     1 70  50
#> 23  3   1.0 NA    NA   NA   NA     1     1 70  50
#> 24  3   2.0 NA    NA   NA   NA     1     1 70  50
#> 25  3   6.0 NA    NA   NA   NA     1     1 70  50
#> 26  3  10.0 NA    NA   NA   20     1     1 70  50
#> 27  3  24.0 NA    NA   NA   NA     1     2 70  50
#> 28  3  36.0 NA    NA   NA   NA     1     2 70  50
#> 29  3  72.0 NA    NA   NA   NA     1     2 70  50
#> 30  3 120.0 NA    NA   NA   NA     1     2 70  50
#> 31  4   0.5 NA    NA   NA 1000     2     1 60  45
#> 32  4   0.5 NA    NA   NA   NA     2     1 60  45
#> 33  4   1.0 NA    NA   NA   NA     2     1 60  45
#> 34  4   2.0 NA    NA   NA   NA     2     1 60  45
#> 35  4   6.0 NA    NA   NA   NA     2     1 60  45
#> 36  4  10.0 NA    NA   NA   20     2     1 60  45
#> 37  4  24.0 NA    NA   NA   NA     2     2 60  45
#> 38  4  36.0 NA    NA   NA   NA     2     2 60  45
#> 39  4  72.0 NA    NA   NA   NA     2     2 60  45
#> 40  4 120.0 NA    NA   NA   NA     2     2 60  45
#> 41  5   0.5 NA    NA   NA 1000     2     1 60  45
#> 42  5   0.5 NA    NA   NA   NA     2     1 60  45
#> 43  5   1.0 NA    NA   NA   NA     2     1 60  45
#> 44  5   2.0 NA    NA   NA   NA     2     1 60  45
#> 45  5   6.0 NA    NA   NA   NA     2     1 60  45
#> 46  5  10.0 NA    NA   NA   20     2     1 60  45
#> 47  5  24.0 NA    NA   NA   NA     2     2 60  45
#> 48  5  36.0 NA    NA   NA   NA     2     2 60  45
#> 49  5  72.0 NA    NA   NA   NA     2     2 60  45
#> 50  5 120.0 NA    NA   NA   NA     2     2 60  45
#> 51  6   0.5 NA    NA   NA 1000     2     1 60  45
#> 52  6   0.5 NA    NA   NA   NA     2     1 60  45
#> 53  6   1.0 NA    NA   NA   NA     2     1 60  45
#> 54  6   2.0 NA    NA   NA   NA     2     1 60  45
#> 55  6   6.0 NA    NA   NA   NA     2     1 60  45
#> 56  6  10.0 NA    NA   NA   20     2     1 60  45
#> 57  6  24.0 NA    NA   NA   NA     2     2 60  45
#> 58  6  36.0 NA    NA   NA   NA     2     2 60  45
#> 59  6  72.0 NA    NA   NA   NA     2     2 60  45
#> 60  6 120.0 NA    NA   NA   NA     2     2 60  45
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = 10)
#>     ID       Time DV IPRED PRED  AMT Group Model WT AGE
#> 1    1   0.500000 NA    NA   NA 1000     1     1 70  50
#> 2    1   0.500000 NA    NA   NA   NA     1     1 70  50
#> 3    1   1.000000 NA    NA   NA   NA     1     1 70  50
#> 4    1   1.111111 NA    NA   NA   NA     1     1 70  50
#> 5    1   1.722222 NA    NA   NA   NA     1     1 70  50
#> 6    1   2.000000 NA    NA   NA   NA     1     1 70  50
#> 7    1   2.333333 NA    NA   NA   NA     1     1 70  50
#> 8    1   2.944444 NA    NA   NA   NA     1     1 70  50
#> 9    1   3.555556 NA    NA   NA   NA     1     1 70  50
#> 10   1   4.166667 NA    NA   NA   NA     1     1 70  50
#> 11   1   4.777778 NA    NA   NA   NA     1     1 70  50
#> 12   1   5.388889 NA    NA   NA   NA     1     1 70  50
#> 13   1   6.000000 NA    NA   NA   NA     1     1 70  50
#> 14   1  10.000000 NA    NA   NA   20     1     1 70  50
#> 15   1  24.000000 NA    NA   NA   NA     1     2 70  50
#> 16   1  34.666667 NA    NA   NA   NA     1     2 70  50
#> 17   1  36.000000 NA    NA   NA   NA     1     2 70  50
#> 18   1  45.333333 NA    NA   NA   NA     1     2 70  50
#> 19   1  56.000000 NA    NA   NA   NA     1     2 70  50
#> 20   1  66.666667 NA    NA   NA   NA     1     2 70  50
#> 21   1  72.000000 NA    NA   NA   NA     1     2 70  50
#> 22   1  77.333333 NA    NA   NA   NA     1     2 70  50
#> 23   1  88.000000 NA    NA   NA   NA     1     2 70  50
#> 24   1  98.666667 NA    NA   NA   NA     1     2 70  50
#> 25   1 109.333333 NA    NA   NA   NA     1     2 70  50
#> 26   1 120.000000 NA    NA   NA   NA     1     2 70  50
#> 27   2   0.500000 NA    NA   NA 1000     1     1 70  50
#> 28   2   0.500000 NA    NA   NA   NA     1     1 70  50
#> 29   2   1.000000 NA    NA   NA   NA     1     1 70  50
#> 30   2   1.111111 NA    NA   NA   NA     1     1 70  50
#> 31   2   1.722222 NA    NA   NA   NA     1     1 70  50
#> 32   2   2.000000 NA    NA   NA   NA     1     1 70  50
#> 33   2   2.333333 NA    NA   NA   NA     1     1 70  50
#> 34   2   2.944444 NA    NA   NA   NA     1     1 70  50
#> 35   2   3.555556 NA    NA   NA   NA     1     1 70  50
#> 36   2   4.166667 NA    NA   NA   NA     1     1 70  50
#> 37   2   4.777778 NA    NA   NA   NA     1     1 70  50
#> 38   2   5.388889 NA    NA   NA   NA     1     1 70  50
#> 39   2   6.000000 NA    NA   NA   NA     1     1 70  50
#> 40   2  10.000000 NA    NA   NA   20     1     1 70  50
#> 41   2  24.000000 NA    NA   NA   NA     1     2 70  50
#> 42   2  34.666667 NA    NA   NA   NA     1     2 70  50
#> 43   2  36.000000 NA    NA   NA   NA     1     2 70  50
#> 44   2  45.333333 NA    NA   NA   NA     1     2 70  50
#> 45   2  56.000000 NA    NA   NA   NA     1     2 70  50
#> 46   2  66.666667 NA    NA   NA   NA     1     2 70  50
#> 47   2  72.000000 NA    NA   NA   NA     1     2 70  50
#> 48   2  77.333333 NA    NA   NA   NA     1     2 70  50
#> 49   2  88.000000 NA    NA   NA   NA     1     2 70  50
#> 50   2  98.666667 NA    NA   NA   NA     1     2 70  50
#> 51   2 109.333333 NA    NA   NA   NA     1     2 70  50
#> 52   2 120.000000 NA    NA   NA   NA     1     2 70  50
#> 53   3   0.500000 NA    NA   NA 1000     1     1 70  50
#> 54   3   0.500000 NA    NA   NA   NA     1     1 70  50
#> 55   3   1.000000 NA    NA   NA   NA     1     1 70  50
#> 56   3   1.111111 NA    NA   NA   NA     1     1 70  50
#> 57   3   1.722222 NA    NA   NA   NA     1     1 70  50
#> 58   3   2.000000 NA    NA   NA   NA     1     1 70  50
#> 59   3   2.333333 NA    NA   NA   NA     1     1 70  50
#> 60   3   2.944444 NA    NA   NA   NA     1     1 70  50
#> 61   3   3.555556 NA    NA   NA   NA     1     1 70  50
#> 62   3   4.166667 NA    NA   NA   NA     1     1 70  50
#> 63   3   4.777778 NA    NA   NA   NA     1     1 70  50
#> 64   3   5.388889 NA    NA   NA   NA     1     1 70  50
#> 65   3   6.000000 NA    NA   NA   NA     1     1 70  50
#> 66   3  10.000000 NA    NA   NA   20     1     1 70  50
#> 67   3  24.000000 NA    NA   NA   NA     1     2 70  50
#> 68   3  34.666667 NA    NA   NA   NA     1     2 70  50
#> 69   3  36.000000 NA    NA   NA   NA     1     2 70  50
#> 70   3  45.333333 NA    NA   NA   NA     1     2 70  50
#> 71   3  56.000000 NA    NA   NA   NA     1     2 70  50
#> 72   3  66.666667 NA    NA   NA   NA     1     2 70  50
#> 73   3  72.000000 NA    NA   NA   NA     1     2 70  50
#> 74   3  77.333333 NA    NA   NA   NA     1     2 70  50
#> 75   3  88.000000 NA    NA   NA   NA     1     2 70  50
#> 76   3  98.666667 NA    NA   NA   NA     1     2 70  50
#> 77   3 109.333333 NA    NA   NA   NA     1     2 70  50
#> 78   3 120.000000 NA    NA   NA   NA     1     2 70  50
#> 79   4   0.500000 NA    NA   NA 1000     2     1 60  45
#> 80   4   0.500000 NA    NA   NA   NA     2     1 60  45
#> 81   4   1.000000 NA    NA   NA   NA     2     1 60  45
#> 82   4   1.111111 NA    NA   NA   NA     2     1 60  45
#> 83   4   1.722222 NA    NA   NA   NA     2     1 60  45
#> 84   4   2.000000 NA    NA   NA   NA     2     1 60  45
#> 85   4   2.333333 NA    NA   NA   NA     2     1 60  45
#> 86   4   2.944444 NA    NA   NA   NA     2     1 60  45
#> 87   4   3.555556 NA    NA   NA   NA     2     1 60  45
#> 88   4   4.166667 NA    NA   NA   NA     2     1 60  45
#> 89   4   4.777778 NA    NA   NA   NA     2     1 60  45
#> 90   4   5.388889 NA    NA   NA   NA     2     1 60  45
#> 91   4   6.000000 NA    NA   NA   NA     2     1 60  45
#> 92   4  10.000000 NA    NA   NA   20     2     1 60  45
#> 93   4  24.000000 NA    NA   NA   NA     2     2 60  45
#> 94   4  34.666667 NA    NA   NA   NA     2     2 60  45
#> 95   4  36.000000 NA    NA   NA   NA     2     2 60  45
#> 96   4  45.333333 NA    NA   NA   NA     2     2 60  45
#> 97   4  56.000000 NA    NA   NA   NA     2     2 60  45
#> 98   4  66.666667 NA    NA   NA   NA     2     2 60  45
#> 99   4  72.000000 NA    NA   NA   NA     2     2 60  45
#> 100  4  77.333333 NA    NA   NA   NA     2     2 60  45
#> 101  4  88.000000 NA    NA   NA   NA     2     2 60  45
#> 102  4  98.666667 NA    NA   NA   NA     2     2 60  45
#> 103  4 109.333333 NA    NA   NA   NA     2     2 60  45
#> 104  4 120.000000 NA    NA   NA   NA     2     2 60  45
#> 105  5   0.500000 NA    NA   NA 1000     2     1 60  45
#> 106  5   0.500000 NA    NA   NA   NA     2     1 60  45
#> 107  5   1.000000 NA    NA   NA   NA     2     1 60  45
#> 108  5   1.111111 NA    NA   NA   NA     2     1 60  45
#> 109  5   1.722222 NA    NA   NA   NA     2     1 60  45
#> 110  5   2.000000 NA    NA   NA   NA     2     1 60  45
#> 111  5   2.333333 NA    NA   NA   NA     2     1 60  45
#> 112  5   2.944444 NA    NA   NA   NA     2     1 60  45
#> 113  5   3.555556 NA    NA   NA   NA     2     1 60  45
#> 114  5   4.166667 NA    NA   NA   NA     2     1 60  45
#> 115  5   4.777778 NA    NA   NA   NA     2     1 60  45
#> 116  5   5.388889 NA    NA   NA   NA     2     1 60  45
#> 117  5   6.000000 NA    NA   NA   NA     2     1 60  45
#> 118  5  10.000000 NA    NA   NA   20     2     1 60  45
#> 119  5  24.000000 NA    NA   NA   NA     2     2 60  45
#> 120  5  34.666667 NA    NA   NA   NA     2     2 60  45
#> 121  5  36.000000 NA    NA   NA   NA     2     2 60  45
#> 122  5  45.333333 NA    NA   NA   NA     2     2 60  45
#> 123  5  56.000000 NA    NA   NA   NA     2     2 60  45
#> 124  5  66.666667 NA    NA   NA   NA     2     2 60  45
#> 125  5  72.000000 NA    NA   NA   NA     2     2 60  45
#> 126  5  77.333333 NA    NA   NA   NA     2     2 60  45
#> 127  5  88.000000 NA    NA   NA   NA     2     2 60  45
#> 128  5  98.666667 NA    NA   NA   NA     2     2 60  45
#> 129  5 109.333333 NA    NA   NA   NA     2     2 60  45
#> 130  5 120.000000 NA    NA   NA   NA     2     2 60  45
#> 131  6   0.500000 NA    NA   NA 1000     2     1 60  45
#> 132  6   0.500000 NA    NA   NA   NA     2     1 60  45
#> 133  6   1.000000 NA    NA   NA   NA     2     1 60  45
#> 134  6   1.111111 NA    NA   NA   NA     2     1 60  45
#> 135  6   1.722222 NA    NA   NA   NA     2     1 60  45
#> 136  6   2.000000 NA    NA   NA   NA     2     1 60  45
#> 137  6   2.333333 NA    NA   NA   NA     2     1 60  45
#> 138  6   2.944444 NA    NA   NA   NA     2     1 60  45
#> 139  6   3.555556 NA    NA   NA   NA     2     1 60  45
#> 140  6   4.166667 NA    NA   NA   NA     2     1 60  45
#> 141  6   4.777778 NA    NA   NA   NA     2     1 60  45
#> 142  6   5.388889 NA    NA   NA   NA     2     1 60  45
#> 143  6   6.000000 NA    NA   NA   NA     2     1 60  45
#> 144  6  10.000000 NA    NA   NA   20     2     1 60  45
#> 145  6  24.000000 NA    NA   NA   NA     2     2 60  45
#> 146  6  34.666667 NA    NA   NA   NA     2     2 60  45
#> 147  6  36.000000 NA    NA   NA   NA     2     2 60  45
#> 148  6  45.333333 NA    NA   NA   NA     2     2 60  45
#> 149  6  56.000000 NA    NA   NA   NA     2     2 60  45
#> 150  6  66.666667 NA    NA   NA   NA     2     2 60  45
#> 151  6  72.000000 NA    NA   NA   NA     2     2 60  45
#> 152  6  77.333333 NA    NA   NA   NA     2     2 60  45
#> 153  6  88.000000 NA    NA   NA   NA     2     2 60  45
#> 154  6  98.666667 NA    NA   NA   NA     2     2 60  45
#> 155  6 109.333333 NA    NA   NA   NA     2     2 60  45
#> 156  6 120.000000 NA    NA   NA   NA     2     2 60  45
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = 10,
                 model_minxt=10,model_maxxt=100)
#>     ID  Time DV IPRED PRED  AMT Group Model WT AGE
#> 1    1   0.5 NA    NA   NA 1000     1     1 70  50
#> 2    1   0.5 NA    NA   NA   NA     1     1 70  50
#> 3    1   1.0 NA    NA   NA   NA     1     1 70  50
#> 4    1   2.0 NA    NA   NA   NA     1     1 70  50
#> 5    1   6.0 NA    NA   NA   NA     1     1 70  50
#> 6    1  10.0 NA    NA   NA   20     1     1 70  50
#> 7    1  10.0 NA    NA   NA   NA     1     1 70  50
#> 8    1  10.0 NA    NA   NA   NA     1     2 70  50
#> 9    1  20.0 NA    NA   NA   NA     1     1 70  50
#> 10   1  20.0 NA    NA   NA   NA     1     2 70  50
#> 11   1  24.0 NA    NA   NA   NA     1     2 70  50
#> 12   1  30.0 NA    NA   NA   NA     1     1 70  50
#> 13   1  30.0 NA    NA   NA   NA     1     2 70  50
#> 14   1  36.0 NA    NA   NA   NA     1     2 70  50
#> 15   1  40.0 NA    NA   NA   NA     1     1 70  50
#> 16   1  40.0 NA    NA   NA   NA     1     2 70  50
#> 17   1  50.0 NA    NA   NA   NA     1     1 70  50
#> 18   1  50.0 NA    NA   NA   NA     1     2 70  50
#> 19   1  60.0 NA    NA   NA   NA     1     1 70  50
#> 20   1  60.0 NA    NA   NA   NA     1     2 70  50
#> 21   1  70.0 NA    NA   NA   NA     1     1 70  50
#> 22   1  70.0 NA    NA   NA   NA     1     2 70  50
#> 23   1  72.0 NA    NA   NA   NA     1     2 70  50
#> 24   1  80.0 NA    NA   NA   NA     1     1 70  50
#> 25   1  80.0 NA    NA   NA   NA     1     2 70  50
#> 26   1  90.0 NA    NA   NA   NA     1     1 70  50
#> 27   1  90.0 NA    NA   NA   NA     1     2 70  50
#> 28   1 100.0 NA    NA   NA   NA     1     1 70  50
#> 29   1 100.0 NA    NA   NA   NA     1     2 70  50
#> 30   1 120.0 NA    NA   NA   NA     1     2 70  50
#> 31   2   0.5 NA    NA   NA 1000     1     1 70  50
#> 32   2   0.5 NA    NA   NA   NA     1     1 70  50
#> 33   2   1.0 NA    NA   NA   NA     1     1 70  50
#> 34   2   2.0 NA    NA   NA   NA     1     1 70  50
#> 35   2   6.0 NA    NA   NA   NA     1     1 70  50
#> 36   2  10.0 NA    NA   NA   20     1     1 70  50
#> 37   2  10.0 NA    NA   NA   NA     1     1 70  50
#> 38   2  10.0 NA    NA   NA   NA     1     2 70  50
#> 39   2  20.0 NA    NA   NA   NA     1     1 70  50
#> 40   2  20.0 NA    NA   NA   NA     1     2 70  50
#> 41   2  24.0 NA    NA   NA   NA     1     2 70  50
#> 42   2  30.0 NA    NA   NA   NA     1     1 70  50
#> 43   2  30.0 NA    NA   NA   NA     1     2 70  50
#> 44   2  36.0 NA    NA   NA   NA     1     2 70  50
#> 45   2  40.0 NA    NA   NA   NA     1     1 70  50
#> 46   2  40.0 NA    NA   NA   NA     1     2 70  50
#> 47   2  50.0 NA    NA   NA   NA     1     1 70  50
#> 48   2  50.0 NA    NA   NA   NA     1     2 70  50
#> 49   2  60.0 NA    NA   NA   NA     1     1 70  50
#> 50   2  60.0 NA    NA   NA   NA     1     2 70  50
#> 51   2  70.0 NA    NA   NA   NA     1     1 70  50
#> 52   2  70.0 NA    NA   NA   NA     1     2 70  50
#> 53   2  72.0 NA    NA   NA   NA     1     2 70  50
#> 54   2  80.0 NA    NA   NA   NA     1     1 70  50
#> 55   2  80.0 NA    NA   NA   NA     1     2 70  50
#> 56   2  90.0 NA    NA   NA   NA     1     1 70  50
#> 57   2  90.0 NA    NA   NA   NA     1     2 70  50
#> 58   2 100.0 NA    NA   NA   NA     1     1 70  50
#> 59   2 100.0 NA    NA   NA   NA     1     2 70  50
#> 60   2 120.0 NA    NA   NA   NA     1     2 70  50
#> 61   3   0.5 NA    NA   NA 1000     1     1 70  50
#> 62   3   0.5 NA    NA   NA   NA     1     1 70  50
#> 63   3   1.0 NA    NA   NA   NA     1     1 70  50
#> 64   3   2.0 NA    NA   NA   NA     1     1 70  50
#> 65   3   6.0 NA    NA   NA   NA     1     1 70  50
#> 66   3  10.0 NA    NA   NA   20     1     1 70  50
#> 67   3  10.0 NA    NA   NA   NA     1     1 70  50
#> 68   3  10.0 NA    NA   NA   NA     1     2 70  50
#> 69   3  20.0 NA    NA   NA   NA     1     1 70  50
#> 70   3  20.0 NA    NA   NA   NA     1     2 70  50
#> 71   3  24.0 NA    NA   NA   NA     1     2 70  50
#> 72   3  30.0 NA    NA   NA   NA     1     1 70  50
#> 73   3  30.0 NA    NA   NA   NA     1     2 70  50
#> 74   3  36.0 NA    NA   NA   NA     1     2 70  50
#> 75   3  40.0 NA    NA   NA   NA     1     1 70  50
#> 76   3  40.0 NA    NA   NA   NA     1     2 70  50
#> 77   3  50.0 NA    NA   NA   NA     1     1 70  50
#> 78   3  50.0 NA    NA   NA   NA     1     2 70  50
#> 79   3  60.0 NA    NA   NA   NA     1     1 70  50
#> 80   3  60.0 NA    NA   NA   NA     1     2 70  50
#> 81   3  70.0 NA    NA   NA   NA     1     1 70  50
#> 82   3  70.0 NA    NA   NA   NA     1     2 70  50
#> 83   3  72.0 NA    NA   NA   NA     1     2 70  50
#> 84   3  80.0 NA    NA   NA   NA     1     1 70  50
#> 85   3  80.0 NA    NA   NA   NA     1     2 70  50
#> 86   3  90.0 NA    NA   NA   NA     1     1 70  50
#> 87   3  90.0 NA    NA   NA   NA     1     2 70  50
#> 88   3 100.0 NA    NA   NA   NA     1     1 70  50
#> 89   3 100.0 NA    NA   NA   NA     1     2 70  50
#> 90   3 120.0 NA    NA   NA   NA     1     2 70  50
#> 91   4   0.5 NA    NA   NA 1000     2     1 60  45
#> 92   4   0.5 NA    NA   NA   NA     2     1 60  45
#> 93   4   1.0 NA    NA   NA   NA     2     1 60  45
#> 94   4   2.0 NA    NA   NA   NA     2     1 60  45
#> 95   4   6.0 NA    NA   NA   NA     2     1 60  45
#> 96   4  10.0 NA    NA   NA   20     2     1 60  45
#> 97   4  10.0 NA    NA   NA   NA     2     1 60  45
#> 98   4  10.0 NA    NA   NA   NA     2     2 60  45
#> 99   4  20.0 NA    NA   NA   NA     2     1 60  45
#> 100  4  20.0 NA    NA   NA   NA     2     2 60  45
#> 101  4  24.0 NA    NA   NA   NA     2     2 60  45
#> 102  4  30.0 NA    NA   NA   NA     2     1 60  45
#> 103  4  30.0 NA    NA   NA   NA     2     2 60  45
#> 104  4  36.0 NA    NA   NA   NA     2     2 60  45
#> 105  4  40.0 NA    NA   NA   NA     2     1 60  45
#> 106  4  40.0 NA    NA   NA   NA     2     2 60  45
#> 107  4  50.0 NA    NA   NA   NA     2     1 60  45
#> 108  4  50.0 NA    NA   NA   NA     2     2 60  45
#> 109  4  60.0 NA    NA   NA   NA     2     1 60  45
#> 110  4  60.0 NA    NA   NA   NA     2     2 60  45
#> 111  4  70.0 NA    NA   NA   NA     2     1 60  45
#> 112  4  70.0 NA    NA   NA   NA     2     2 60  45
#> 113  4  72.0 NA    NA   NA   NA     2     2 60  45
#> 114  4  80.0 NA    NA   NA   NA     2     1 60  45
#> 115  4  80.0 NA    NA   NA   NA     2     2 60  45
#> 116  4  90.0 NA    NA   NA   NA     2     1 60  45
#> 117  4  90.0 NA    NA   NA   NA     2     2 60  45
#> 118  4 100.0 NA    NA   NA   NA     2     1 60  45
#> 119  4 100.0 NA    NA   NA   NA     2     2 60  45
#> 120  4 120.0 NA    NA   NA   NA     2     2 60  45
#> 121  5   0.5 NA    NA   NA 1000     2     1 60  45
#> 122  5   0.5 NA    NA   NA   NA     2     1 60  45
#> 123  5   1.0 NA    NA   NA   NA     2     1 60  45
#> 124  5   2.0 NA    NA   NA   NA     2     1 60  45
#> 125  5   6.0 NA    NA   NA   NA     2     1 60  45
#> 126  5  10.0 NA    NA   NA   20     2     1 60  45
#> 127  5  10.0 NA    NA   NA   NA     2     1 60  45
#> 128  5  10.0 NA    NA   NA   NA     2     2 60  45
#> 129  5  20.0 NA    NA   NA   NA     2     1 60  45
#> 130  5  20.0 NA    NA   NA   NA     2     2 60  45
#> 131  5  24.0 NA    NA   NA   NA     2     2 60  45
#> 132  5  30.0 NA    NA   NA   NA     2     1 60  45
#> 133  5  30.0 NA    NA   NA   NA     2     2 60  45
#> 134  5  36.0 NA    NA   NA   NA     2     2 60  45
#> 135  5  40.0 NA    NA   NA   NA     2     1 60  45
#> 136  5  40.0 NA    NA   NA   NA     2     2 60  45
#> 137  5  50.0 NA    NA   NA   NA     2     1 60  45
#> 138  5  50.0 NA    NA   NA   NA     2     2 60  45
#> 139  5  60.0 NA    NA   NA   NA     2     1 60  45
#> 140  5  60.0 NA    NA   NA   NA     2     2 60  45
#> 141  5  70.0 NA    NA   NA   NA     2     1 60  45
#> 142  5  70.0 NA    NA   NA   NA     2     2 60  45
#> 143  5  72.0 NA    NA   NA   NA     2     2 60  45
#> 144  5  80.0 NA    NA   NA   NA     2     1 60  45
#> 145  5  80.0 NA    NA   NA   NA     2     2 60  45
#> 146  5  90.0 NA    NA   NA   NA     2     1 60  45
#> 147  5  90.0 NA    NA   NA   NA     2     2 60  45
#> 148  5 100.0 NA    NA   NA   NA     2     1 60  45
#> 149  5 100.0 NA    NA   NA   NA     2     2 60  45
#> 150  5 120.0 NA    NA   NA   NA     2     2 60  45
#> 151  6   0.5 NA    NA   NA 1000     2     1 60  45
#> 152  6   0.5 NA    NA   NA   NA     2     1 60  45
#> 153  6   1.0 NA    NA   NA   NA     2     1 60  45
#> 154  6   2.0 NA    NA   NA   NA     2     1 60  45
#> 155  6   6.0 NA    NA   NA   NA     2     1 60  45
#> 156  6  10.0 NA    NA   NA   20     2     1 60  45
#> 157  6  10.0 NA    NA   NA   NA     2     1 60  45
#> 158  6  10.0 NA    NA   NA   NA     2     2 60  45
#> 159  6  20.0 NA    NA   NA   NA     2     1 60  45
#> 160  6  20.0 NA    NA   NA   NA     2     2 60  45
#> 161  6  24.0 NA    NA   NA   NA     2     2 60  45
#> 162  6  30.0 NA    NA   NA   NA     2     1 60  45
#> 163  6  30.0 NA    NA   NA   NA     2     2 60  45
#> 164  6  36.0 NA    NA   NA   NA     2     2 60  45
#> 165  6  40.0 NA    NA   NA   NA     2     1 60  45
#> 166  6  40.0 NA    NA   NA   NA     2     2 60  45
#> 167  6  50.0 NA    NA   NA   NA     2     1 60  45
#> 168  6  50.0 NA    NA   NA   NA     2     2 60  45
#> 169  6  60.0 NA    NA   NA   NA     2     1 60  45
#> 170  6  60.0 NA    NA   NA   NA     2     2 60  45
#> 171  6  70.0 NA    NA   NA   NA     2     1 60  45
#> 172  6  70.0 NA    NA   NA   NA     2     2 60  45
#> 173  6  72.0 NA    NA   NA   NA     2     2 60  45
#> 174  6  80.0 NA    NA   NA   NA     2     1 60  45
#> 175  6  80.0 NA    NA   NA   NA     2     2 60  45
#> 176  6  90.0 NA    NA   NA   NA     2     1 60  45
#> 177  6  90.0 NA    NA   NA   NA     2     2 60  45
#> 178  6 100.0 NA    NA   NA   NA     2     1 60  45
#> 179  6 100.0 NA    NA   NA   NA     2     2 60  45
#> 180  6 120.0 NA    NA   NA   NA     2     2 60  45
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = 10,
                 model_minxt=c(20,20),model_maxxt=c(100,100))
#>     ID      Time DV IPRED PRED  AMT Group Model WT AGE
#> 1    1   0.50000 NA    NA   NA 1000     1     1 70  50
#> 2    1   0.50000 NA    NA   NA   NA     1     1 70  50
#> 3    1   1.00000 NA    NA   NA   NA     1     1 70  50
#> 4    1   2.00000 NA    NA   NA   NA     1     1 70  50
#> 5    1   6.00000 NA    NA   NA   NA     1     1 70  50
#> 6    1  10.00000 NA    NA   NA   20     1     1 70  50
#> 7    1  20.00000 NA    NA   NA   NA     1     1 70  50
#> 8    1  20.00000 NA    NA   NA   NA     1     2 70  50
#> 9    1  24.00000 NA    NA   NA   NA     1     2 70  50
#> 10   1  28.88889 NA    NA   NA   NA     1     1 70  50
#> 11   1  28.88889 NA    NA   NA   NA     1     2 70  50
#> 12   1  36.00000 NA    NA   NA   NA     1     2 70  50
#> 13   1  37.77778 NA    NA   NA   NA     1     1 70  50
#> 14   1  37.77778 NA    NA   NA   NA     1     2 70  50
#> 15   1  46.66667 NA    NA   NA   NA     1     1 70  50
#> 16   1  46.66667 NA    NA   NA   NA     1     2 70  50
#> 17   1  55.55556 NA    NA   NA   NA     1     1 70  50
#> 18   1  55.55556 NA    NA   NA   NA     1     2 70  50
#> 19   1  64.44444 NA    NA   NA   NA     1     1 70  50
#> 20   1  64.44444 NA    NA   NA   NA     1     2 70  50
#> 21   1  72.00000 NA    NA   NA   NA     1     2 70  50
#> 22   1  73.33333 NA    NA   NA   NA     1     1 70  50
#> 23   1  73.33333 NA    NA   NA   NA     1     2 70  50
#> 24   1  82.22222 NA    NA   NA   NA     1     1 70  50
#> 25   1  82.22222 NA    NA   NA   NA     1     2 70  50
#> 26   1  91.11111 NA    NA   NA   NA     1     1 70  50
#> 27   1  91.11111 NA    NA   NA   NA     1     2 70  50
#> 28   1 100.00000 NA    NA   NA   NA     1     1 70  50
#> 29   1 100.00000 NA    NA   NA   NA     1     2 70  50
#> 30   1 120.00000 NA    NA   NA   NA     1     2 70  50
#> 31   2   0.50000 NA    NA   NA 1000     1     1 70  50
#> 32   2   0.50000 NA    NA   NA   NA     1     1 70  50
#> 33   2   1.00000 NA    NA   NA   NA     1     1 70  50
#> 34   2   2.00000 NA    NA   NA   NA     1     1 70  50
#> 35   2   6.00000 NA    NA   NA   NA     1     1 70  50
#> 36   2  10.00000 NA    NA   NA   20     1     1 70  50
#> 37   2  20.00000 NA    NA   NA   NA     1     1 70  50
#> 38   2  20.00000 NA    NA   NA   NA     1     2 70  50
#> 39   2  24.00000 NA    NA   NA   NA     1     2 70  50
#> 40   2  28.88889 NA    NA   NA   NA     1     1 70  50
#> 41   2  28.88889 NA    NA   NA   NA     1     2 70  50
#> 42   2  36.00000 NA    NA   NA   NA     1     2 70  50
#> 43   2  37.77778 NA    NA   NA   NA     1     1 70  50
#> 44   2  37.77778 NA    NA   NA   NA     1     2 70  50
#> 45   2  46.66667 NA    NA   NA   NA     1     1 70  50
#> 46   2  46.66667 NA    NA   NA   NA     1     2 70  50
#> 47   2  55.55556 NA    NA   NA   NA     1     1 70  50
#> 48   2  55.55556 NA    NA   NA   NA     1     2 70  50
#> 49   2  64.44444 NA    NA   NA   NA     1     1 70  50
#> 50   2  64.44444 NA    NA   NA   NA     1     2 70  50
#> 51   2  72.00000 NA    NA   NA   NA     1     2 70  50
#> 52   2  73.33333 NA    NA   NA   NA     1     1 70  50
#> 53   2  73.33333 NA    NA   NA   NA     1     2 70  50
#> 54   2  82.22222 NA    NA   NA   NA     1     1 70  50
#> 55   2  82.22222 NA    NA   NA   NA     1     2 70  50
#> 56   2  91.11111 NA    NA   NA   NA     1     1 70  50
#> 57   2  91.11111 NA    NA   NA   NA     1     2 70  50
#> 58   2 100.00000 NA    NA   NA   NA     1     1 70  50
#> 59   2 100.00000 NA    NA   NA   NA     1     2 70  50
#> 60   2 120.00000 NA    NA   NA   NA     1     2 70  50
#> 61   3   0.50000 NA    NA   NA 1000     1     1 70  50
#> 62   3   0.50000 NA    NA   NA   NA     1     1 70  50
#> 63   3   1.00000 NA    NA   NA   NA     1     1 70  50
#> 64   3   2.00000 NA    NA   NA   NA     1     1 70  50
#> 65   3   6.00000 NA    NA   NA   NA     1     1 70  50
#> 66   3  10.00000 NA    NA   NA   20     1     1 70  50
#> 67   3  20.00000 NA    NA   NA   NA     1     1 70  50
#> 68   3  20.00000 NA    NA   NA   NA     1     2 70  50
#> 69   3  24.00000 NA    NA   NA   NA     1     2 70  50
#> 70   3  28.88889 NA    NA   NA   NA     1     1 70  50
#> 71   3  28.88889 NA    NA   NA   NA     1     2 70  50
#> 72   3  36.00000 NA    NA   NA   NA     1     2 70  50
#> 73   3  37.77778 NA    NA   NA   NA     1     1 70  50
#> 74   3  37.77778 NA    NA   NA   NA     1     2 70  50
#> 75   3  46.66667 NA    NA   NA   NA     1     1 70  50
#> 76   3  46.66667 NA    NA   NA   NA     1     2 70  50
#> 77   3  55.55556 NA    NA   NA   NA     1     1 70  50
#> 78   3  55.55556 NA    NA   NA   NA     1     2 70  50
#> 79   3  64.44444 NA    NA   NA   NA     1     1 70  50
#> 80   3  64.44444 NA    NA   NA   NA     1     2 70  50
#> 81   3  72.00000 NA    NA   NA   NA     1     2 70  50
#> 82   3  73.33333 NA    NA   NA   NA     1     1 70  50
#> 83   3  73.33333 NA    NA   NA   NA     1     2 70  50
#> 84   3  82.22222 NA    NA   NA   NA     1     1 70  50
#> 85   3  82.22222 NA    NA   NA   NA     1     2 70  50
#> 86   3  91.11111 NA    NA   NA   NA     1     1 70  50
#> 87   3  91.11111 NA    NA   NA   NA     1     2 70  50
#> 88   3 100.00000 NA    NA   NA   NA     1     1 70  50
#> 89   3 100.00000 NA    NA   NA   NA     1     2 70  50
#> 90   3 120.00000 NA    NA   NA   NA     1     2 70  50
#> 91   4   0.50000 NA    NA   NA 1000     2     1 60  45
#> 92   4   0.50000 NA    NA   NA   NA     2     1 60  45
#> 93   4   1.00000 NA    NA   NA   NA     2     1 60  45
#> 94   4   2.00000 NA    NA   NA   NA     2     1 60  45
#> 95   4   6.00000 NA    NA   NA   NA     2     1 60  45
#> 96   4  10.00000 NA    NA   NA   20     2     1 60  45
#> 97   4  20.00000 NA    NA   NA   NA     2     1 60  45
#> 98   4  20.00000 NA    NA   NA   NA     2     2 60  45
#> 99   4  24.00000 NA    NA   NA   NA     2     2 60  45
#> 100  4  28.88889 NA    NA   NA   NA     2     1 60  45
#> 101  4  28.88889 NA    NA   NA   NA     2     2 60  45
#> 102  4  36.00000 NA    NA   NA   NA     2     2 60  45
#> 103  4  37.77778 NA    NA   NA   NA     2     1 60  45
#> 104  4  37.77778 NA    NA   NA   NA     2     2 60  45
#> 105  4  46.66667 NA    NA   NA   NA     2     1 60  45
#> 106  4  46.66667 NA    NA   NA   NA     2     2 60  45
#> 107  4  55.55556 NA    NA   NA   NA     2     1 60  45
#> 108  4  55.55556 NA    NA   NA   NA     2     2 60  45
#> 109  4  64.44444 NA    NA   NA   NA     2     1 60  45
#> 110  4  64.44444 NA    NA   NA   NA     2     2 60  45
#> 111  4  72.00000 NA    NA   NA   NA     2     2 60  45
#> 112  4  73.33333 NA    NA   NA   NA     2     1 60  45
#> 113  4  73.33333 NA    NA   NA   NA     2     2 60  45
#> 114  4  82.22222 NA    NA   NA   NA     2     1 60  45
#> 115  4  82.22222 NA    NA   NA   NA     2     2 60  45
#> 116  4  91.11111 NA    NA   NA   NA     2     1 60  45
#> 117  4  91.11111 NA    NA   NA   NA     2     2 60  45
#> 118  4 100.00000 NA    NA   NA   NA     2     1 60  45
#> 119  4 100.00000 NA    NA   NA   NA     2     2 60  45
#> 120  4 120.00000 NA    NA   NA   NA     2     2 60  45
#> 121  5   0.50000 NA    NA   NA 1000     2     1 60  45
#> 122  5   0.50000 NA    NA   NA   NA     2     1 60  45
#> 123  5   1.00000 NA    NA   NA   NA     2     1 60  45
#> 124  5   2.00000 NA    NA   NA   NA     2     1 60  45
#> 125  5   6.00000 NA    NA   NA   NA     2     1 60  45
#> 126  5  10.00000 NA    NA   NA   20     2     1 60  45
#> 127  5  20.00000 NA    NA   NA   NA     2     1 60  45
#> 128  5  20.00000 NA    NA   NA   NA     2     2 60  45
#> 129  5  24.00000 NA    NA   NA   NA     2     2 60  45
#> 130  5  28.88889 NA    NA   NA   NA     2     1 60  45
#> 131  5  28.88889 NA    NA   NA   NA     2     2 60  45
#> 132  5  36.00000 NA    NA   NA   NA     2     2 60  45
#> 133  5  37.77778 NA    NA   NA   NA     2     1 60  45
#> 134  5  37.77778 NA    NA   NA   NA     2     2 60  45
#> 135  5  46.66667 NA    NA   NA   NA     2     1 60  45
#> 136  5  46.66667 NA    NA   NA   NA     2     2 60  45
#> 137  5  55.55556 NA    NA   NA   NA     2     1 60  45
#> 138  5  55.55556 NA    NA   NA   NA     2     2 60  45
#> 139  5  64.44444 NA    NA   NA   NA     2     1 60  45
#> 140  5  64.44444 NA    NA   NA   NA     2     2 60  45
#> 141  5  72.00000 NA    NA   NA   NA     2     2 60  45
#> 142  5  73.33333 NA    NA   NA   NA     2     1 60  45
#> 143  5  73.33333 NA    NA   NA   NA     2     2 60  45
#> 144  5  82.22222 NA    NA   NA   NA     2     1 60  45
#> 145  5  82.22222 NA    NA   NA   NA     2     2 60  45
#> 146  5  91.11111 NA    NA   NA   NA     2     1 60  45
#> 147  5  91.11111 NA    NA   NA   NA     2     2 60  45
#> 148  5 100.00000 NA    NA   NA   NA     2     1 60  45
#> 149  5 100.00000 NA    NA   NA   NA     2     2 60  45
#> 150  5 120.00000 NA    NA   NA   NA     2     2 60  45
#> 151  6   0.50000 NA    NA   NA 1000     2     1 60  45
#> 152  6   0.50000 NA    NA   NA   NA     2     1 60  45
#> 153  6   1.00000 NA    NA   NA   NA     2     1 60  45
#> 154  6   2.00000 NA    NA   NA   NA     2     1 60  45
#> 155  6   6.00000 NA    NA   NA   NA     2     1 60  45
#> 156  6  10.00000 NA    NA   NA   20     2     1 60  45
#> 157  6  20.00000 NA    NA   NA   NA     2     1 60  45
#> 158  6  20.00000 NA    NA   NA   NA     2     2 60  45
#> 159  6  24.00000 NA    NA   NA   NA     2     2 60  45
#> 160  6  28.88889 NA    NA   NA   NA     2     1 60  45
#> 161  6  28.88889 NA    NA   NA   NA     2     2 60  45
#> 162  6  36.00000 NA    NA   NA   NA     2     2 60  45
#> 163  6  37.77778 NA    NA   NA   NA     2     1 60  45
#> 164  6  37.77778 NA    NA   NA   NA     2     2 60  45
#> 165  6  46.66667 NA    NA   NA   NA     2     1 60  45
#> 166  6  46.66667 NA    NA   NA   NA     2     2 60  45
#> 167  6  55.55556 NA    NA   NA   NA     2     1 60  45
#> 168  6  55.55556 NA    NA   NA   NA     2     2 60  45
#> 169  6  64.44444 NA    NA   NA   NA     2     1 60  45
#> 170  6  64.44444 NA    NA   NA   NA     2     2 60  45
#> 171  6  72.00000 NA    NA   NA   NA     2     2 60  45
#> 172  6  73.33333 NA    NA   NA   NA     2     1 60  45
#> 173  6  73.33333 NA    NA   NA   NA     2     2 60  45
#> 174  6  82.22222 NA    NA   NA   NA     2     1 60  45
#> 175  6  82.22222 NA    NA   NA   NA     2     2 60  45
#> 176  6  91.11111 NA    NA   NA   NA     2     1 60  45
#> 177  6  91.11111 NA    NA   NA   NA     2     2 60  45
#> 178  6 100.00000 NA    NA   NA   NA     2     1 60  45
#> 179  6 100.00000 NA    NA   NA   NA     2     2 60  45
#> 180  6 120.00000 NA    NA   NA   NA     2     2 60  45
model_prediction(design=design_4,DV=TRUE,dosing=dosing_4,model_num_points = c(10,10),
                 model_minxt=c(20,20),model_maxxt=c(100,100))
#>     ID      Time DV IPRED PRED  AMT Group Model WT AGE
#> 1    1   0.50000 NA    NA   NA 1000     1     1 70  50
#> 2    1   0.50000 NA    NA   NA   NA     1     1 70  50
#> 3    1   1.00000 NA    NA   NA   NA     1     1 70  50
#> 4    1   2.00000 NA    NA   NA   NA     1     1 70  50
#> 5    1   6.00000 NA    NA   NA   NA     1     1 70  50
#> 6    1  10.00000 NA    NA   NA   20     1     1 70  50
#> 7    1  20.00000 NA    NA   NA   NA     1     1 70  50
#> 8    1  20.00000 NA    NA   NA   NA     1     2 70  50
#> 9    1  24.00000 NA    NA   NA   NA     1     2 70  50
#> 10   1  28.88889 NA    NA   NA   NA     1     1 70  50
#> 11   1  28.88889 NA    NA   NA   NA     1     2 70  50
#> 12   1  36.00000 NA    NA   NA   NA     1     2 70  50
#> 13   1  37.77778 NA    NA   NA   NA     1     1 70  50
#> 14   1  37.77778 NA    NA   NA   NA     1     2 70  50
#> 15   1  46.66667 NA    NA   NA   NA     1     1 70  50
#> 16   1  46.66667 NA    NA   NA   NA     1     2 70  50
#> 17   1  55.55556 NA    NA   NA   NA     1     1 70  50
#> 18   1  55.55556 NA    NA   NA   NA     1     2 70  50
#> 19   1  64.44444 NA    NA   NA   NA     1     1 70  50
#> 20   1  64.44444 NA    NA   NA   NA     1     2 70  50
#> 21   1  72.00000 NA    NA   NA   NA     1     2 70  50
#> 22   1  73.33333 NA    NA   NA   NA     1     1 70  50
#> 23   1  73.33333 NA    NA   NA   NA     1     2 70  50
#> 24   1  82.22222 NA    NA   NA   NA     1     1 70  50
#> 25   1  82.22222 NA    NA   NA   NA     1     2 70  50
#> 26   1  91.11111 NA    NA   NA   NA     1     1 70  50
#> 27   1  91.11111 NA    NA   NA   NA     1     2 70  50
#> 28   1 100.00000 NA    NA   NA   NA     1     1 70  50
#> 29   1 100.00000 NA    NA   NA   NA     1     2 70  50
#> 30   1 120.00000 NA    NA   NA   NA     1     2 70  50
#> 31   2   0.50000 NA    NA   NA 1000     1     1 70  50
#> 32   2   0.50000 NA    NA   NA   NA     1     1 70  50
#> 33   2   1.00000 NA    NA   NA   NA     1     1 70  50
#> 34   2   2.00000 NA    NA   NA   NA     1     1 70  50
#> 35   2   6.00000 NA    NA   NA   NA     1     1 70  50
#> 36   2  10.00000 NA    NA   NA   20     1     1 70  50
#> 37   2  20.00000 NA    NA   NA   NA     1     1 70  50
#> 38   2  20.00000 NA    NA   NA   NA     1     2 70  50
#> 39   2  24.00000 NA    NA   NA   NA     1     2 70  50
#> 40   2  28.88889 NA    NA   NA   NA     1     1 70  50
#> 41   2  28.88889 NA    NA   NA   NA     1     2 70  50
#> 42   2  36.00000 NA    NA   NA   NA     1     2 70  50
#> 43   2  37.77778 NA    NA   NA   NA     1     1 70  50
#> 44   2  37.77778 NA    NA   NA   NA     1     2 70  50
#> 45   2  46.66667 NA    NA   NA   NA     1     1 70  50
#> 46   2  46.66667 NA    NA   NA   NA     1     2 70  50
#> 47   2  55.55556 NA    NA   NA   NA     1     1 70  50
#> 48   2  55.55556 NA    NA   NA   NA     1     2 70  50
#> 49   2  64.44444 NA    NA   NA   NA     1     1 70  50
#> 50   2  64.44444 NA    NA   NA   NA     1     2 70  50
#> 51   2  72.00000 NA    NA   NA   NA     1     2 70  50
#> 52   2  73.33333 NA    NA   NA   NA     1     1 70  50
#> 53   2  73.33333 NA    NA   NA   NA     1     2 70  50
#> 54   2  82.22222 NA    NA   NA   NA     1     1 70  50
#> 55   2  82.22222 NA    NA   NA   NA     1     2 70  50
#> 56   2  91.11111 NA    NA   NA   NA     1     1 70  50
#> 57   2  91.11111 NA    NA   NA   NA     1     2 70  50
#> 58   2 100.00000 NA    NA   NA   NA     1     1 70  50
#> 59   2 100.00000 NA    NA   NA   NA     1     2 70  50
#> 60   2 120.00000 NA    NA   NA   NA     1     2 70  50
#> 61   3   0.50000 NA    NA   NA 1000     1     1 70  50
#> 62   3   0.50000 NA    NA   NA   NA     1     1 70  50
#> 63   3   1.00000 NA    NA   NA   NA     1     1 70  50
#> 64   3   2.00000 NA    NA   NA   NA     1     1 70  50
#> 65   3   6.00000 NA    NA   NA   NA     1     1 70  50
#> 66   3  10.00000 NA    NA   NA   20     1     1 70  50
#> 67   3  20.00000 NA    NA   NA   NA     1     1 70  50
#> 68   3  20.00000 NA    NA   NA   NA     1     2 70  50
#> 69   3  24.00000 NA    NA   NA   NA     1     2 70  50
#> 70   3  28.88889 NA    NA   NA   NA     1     1 70  50
#> 71   3  28.88889 NA    NA   NA   NA     1     2 70  50
#> 72   3  36.00000 NA    NA   NA   NA     1     2 70  50
#> 73   3  37.77778 NA    NA   NA   NA     1     1 70  50
#> 74   3  37.77778 NA    NA   NA   NA     1     2 70  50
#> 75   3  46.66667 NA    NA   NA   NA     1     1 70  50
#> 76   3  46.66667 NA    NA   NA   NA     1     2 70  50
#> 77   3  55.55556 NA    NA   NA   NA     1     1 70  50
#> 78   3  55.55556 NA    NA   NA   NA     1     2 70  50
#> 79   3  64.44444 NA    NA   NA   NA     1     1 70  50
#> 80   3  64.44444 NA    NA   NA   NA     1     2 70  50
#> 81   3  72.00000 NA    NA   NA   NA     1     2 70  50
#> 82   3  73.33333 NA    NA   NA   NA     1     1 70  50
#> 83   3  73.33333 NA    NA   NA   NA     1     2 70  50
#> 84   3  82.22222 NA    NA   NA   NA     1     1 70  50
#> 85   3  82.22222 NA    NA   NA   NA     1     2 70  50
#> 86   3  91.11111 NA    NA   NA   NA     1     1 70  50
#> 87   3  91.11111 NA    NA   NA   NA     1     2 70  50
#> 88   3 100.00000 NA    NA   NA   NA     1     1 70  50
#> 89   3 100.00000 NA    NA   NA   NA     1     2 70  50
#> 90   3 120.00000 NA    NA   NA   NA     1     2 70  50
#> 91   4   0.50000 NA    NA   NA 1000     2     1 60  45
#> 92   4   0.50000 NA    NA   NA   NA     2     1 60  45
#> 93   4   1.00000 NA    NA   NA   NA     2     1 60  45
#> 94   4   2.00000 NA    NA   NA   NA     2     1 60  45
#> 95   4   6.00000 NA    NA   NA   NA     2     1 60  45
#> 96   4  10.00000 NA    NA   NA   20     2     1 60  45
#> 97   4  20.00000 NA    NA   NA   NA     2     1 60  45
#> 98   4  20.00000 NA    NA   NA   NA     2     2 60  45
#> 99   4  24.00000 NA    NA   NA   NA     2     2 60  45
#> 100  4  28.88889 NA    NA   NA   NA     2     1 60  45
#> 101  4  28.88889 NA    NA   NA   NA     2     2 60  45
#> 102  4  36.00000 NA    NA   NA   NA     2     2 60  45
#> 103  4  37.77778 NA    NA   NA   NA     2     1 60  45
#> 104  4  37.77778 NA    NA   NA   NA     2     2 60  45
#> 105  4  46.66667 NA    NA   NA   NA     2     1 60  45
#> 106  4  46.66667 NA    NA   NA   NA     2     2 60  45
#> 107  4  55.55556 NA    NA   NA   NA     2     1 60  45
#> 108  4  55.55556 NA    NA   NA   NA     2     2 60  45
#> 109  4  64.44444 NA    NA   NA   NA     2     1 60  45
#> 110  4  64.44444 NA    NA   NA   NA     2     2 60  45
#> 111  4  72.00000 NA    NA   NA   NA     2     2 60  45
#> 112  4  73.33333 NA    NA   NA   NA     2     1 60  45
#> 113  4  73.33333 NA    NA   NA   NA     2     2 60  45
#> 114  4  82.22222 NA    NA   NA   NA     2     1 60  45
#> 115  4  82.22222 NA    NA   NA   NA     2     2 60  45
#> 116  4  91.11111 NA    NA   NA   NA     2     1 60  45
#> 117  4  91.11111 NA    NA   NA   NA     2     2 60  45
#> 118  4 100.00000 NA    NA   NA   NA     2     1 60  45
#> 119  4 100.00000 NA    NA   NA   NA     2     2 60  45
#> 120  4 120.00000 NA    NA   NA   NA     2     2 60  45
#> 121  5   0.50000 NA    NA   NA 1000     2     1 60  45
#> 122  5   0.50000 NA    NA   NA   NA     2     1 60  45
#> 123  5   1.00000 NA    NA   NA   NA     2     1 60  45
#> 124  5   2.00000 NA    NA   NA   NA     2     1 60  45
#> 125  5   6.00000 NA    NA   NA   NA     2     1 60  45
#> 126  5  10.00000 NA    NA   NA   20     2     1 60  45
#> 127  5  20.00000 NA    NA   NA   NA     2     1 60  45
#> 128  5  20.00000 NA    NA   NA   NA     2     2 60  45
#> 129  5  24.00000 NA    NA   NA   NA     2     2 60  45
#> 130  5  28.88889 NA    NA   NA   NA     2     1 60  45
#> 131  5  28.88889 NA    NA   NA   NA     2     2 60  45
#> 132  5  36.00000 NA    NA   NA   NA     2     2 60  45
#> 133  5  37.77778 NA    NA   NA   NA     2     1 60  45
#> 134  5  37.77778 NA    NA   NA   NA     2     2 60  45
#> 135  5  46.66667 NA    NA   NA   NA     2     1 60  45
#> 136  5  46.66667 NA    NA   NA   NA     2     2 60  45
#> 137  5  55.55556 NA    NA   NA   NA     2     1 60  45
#> 138  5  55.55556 NA    NA   NA   NA     2     2 60  45
#> 139  5  64.44444 NA    NA   NA   NA     2     1 60  45
#> 140  5  64.44444 NA    NA   NA   NA     2     2 60  45
#> 141  5  72.00000 NA    NA   NA   NA     2     2 60  45
#> 142  5  73.33333 NA    NA   NA   NA     2     1 60  45
#> 143  5  73.33333 NA    NA   NA   NA     2     2 60  45
#> 144  5  82.22222 NA    NA   NA   NA     2     1 60  45
#> 145  5  82.22222 NA    NA   NA   NA     2     2 60  45
#> 146  5  91.11111 NA    NA   NA   NA     2     1 60  45
#> 147  5  91.11111 NA    NA   NA   NA     2     2 60  45
#> 148  5 100.00000 NA    NA   NA   NA     2     1 60  45
#> 149  5 100.00000 NA    NA   NA   NA     2     2 60  45
#> 150  5 120.00000 NA    NA   NA   NA     2     2 60  45
#> 151  6   0.50000 NA    NA   NA 1000     2     1 60  45
#> 152  6   0.50000 NA    NA   NA   NA     2     1 60  45
#> 153  6   1.00000 NA    NA   NA   NA     2     1 60  45
#> 154  6   2.00000 NA    NA   NA   NA     2     1 60  45
#> 155  6   6.00000 NA    NA   NA   NA     2     1 60  45
#> 156  6  10.00000 NA    NA   NA   20     2     1 60  45
#> 157  6  20.00000 NA    NA   NA   NA     2     1 60  45
#> 158  6  20.00000 NA    NA   NA   NA     2     2 60  45
#> 159  6  24.00000 NA    NA   NA   NA     2     2 60  45
#> 160  6  28.88889 NA    NA   NA   NA     2     1 60  45
#> 161  6  28.88889 NA    NA   NA   NA     2     2 60  45
#> 162  6  36.00000 NA    NA   NA   NA     2     2 60  45
#> 163  6  37.77778 NA    NA   NA   NA     2     1 60  45
#> 164  6  37.77778 NA    NA   NA   NA     2     2 60  45
#> 165  6  46.66667 NA    NA   NA   NA     2     1 60  45
#> 166  6  46.66667 NA    NA   NA   NA     2     2 60  45
#> 167  6  55.55556 NA    NA   NA   NA     2     1 60  45
#> 168  6  55.55556 NA    NA   NA   NA     2     2 60  45
#> 169  6  64.44444 NA    NA   NA   NA     2     1 60  45
#> 170  6  64.44444 NA    NA   NA   NA     2     2 60  45
#> 171  6  72.00000 NA    NA   NA   NA     2     2 60  45
#> 172  6  73.33333 NA    NA   NA   NA     2     1 60  45
#> 173  6  73.33333 NA    NA   NA   NA     2     2 60  45
#> 174  6  82.22222 NA    NA   NA   NA     2     1 60  45
#> 175  6  82.22222 NA    NA   NA   NA     2     2 60  45
#> 176  6  91.11111 NA    NA   NA   NA     2     1 60  45
#> 177  6  91.11111 NA    NA   NA   NA     2     2 60  45
#> 178  6 100.00000 NA    NA   NA   NA     2     1 60  45
#> 179  6 100.00000 NA    NA   NA   NA     2     2 60  45
#> 180  6 120.00000 NA    NA   NA   NA     2     2 60  45
```
