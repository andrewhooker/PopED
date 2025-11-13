# Evaluate the expectation of the Fisher Information Matrix (FIM) and the expectation of the OFV(FIM).

Compute the expectation of the FIM given the model, parameters,
distributions of parameter uncertainty, design and methods defined in
the PopED database.

## Usage

``` r
ed_mftot(
  model_switch,
  groupsize,
  ni,
  xtoptn,
  xoptn,
  aoptn,
  bpopdescr,
  ddescr,
  covd,
  sigma,
  docc,
  poped.db,
  calc_fim = TRUE,
  ...
)
```

## Arguments

- model_switch:

  A matrix that is the same size as xt, specifying which model each
  sample belongs to.

- groupsize:

  A vector of the number of individuals in each group.

- ni:

  A vector of the number of samples in each group.

- xtoptn:

  The xtoptn value

- xoptn:

  The xoptn

- aoptn:

  The aoptn value

- bpopdescr:

  Matrix defining the fixed effects, per row (row number =
  parameter_number) we should have:

  - column 1 the type of the distribution for E-family designs (0 =
    Fixed, 1 = Normal, 2 = Uniform, 3 = User Defined Distribution, 4 =
    lognormal and 5 = truncated normal)

  - column 2 defines the mean.

  - column 3 defines the variance of the distribution (or length of
    uniform distribution).

- ddescr:

  Matrix defining the diagonals of the IIV (same logic as for the
  `bpopdescr`).

- covd:

  Column major vector defining the covariances of the IIV variances.
  That is, from your full IIV matrix `covd <- IIV[lower.tri(IIV)]`.

- sigma:

  Matrix defining the variances can covariances of the residual
  variability terms of the model. can also just supply the diagonal
  parameter values (variances) as a
  [`c()`](https://rdrr.io/r/base/c.html).

- docc:

  Matrix defining the IOV, the IOV variances and the IOV distribution as
  for d and bpop.

- poped.db:

  A PopED database.

- calc_fim:

  Should the FIM be calculated or should we just use the user defined
  ed_penalty_pointer.

- ...:

  Other arguments passed to the function.

## Value

A list containing the E(FIM) and E(OFV(FIM)) and the a poped.db.

## See also

Other FIM:
[`LinMatrixH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixH.md),
[`LinMatrixLH()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixLH.md),
[`LinMatrixL_occ()`](https://andrewhooker.github.io/PopED/dev/reference/LinMatrixL_occ.md),
[`calc_ofv_and_fim()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_fim.md),
[`ed_laplace_ofv()`](https://andrewhooker.github.io/PopED/dev/reference/ed_laplace_ofv.md),
[`efficiency()`](https://andrewhooker.github.io/PopED/dev/reference/efficiency.md),
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md),
[`evaluate.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.fim.md),
[`gradf_eps()`](https://andrewhooker.github.io/PopED/dev/reference/gradf_eps.md),
[`mf3()`](https://andrewhooker.github.io/PopED/dev/reference/mf3.md),
[`mf7()`](https://andrewhooker.github.io/PopED/dev/reference/mf7.md),
[`mftot()`](https://andrewhooker.github.io/PopED/dev/reference/mftot.md),
[`ofv_criterion()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_criterion.md),
[`ofv_fim()`](https://andrewhooker.github.io/PopED/dev/reference/ofv_fim.md)

Other E-family:
[`calc_ofv_and_fim()`](https://andrewhooker.github.io/PopED/dev/reference/calc_ofv_and_fim.md),
[`ed_laplace_ofv()`](https://andrewhooker.github.io/PopED/dev/reference/ed_laplace_ofv.md),
[`evaluate.e.ofv.fim()`](https://andrewhooker.github.io/PopED/dev/reference/evaluate.e.ofv.fim.md)

## Examples

``` r
library(PopED)

############# START #################
## Create PopED database
## (warfarin model for optimization
##  with parameter uncertainty)
#####################################

## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

## Optimization using an additive + proportional reidual error
## to avoid sample times at very low concentrations (time 0 or very late samoples).

## find the parameters that are needed to define from the structural model
ff.PK.1.comp.oral.sd.CL
#> function (model_switch, xt, parameters, poped.db) 
#> {
#>     with(as.list(parameters), {
#>         y = xt
#>         y = (DOSE * Favail * KA/(V * (KA - CL/V))) * (exp(-CL/V * 
#>             xt) - exp(-KA * xt))
#>         return(list(y = y, poped.db = poped.db))
#>     })
#> }
#> <bytecode: 0x560952eadfa8>
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

# Adding 10% log-normal Uncertainty to fixed effects (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                         bpop_vals,
                         ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
bpop_vals_ed_ln
#>          bpop_vals         
#> CL     4      0.15 0.000225
#> V      4      8.00 0.640000
#> KA     4      1.00 0.010000
#> Favail 0      1.00 0.000000

## -- Define initial design  and design space
poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                  fg_fun=sfg,
                                  fError_fun=feps.add.prop,
                                  bpop=bpop_vals_ed_ln, 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=c( 0.5,1,2,6,24,36,72,120),
                                  minxt=0,
                                  maxxt=120,
                                  a=70,
                                  mina=0,
                                  maxa=100)

############# END ###################
## Create PopED database
## (warfarin model for optimization
##  with parameter uncertainty)
#####################################


# very few samples
poped.db$settings$ED_samp_size=10
ed_mftot(model_switch=poped.db$design$model_switch,
         groupsize=poped.db$design$groupsize,
         ni=poped.db$design$ni,
         xtoptn=poped.db$design$xt,
         xoptn=poped.db$design$x,
         aoptn=poped.db$design$a,
         bpopdescr=poped.db$parameters$bpop,
         ddescr=poped.db$parameters$d,
         covd=poped.db$parameters$covd,
         sigma=poped.db$parameters$sigma,
         docc=poped.db$parameters$docc, 
         poped.db)["ED_ofv"]
#> $ED_ofv
#> [1] 55.44415
#> 

  
```
