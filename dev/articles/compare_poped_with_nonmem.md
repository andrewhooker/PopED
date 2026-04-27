# Compare Uncertaitny Predictions from PopED with NONMEM

Packages used for this vignette

``` r
library(PopED)
library(tidyverse)
```

## 1 Population models

In this work we are using population, or nonlinear mixed-effect (NLME)
models. Here we define $y_{ij}$ for the $j^{th}$ observation of the
$i^{th}$ individual in a population as:

$$y_{ij} = f(t_{ij},{\overset{\rightarrow}{a}}_{i},{\overset{\rightarrow}{\theta}}_{i}) + h(t_{ij},{\overset{\rightarrow}{a}}_{i},{\overset{\rightarrow}{\theta}}_{i},{\overset{\rightarrow}{\varepsilon}}_{ij})\qquad(1)$$

Where $t_{ij}$ are the measurement times
${\overset{\rightarrow}{a}}_{i}$ is a vector of covariates (doses of a
drug, weight, age, concentration of a drug in blood plasma, etc.),
${\overset{\rightarrow}{\theta}}_{i}$ is a vector of model parameter
values, $h(.)$ is a model for the residual error (also referred to as
residual unexplained variability, or RUV) in the model and
${\overset{\rightarrow}{\varepsilon}}_{ij}$ is a vector of of random
variables describing data-level deviations from the model. Often, as is
the case in the models described here, the elements of
${\overset{\rightarrow}{\varepsilon}}_{ij}$ are assumed to come from
normal distributions with means of zero and a covariance matrix of
$\mathbf{\Sigma}$ (elements of $\sigma_{lm}^{2}$), where
$\mathbf{\Sigma}$ is typically diagonal.

Population effects are modeled on the parameter level, where individual
parameter values are derived from typical (or population) parameters
$\overset{\rightarrow}{\theta}$, individual deviations due to covariates
${\overset{\rightarrow}{a}}_{i}$, and random individual deviations
${\overset{\rightarrow}{\eta}}_{i}$ (referred to as a between-subject
variability, or BSV, term).

$${\overset{\rightarrow}{\theta}}_{i} = g(\overset{\rightarrow}{\theta},{\overset{\rightarrow}{a}}_{i},{\overset{\rightarrow}{\eta}}_{i})\qquad(2)$$

Extensions, where deviations are on other scales, are possible as well
(such as parameter deviations between occasions within an individual’s
study, center level deviations, study level deviations, etc.). Often, as
is the case in the models described here, the elements of
${\overset{\rightarrow}{\eta}}_{i}$ are assumed to come from normal
distributions with means of zero and a covariance matrix of
$\mathbf{\Omega}$ (elements of $\omega_{pq}^{2}$).

## 2 Simple Population PK model in PopED

### 2.1 Structural model

Here we define a one-compartment pharmacokinetic model with linear
absorption and a single drug dose using an analytical solution.

$$f(t_{ij},D_{i},{\overset{\rightarrow}{\theta}}_{i}) = \frac{D_{i} \cdot F \cdot Ka_{i}}{V_{i} \cdot (Ka_{i} - CL_{i}/V_{i})} \cdot (e^{-\frac{CL_{i}}{V_{i}}t_{ij}} - e^{-Ka_{i} \cdot t_{ij}})\qquad(3)$$

Where $D_{i}$ is the dose of drug given to individual $i$, $F$ is the
bioavailability, $Ka_{i}$ is the absorption rate constant for individual
$i$, $V_{i}$ is the volume of distribution for individual $i$, and
$CL_{i}$ is the clearance for individual $i$.

Defining this in PopED we have:

``` r
##-- Model: One comp first order absorption
ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list(y=y,poped.db=poped.db))
  })
}
```

### 2.2 Random effects model

Taking the random effect structure and parameter values from the
Warfarin example from software comparison in Nyberg et al., “Methods and
software tools for design evaluation for population
pharmacokinetics-pharmacodynamics studies”, Br. J. Clin. Pharm., 2014
([Nyberg et al. 2015](#ref-nybergMethodsSoftwareTools2015)), we have a
proportional residual error model, with a coefficient of variation of
10%, and exponential random effects are assumed for CL, V and Ka.

$$\begin{aligned}
y_{ij} & {= f(t_{ij},D_{i},{\overset{\rightarrow}{\theta}}_{i}) \cdot (1 + {\overset{\rightarrow}{\varepsilon}}_{ij})} \\
 & \\
{CL_{i}} & {= \theta_{CL} \cdot e^{\eta_{CL,i}}} \\
V_{i} & {= \theta_{V} \cdot e^{\eta_{V,i}}} \\
{Ka_{i}} & {= \theta_{Ka} \cdot e^{\eta_{Ka,i}}} \\
\mathbf{\Omega} & {= \begin{bmatrix}
\omega_{CL}^{2} & 0 & 0 \\
0 & \omega_{V}^{2} & 0 \\
0 & 0 & \omega_{Ka}^{2}
\end{bmatrix}} \\
\mathbf{\Sigma} & {= \begin{bmatrix}
\sigma_{\text{prop}}^{2}
\end{bmatrix}}
\end{aligned}\qquad(4)$$

This is defined in PopED with the following:

``` r
## -- parameter definition function 
sfg <- function(x,a,bpop,b,bocc){
  parameters=c(CL=bpop[1]*exp(b[1]),
               V=bpop[2]*exp(b[2]),
               KA=bpop[3]*exp(b[3]),
               Favail=bpop[4],
               DOSE=a[1])
  return(parameters) 
}

## -- Residual Error function
feps <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]] 
  y = y*(1+epsi[,1])  
  return(list(y=y,poped.db=poped.db)) 
}
```

### 2.3 Parameter values and study design

We create a poped database to link the model defined above with a set of
model parameters, the initial design and a design space for
optimization.

In this example, the parameter values are defined for the fixed effects
(`bpop`), the between-subject variability variances (`d`) and the
residual variability variances (`sigma`). We also fix the parameter
`Favail` using `notfixed_bpop`, since we have only oral dosing and the
parameter is not identifiable. Fixing a parameter means that we assume
the parameter will not be estimated (and is know without uncertainty).

For the initial design, we define 1 group (`m=1`) of 32 individuals
(`groupsize=32`), with a single dose of 70 mg (`a`). The initial design
has 8 sample times per individual (`xt`).

For the design space, which can be searched during design optimization,
we define a potential dose range of between 0 and 100 mg (`mina` and
`maxa`), and a range of potential sample times between 0 and 120 hours
(`minxt` and `maxxt`).

``` r
## -- Define initial design  and design space
poped.db <- 
  create.poped.database(
    # Model
    ff_fun = ff,
    fg_fun = sfg,
    fError_fun = feps,
                                  
    # Parameters 
    bpop = c(CL=0.15, V=8, KA=1.0, Favail=1), 
    notfixed_bpop = c(1,1,1,0),
    d = c(CL=0.07, V=0.02, KA=0.6), 
    sigma = c(0.01),
                                  
    # Design
    groupsize = 32,
    xt = c(0.5,1,2,6,24,36,72,120),
    a = 70,
    
    # Design space
    minxt = 0,
    maxxt = 120,
    mina = 0,
    maxa = 100)
```

### 2.4 Design simulation

First it may make sense to check your model and design to make sure you
get what you expect when simulating data. Here we plot the model typical
values:

``` r
plot_model_prediction(poped.db, model_num_points = 300)
```

![](compare_poped_with_nonmem_files/figure-html/simulate_without_BSV-1.png)

Next, we plot the expected prediction interval (by default a 95% PI) of
the data taking into account the BSV and RUV using the option `PI=TRUE`.
This option makes predictions based on first-order approximations to the
model variance and a normality assumption of that variance. Better (and
slower) computations are possible with the `DV=T`, `IPRED=T` and
`groupsize_sim = some large number` options.

``` r
plot_model_prediction(poped.db, 
                      PI=TRUE, 
                      model_num_points = 300, 
                      sample.times = FALSE)
```

![](compare_poped_with_nonmem_files/figure-html/simulate_with_BSV-1.png)

We can get these predictions numerically as well:

``` r
dat <- model_prediction(poped.db,DV=TRUE)
head(dat,n=5);tail(dat,n=5)
```

``` output
  ID Time       DV    IPRED     PRED Group Model a_i
1  1  0.5 2.886426 2.926416 3.425436     1     1  70
2  1  1.0 4.319344 4.780023 5.471104     1     1  70
3  1  2.0 5.461513 6.669522 7.382183     1     1  70
4  1  6.0 6.755186 7.527285 7.946280     1     1  70
5  1 24.0 6.104940 5.616734 5.685856     1     1  70
```

``` output
    ID Time        DV     IPRED      PRED Group Model a_i
252 32    6 7.2368395 8.7990243 7.9462805     1     1  70
253 32   24 3.9584873 4.7926859 5.6858561     1     1  70
254 32   36 3.1444957 3.1886025 4.5402483     1     1  70
255 32   72 0.9887762 0.9389979 2.3116966     1     1  70
256 32  120 0.1803847 0.1839714 0.9398657     1     1  70
```

### 2.5 Design evaluation

Next, we evaluate the initial design

``` r
(ds1 <- evaluate_design(poped.db))
```

``` output
$ofv
[1] 52.44799

$fim
                    CL          V        KA         d_CL          d_V
CL         19821.28445 -21.836551 -8.622140 0.000000e+00     0.000000
V            -21.83655  20.656071 -1.807099 0.000000e+00     0.000000
KA            -8.62214  -1.807099 51.729039 0.000000e+00     0.000000
d_CL           0.00000   0.000000  0.000000 3.107768e+03    10.728786
d_V            0.00000   0.000000  0.000000 1.072879e+01 27307.089308
d_KA           0.00000   0.000000  0.000000 2.613561e-02     3.265608
SIGMA[1,1]     0.00000   0.000000  0.000000 5.215403e+02 11214.210707
                  d_KA   SIGMA[1,1]
CL          0.00000000      0.00000
V           0.00000000      0.00000
KA          0.00000000      0.00000
d_CL        0.02613561    521.54030
d_V         3.26560786  11214.21071
d_KA       41.81083599     71.08764
SIGMA[1,1] 71.08763902 806176.95068

$rse
        CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1]
  4.738266   2.756206  13.925829  25.627205  30.344316  25.777327  11.170784 
```

We see that all parameters relative standard errors (rse above) in % are
predicted to be relatively well estimated.

## 3 Evaluation of design in NONMEM

### 3.1 Create dataset

Now we create a dataset simulated from this model and design

``` r
set.seed(12398)
dosing_1 <- list(list(AMT=70,Time=0))
dat <- model_prediction(poped.db, dosing = dosing_1, DV=T,filename = "./compare_poped_with_nonmem_data/temp.csv") |> 
  dplyr::select(-c(IPRED,PRED))
readr::write_csv(dat,file = "./compare_poped_with_nonmem_data/warfarin.csv",na = ".")
dat[1:20,]
```

``` output
   ID  Time       DV AMT Group Model a_i
1   1   0.0       NA  70     1     1  70
2   1   0.5 3.198757  NA     1     1  70
3   1   1.0 5.319627  NA     1     1  70
4   1   2.0 6.977056  NA     1     1  70
5   1   6.0 9.502303  NA     1     1  70
6   1  24.0 7.424704  NA     1     1  70
7   1  36.0 6.072961  NA     1     1  70
8   1  72.0 3.646986  NA     1     1  70
9   1 120.0 1.595916  NA     1     1  70
10  2   0.0       NA  70     1     1  70
11  2   0.5 4.997493  NA     1     1  70
12  2   1.0 6.519595  NA     1     1  70
13  2   2.0 8.211107  NA     1     1  70
14  2   6.0 7.138961  NA     1     1  70
15  2  24.0 5.664554  NA     1     1  70
16  2  36.0 5.059473  NA     1     1  70
17  2  72.0 2.647690  NA     1     1  70
18  2 120.0 1.227866  NA     1     1  70
19  3   0.0       NA  70     1     1  70
20  3   0.5 8.317745  NA     1     1  70
```

The simulated data look like this:

``` r
library(ggplot2)
p <- ggplot(dat,aes(x=Time,y=DV,group=ID)) + 
    geom_line(color="grey") + 
    geom_point(size=3,alpha=0.5,
               aes(colour = ID) ) + 
  theme(legend.position = "none")#+
    # labs(color="ID")
p
```

![](compare_poped_with_nonmem_files/figure-html/unnamed-chunk-8-1.png)

### 3.2 NONMEM evaluation

We create a NONMEM model file.

``` r
nm_mod <- paste(
"$PROBLEM
$INPUT ID TIME DV AMT GROUP MODEL DOSE
$DATA ../compare_poped_with_nonmem_data/warfarin.csv IGNORE=@
$SUBROUTINES ADVAN2  TRANS2 
$PK
  CL = THETA(1) * EXP(ETA(1))
  VC = THETA(2) * EXP(ETA(2))
  KA = THETA(3) * EXP(ETA(3))
  V = VC

$ERROR
  CONC = A(2)/VC
  Y = CONC + CONC * EPS(1)

$ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 POSTHOC 
$COVARIANCE

$THETA (0, 0.15) ; POP_CL
$THETA (0, 8) ; POP_VC
$THETA (0, 1) ; POP_KA

$OMEGA 0.07 ; IIV_CL
$OMEGA 0.02 ; IIV_VC
$OMEGA 0.6 ; IIV_KA

$SIGMA 0.01; RUV_PROP

$TABLE ID VC CL KA ETA(1) ETA(2) ETA(3) FILE=sdtab1 NOAPPEND NOPRINT
") |> glue::as_glue()

write_file(nm_mod,file="./compare_poped_with_nonmem_models/run1.mod")
#nm_mod
```

Then we estimate the model from the data in NONMEM. In PsN this would be

`execute ./compare_poped_with_nonmem_models/run1.mod`

and

`sumo ./compare_poped_with_nonmem_models/run1.lst`

For NONMEM, parameter RSE is:

``` r
ext.val <- read_table("./compare_poped_with_nonmem_models/run1.ext",skip=1,show_col_types = FALSE)
par_est_nonmem <- ext.val |> filter(ITERATION==-1000000000) |> 
  select(-c(OBJ,ITERATION,"OMEGA(2,1)","OMEGA(3,1)","OMEGA(3,2)")) 
par_se_nonmem <- ext.val |> filter(ITERATION==-1000000001) |>
  select(-c(OBJ,ITERATION,"OMEGA(2,1)","OMEGA(3,1)","OMEGA(3,2)")) 
par_rse_nonmem <- par_se_nonmem/par_est_nonmem*100
names(par_rse_nonmem) <- c("CL","V","KA","SIGMA[1,1]","d_CL","d_V","d_KA")
par_rse_nonmem <- par_rse_nonmem |> relocate("SIGMA[1,1]", .after = last_col())
round(par_rse_nonmem,2)
```

``` output
    CL    V    KA  d_CL   d_V  d_KA SIGMA[1,1]
1 4.58 2.28 16.63 26.81 38.72 24.14       9.02
```

## 4 Compare PopED and NONMEM

We can compare between the prediction in PopED and the covariance matrix
in NONMEM

``` r
result <- rbind(par_rse_nonmem,ds1$rse)|> 
  mutate("Method"=c("NONMEM","PopED"),.before=CL)
result
```

``` output
  Method       CL        V       KA     d_CL      d_V     d_KA SIGMA[1,1]
1 NONMEM 4.584171 2.283751 16.63198 26.80822 38.71585 24.14017   9.019811
2  PopED 4.738266 2.756206 13.92583 25.62720 30.34432 25.77733  11.170784
```

``` r
result |> 
  pivot_longer(-Method,names_to = "Parameter",values_to = "RSE") |> 
  ggplot(aes(x=Parameter,y=RSE,fill=Method)) + 
  geom_col(position = position_dodge()) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = names(result[-1]))
```

![](compare_poped_with_nonmem_files/figure-html/unnamed-chunk-12-1.png)

The NONMEM and PopED results are similar, though they differ in how they
are obtained. The NONMEM results are based on a single simulation of the
design and estimation of the model, whereas the PopED results reflect an
evaluation of the expected uncertainties, on average, across all
realizations of the data. Additionally, the PopED results rely on the
Fisher Information Matrix, which is an asymptotic approximation and
therefore represents an upper bound on the overall information. As a
result, the predicted RSE values from PopED are expected to be lower on
average (Cramér-Rao Lower Bound).

We could also compare the results of an SSE in NONMEM to the PopED
results, which would be more comparable since both would be based on
many realizations of the data. This was done in Nyberg et al.
([2015](#ref-nybergMethodsSoftwareTools2015)).

## References

Nyberg, Joakim, Caroline Bazzoli, Kay Ogungbenro, et al. 2015. “Methods
and Software Tools for Design Evaluation in Population
Pharmacokinetics-Pharmacodynamics Studies.” *British Journal of Clinical
Pharmacology* 79 (1): 6–17. <https://doi.org/10.1111/bcp.12352>.
