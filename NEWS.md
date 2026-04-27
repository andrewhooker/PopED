# PopED (development version)

# PopED 0.7.0

* `create.poped.database()` now uses a better method of identifying the total number of parameters of each type (bpop, d, sigma, etc.) in a user defined model parameter function (the `ff_fun` argument in `create.poped.database()`) (#73).

* `create.poped.database()` has a new option `reorder_parameter_vectors`, which is turned off by default. When turned on, if you use named arguments in `bpop` or `d` then PopED will try to figure out the order of the parameters based on what is found in the `fg_fun`. See the resulting `poped_db$parameters` and make sure the order matches with `fg_fun`.

* `start_parallel()` has a new default for`num_cores`, which is now one less than the number of cores available from `parallel::detectCores()`.

* `model_prediction()` and therefore `plot_model_prediction()` allow for log-normal distributions when using the PI option. This makes sense if you expect the prediction interval of the model will be approximately log-normally distributed, which might often be the case in pop PK models. The new default is now `PI_ln_dist = TRUE`.

* `poped_optim()` now has an explicit argument allowing for the specification of Ds-optimal parameters of interest. The option is `ds_index`.

* Minor bug fixes 


# PopED 0.6.0

* Added the options `allow_replicates=TRUE/FALSE`, `allow_replicates_xt=TRUE/FALSE` and `allow_replicates_a=TRUE/FALSE` to `poped_optim`.  This allows the optimization algorithm to avoid replicates (or not) in the design components.  Currently only works for discrete variable optimization. Future versions will also handle continuous optimization.

* Exported a function for the computation of the Bayesian Fisher information matrix for individual parameters of a population model based on Maximum A Posteriori (MAP) estimation of the empirical Bayes estimates (EBEs) in a population model. See `?evaluate_fim_map` for more details.

* Allowing for no covariates in the function that automatically builds a PopED parameter function from a model function (`?build_sfg`).

* Updates to documentation and package testing.

* Minor bug fixes.


# PopED 0.5.0

* Added the ability to incorporate limit of quantification information into FIM calculations (both upper and lower limits). See the new vignette on the webpage https://andrewhooker.github.io/PopED/articles/handling_loq.html 

* Adding functionality to optimize groupsize and total size of the study. See `?optimize_groupsize`, ``?optimize_n_eff` and `?optimize_n_rse`.  This is also implemented in `poped_optim` through the `opt_inds=T` argument.

* Updating Vignettes, including a new one about using other tools to use as simulators for design computations.  See https://andrewhooker.github.io/PopED/articles/model_def_other_pkgs.html 

* Simplify RxODE syntax in the above vingette (#47, @mattfidler).

* Added the ability to predict and plot model prediction intervals by computing the expected variance (using an FO approximation) and then computing a prediction interval based on an assumption of normality. See `?model_prediciton` and `?plot_model_prediction`. The computation is faster but less accurate compared to using `DV=TRUE` (and `groupsize_sim = 500`) in the two functions.

* Named parameters are now passed to all calculations so that the FIM and RSE output is more readable with parameter names instead of default names.

* Allow for parallel computation in `plot_efficiency_of_windows` (#50).

* Make parallelization work with mrgsolve on windows (#37, #45, #46, #51, @Vincent-AC).

* Updated the function for automatic building of parameter model function (see `build_sfg`).

* Simplify derivative calculations (#34, @martin-gmx).

* Allow for only simulating model_switch > 1 models.

* Change the defult Ds calculation to be on log scale.

* Updated the website at https://andrewhooker.github.io/PopED 

* Remove options for discontinued dplyr commands `rbind_all` and `rbind_list`.

* Minor bug fixes in shrinkage calculations (#44, #39, @martin-gmx).


# PopED 0.4.0

* New and improved vignettes (#30, @giulialestini)!

* Added power evaluation script to test the power of a design to identify a parameter different than an assumed value.  The function also calculates the number of individuals needed in a design to have a specific power. See `?evaluate_power` for more information (#26, @martin-gmx).

* Added function to compute expected shrinkage of a design.  See `?shrinkage` for more information.

* Updated and added new example scripts in `system.file("examples", package="PopED")` ().  This includes an example describing 
how to handle covariate distributions in optimal design, an example on how to incorporate IOV,  an example on how to handle shrinkage, an example with a full covariance matrix and an example with a prior FIM (#30, @giulialestini and @martin-gmx).

* Major overhaul in optimization methods used in `poped_optim()` so that generic optimization routines like `optim()` can be easily used in optimizing  PopED designs.  

* Update speed of FIM calculations (#20, @martin-gmx).

* Update RSE calculations so that prior FIM is handled correctly (#22, @martin-gmx).

* Simplified code and removed duplicated code (#21, #24 and #32, @martin-gmx).

* New ways of handling inverting matricies, should be faster and work better when the matricies are ill-conditioned. See `?inv` for more information (#19, @martin-gmx).

* Updated functionality of IOV calculations.

* Updates to `optim_ARS()` for when to stop search.

* Extended functionality of `plot_model_prediction()` (#23, @martin-gmx).

* Bug fixing.  See https://github.com/andrewhooker/PopED/commits/master for more information.

# PopED 0.3.2

* Exported the `summary` method for the results of `poped_optim` in the PopED NAMESPACE, so that the method can actually be used!  Just use `summary(output)`.

* Fixed some old bugs that used `return` as a varible in functions, a la MATLAB.


# PopED 0.3.1

* Added a vignette to introduce PopED!

* Improved optimization with `poped_optim`, plus all example scripts now running with `poped_optim`.

* Update to more easily allow discrete optimization of xt and a variables.  See the example scripts.

* Added a summary method for the results of `poped_optim`.  Just use `summary(output)`.

* changed handling of seed numbers in optimizations.

* more robust handling of non-population models

* more natural handling of NA values in design vectors

* NAMESPACE: removed ggplot2 from "Depends" and added to "Imports" 

* Added mean line to efficiency plots.

* Update to computation and error handling for Laplace approximation to ED objective function. 

* Added more intuitive cost function input.  See examples in `?poped_optim`

* Various small changes and bug fixes.


# PopED 0.3.0

* Added new optimization methods and tools, see `?poped_optim()`. This function incorporates the new optimization routines `optim_ARS()` and `optim_LS` which are optimized versions of previous optimization algorithms used in PopED. Both can be run with parallelization. `poped_optim()` also incorporates the genetic algorithm from `GA::ga()`, which can also be run with parallelization, and the "L-BFGS-B" method from `stats::optim()`. `poped_optim()` should be more efficient and faster than `poped_optimize()`.

* Changed the default objective function to be the log of the determinant of the FIM.  `create.poped.database(ofv_calc_type=4)`

* Various small changes and bug fixes.


# PopED 0.2.0

* Fixed `plot_efficiency_of_windows()` bug that had wrong headers on each subplot.

* Fixed bug in `plot_model_prediction()` that did not plot the optimized design, but instead the initial design

* Reorganized the database created from `create.poped.database()`.  The output from this function is now a list with 5 sub-lists: design, design_space, model, parameters and settings.  Also removed duplicate entries in the database for easier manipulation.  This will cause some back compatibility issues when referring to elements in a database.

* Added example 10 describing a PKPD design of hepatitis C virus (HCV) kinetics to the `system.file("examples",package="PopED")` directory of the PopED installation.

# PopED 0.1.2

* Updated model_prediction() to allow for creation of NONMEM datasets.  
  Useful for testing of optimized designs via PsN's (http://psn.sf.net) SSE tool, for example.

* Two new functions create_design() and create_design_space() that allow for design and design space creation without the 
  need for a model or parameter values.

* Updated the create.poped.database() function to use create_design() and create_design_space()

* Added examples for evaluation and optimization of a one-target quasi-steady-state 
  target mediated drug disposition model (TMDD) to the system.file("examples",package="PopED") directory of the PopED installation.

* Added a 2-compartment, oral absorption, multiple dose example to the system.file("examples", 
  package="PopED") directory of the PopED installation.

* Updated plot_efficiency_of_windows() to allow for the plotting of the RSE of each parameter on the y-axis.

* Updated error handing for the Laplace approximation of the ED OFV.

* Fixed bug when computing FIM with only one BSV term present in model (calculation gave 
  an error).

* Fixed a bug in plot_model_predictions where an error was returned if not all time 
  values in the xt matrix were to be used for the design calculation 
  (ni is different from size(xt,2), see ?create_poped_database).

* Various small bug fixes.


# PopED 0.1.1

* Updated package author list

* New functionality to compute the ED OFV using the Laplace approximation.
  This can be orders of magnitude faster than the standard MC integration approach.
  See '?ed_laplace_ofv' and '?evaluate.e.ofv.fim' 

* Added a general function to compute the FIM and OFV(FIM) for all available methods in PopED.
  See '?calc_ofv_and_fim'.

* Added a general optimization algorithm 'RS_opt_gen()' that works for both D-family and 
  E-family optimization.

* Added optimization of E-family designs to 'poped_optimize()'.

* Changed distribution tests for package building 

* Fixed bug where correlations between BSV (between subject variability) terms in the model gave an 
  error when creating a PopED database

* Fixed a bug where get_rse failed when a parameter had a value of 3.


# PopED 0.1.0
 
* PopED has been translated to R from MATLAB and this is the initial release.
