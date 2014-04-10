warfarin.design.1.red.input <- function(){


    popedInput <- list()

    ## --------------------------
    ## ---- Labeling and file names
    ## --------------------------
    
    ## -- Filname and path of the model file --
    popedInput$ff_file='ff'
    ## -- Filname and path of the parameter file --
    popedInput$fg_file='sfg'
    ## -- Filname and path of the residual error model file --
    popedInput$fError_file='feps.add.prop'
    ## -- The model title --
    popedInput$modtit='One comp first order absorption'
    
    ## --------------------------
    ## ---- Design variables ----
    ## --------------------------

    #     ## -- Number of discrete design variables --
    #     popedInput$nx=0
    ## -- Number of continuous design variables that are not time (e.g. continuous covariates) --
    #     ## dose in this example
    #     popedInput$na=1 
    ## -- Number of groups/individuals --
    #popedInput$m=4
    #     ## -- Max number of samples per group/individual --
     #    popedInput$maxni=8
    #     ## -- Min number of samples per group/individual --
    #     popedInput$minni=1
    ## -- Vector defining the size of the different groups (num individuals in each group) --
    popedInput$design$groupsize=rbind(8,8,8,8)
    #     ## -- Vector defining the max size of the different groups (max num individuals in each group) --
    #     popedInput$design$maxgroupsize=rbind(20,20,20,20)
    #     ## -- Vector defining the min size of the different groups (min num individuals in each group) --
    #     popedInput$design$mingroupsize = rbind(1, 1, 1, 1)
    #     ## -- The total maximal groupsize over all groups--
    #     popedInput$design$maxtotgroupsize=32
    #     ## -- The total minimal groupsize over all groups--
    #     popedInput$design$mintotgroupsize=4
    #     ## -- Vector defining the number of samples for each group --
    #     popedInput$design$ni=rbind(8,8,8,8)
    ## -- Matrix defining the initial sampling schedule --
    popedInput$design$xt=rbind(c(0.5,1,2,6,24,36,72,120),
                               c(0.5,1,2,6,24,36,72,120),
                               c(0.5,1,2,6,24,36,72,120),
                               c(0.5,1,2,6,24,36,72,120))
    #     ## -- Matrix defining the max value for each sample --
    #     popedInput$design$maxxt=rbind(c(25,25,25,120,120,120,120,120),
    #                                   c(25,25,25,120,120,120,120,120),
    #                                   c(25,25,25,120,120,120,120,120),
    #                                   c(25,25,25,120,120,120,120,120))
    #     ## -- Matrix defining the min value for each sample --
    #     popedInput$design$minxt=rbind(c(0,0,0,0,0,0,0,0),
    #                                   c(0,0,0,0,0,0,0,0),
    #                                   c(0,0,0,0,0,0,0,0),
    #                                   c(0,0,0,0,0,0,0,0))
    #     ## -- Matrix defining the initial discrete values --
    #     popedInput$design$x=matrix(0,4,0)
    #     ## -- Cell defining the discrete variables --
    #     popedInput$design$discrete_x=cell(4,0)
    ## -- Vector defining the initial covariate values --
    popedInput$design$a=rbind(70, 70, 70, 70)
    #     ## -- Vector defining the max value for each covariate --
    #     popedInput$design$maxa=rbind(100, 100, 100, 100)
    #     ## -- Vector defining the min value for each covariate --
    #     popedInput$design$mina=rbind(1,1,1,1)
    #     ## -- Vector defining which response a certain sampling time belongs to --
    #     popedInput$design$model_switch=rbind(c(1,1,1,1,1,1,1,1),
    #                                          c(1,1,1,1,1,1,1,1),
    #                                          c(1,1,1,1,1,1,1,1),
    #                                          c(1,1,1,1,1,1,1,1))
    #     ## -- Use grouped time points (1=TRUE, 0=FALSE) --
    #     popedInput$bUseGrouped_xt=0
    #     ## -- Matrix defining the grouping of sample points --
    #     popedInput$design$G=rbind(c(1,2,3,4,5,6,7,8),
    #                               c(9,10,11,12,13,14,15,16),
    #                               c(17,18,19,20,21,22,23,24),
    #                               c(25,26,27,28,29,30,31,32))
    #     ## -- Use grouped covariates (1=TRUE, 0=FALSE) --
    #     popedInput$bUseGrouped_a=0
    #     ## -- Matrix defining the grouping of covariates --
    #     popedInput$design$Ga=rbind(1, 2, 3, 4)
    #     ## -- Use grouped discrete design variables (1=TRUE, 0=FALSE) --
    #     popedInput$bUseGrouped_x=0
    #     ## -- Matrix defining the grouping of discrete design variables --
    #     popedInput$design$Gx=matrix(0,4,0)
    ## -- Value to interpret as zero in design --
    #popedInput$ourzero=1.0e-5


    ## --------------------------
    ## ---- Criterion specification
    ## --------------------------

    ## -- OFV calculation for FIM (1=Determinant of FIM,4=log determinant of FIM,6=determinant of interesting part of FIM (Ds)) --
    #popedInput$ofv_calc_type=1
    #     ## -- Approximation method, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI --
    #     popedInput$iApproximationMethod=0
    #     ## -- Num indivduals in each step of FOCE --
    #     popedInput$iFOCENumInd=1000    
    ## -- Fisher Information Matrix type
    ## (0=Full FIM,
    ##  1=Reduced FIM,
    ##  2=weighted models,
    ##  3=Loc models,
    ##  4=reduced FIM with derivative of SD of sigma as pfim,
    ##  5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
    ##  6=Calculate one model switch at a time, good for large matrices,
    ##  7=Reduced FIM parameterized with A,B,C matrices & derivative of variance) --
    #popedInput$iFIMCalculationType=4
    #     ## -- D-family design (1) or ED-familty design (0) (with or without parameter uncertainty) --
    #     popedInput$d_switch=1
    #     ## -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
    #     popedInput$iEDCalculationType=0    
    #     ## -- The prior FIM (added to calculated FIM) --
    #     popedInput$prior_fim=matrix(0,0,1)
    #     ## -- Penalty function, empty string means no penalty.  User defined criterion --
    #     popedInput$strEDPenaltyFile=''
    #     ## -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
    #     popedInput$strAutoCorrelationFile=''
    
    ## --------------------------
    ## ---- Model parameters and parmeter distributions and fixed or not
    ## --------------------------

    ## -- Matrix defining the fixed effects, per row (row number = parameter_number),
    ## the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
    ## 3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal).
    ## The second column defines the mean.
    ## The third column defines the variance of the distribution.
    popedInput$design$bpop=rbind(c(4,8,4),
                                 c(4,1,0.0625),
                                 c(4,0.15,0.0014),
                                 c(0,1,0))
    ## -- Matrix defining the diagnonals of the IIV (same logic as for the fixed efects) --
    popedInput$design$d=rbind(c(4,0.02,0.000025),
                              c(4,0.6,0.000225),
                              c(4,0.07,0.00030625))
    #     ## -- Matrix defining the covariances of the IIV variances --
    #     popedInput$design$covd=cbind(0,0,0) 
    ## -- Matrix defining the variances of the residual variability terms --
    popedInput$design$sigma <- diag(c(0.01,0.25))
    #     ## -- Matrix defining the IOV, the IOV variances and the IOV distribution --
    #     popedInput$design$docc=matrix(0,0,3)
    #     ## -- Matrix defining the covariance of the IOV --
    #     popedInput$design$covdocc=matrix(0,1,0)
    ## -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
    popedInput$notfixed_bpop=cbind(1,1,1,0)
    #     ## -- Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) --
     #    popedInput$notfixed_d=cbind(1,1,1)
    #     ## -- Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed) --
     #    popedInput$notfixed_covd=cbind(0,0,0)
    #     ## -- Vector defining if a IOV variance is fixed or not (1=not fixed, 0=fixed) --
    #     popedInput$notfixed_docc=matrix(0,1,0)
    #     ## -- Vector defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) --
    #     popedInput$notfixed_covdocc=matrix(0,1,0)
    ## -- Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) --
    #popedInput$notfixed_sigma <- cbind(1,1)
    ## -- Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed) --
    ##popedInput$notfixed_covsigma=matrix(0,1,0)
    #     ## -- Filname and path for user defined distributions for E-family designs --
    #     popedInput$strUserDistributionFile=''
    ## -- How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
    #popedInput$bLHS=0
    ## -- Sample size for E-family sampling --
    #popedInput$ED_samp_size=45
    ## -- Ds_index, set index to 1 if a parameter is uninteresting, otherwise 0.
    ## size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma --
    #popedInput$CriterionOptions$ds_index <- cbind(0,0,0,0,0,0,0,0)
    ## -- User defined data structure that, for example could be used to send in data to the model --
    #popedInput$user_data=cell(0,0)
    ## -- The length of the g parameter vector --
    #popedInput$ng=5
    ## -- Number of typical values --
    #popedInput$nbpop=4
    ## -- Number of IIV parameters --
    #popedInput$nb=3
    #     ## -- Number of IOV variance parameters --
    #     popedInput$ndocc=0
    #     ## -- Number of occassions --
    #     popedInput$NumOcc=0


    #     ## --------------------------
    #     ## ---- What to optimize
    #     ## --------------------------
    #     
    #     ## -- Vector of optimization tasks (1=TRUE,0=FALSE)
    #     ## (Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group)
    #     ## -- All elements set to zero => only calculate the FIM with current design --
    #     popedInput$optsw = cbind(0,0,0,0,0)
    

    ## --------------------------
    ## ---- Optimization algorithm settings
    ## --------------------------

    #     ## -- Use random search (1=TRUE, 0=FALSE) --
    #     popedInput$bUseRandomSearch=1
    #     ## -- Use Stochastic Gradient search (1=TRUE, 0=FALSE) --
    #     popedInput$bUseStochasticGradient=0
    #     ## -- Use Line search (1=TRUE, 0=FALSE) --
    #     popedInput$bUseLineSearch=0
    #     ## -- Use Exchange algorithm (1=TRUE, 0=FALSE) --
    #     popedInput$bUseExchangeAlgorithm=0
    #     ## -- Use BFGS Minimizer (1=TRUE, 0=FALSE) --
    #     popedInput$bUseBFGSMinimizer=0
    #     ## -- Vector for line search on continuous design variables (1=TRUE,0=FALSE) --
    #     popedInput$line_opta=1
    #     ## -- Vector for line search on discrete design variables (1=TRUE,0=FALSE) --
    #     popedInput$line_optx=matrix(0,0,1)
    #     ## -- The seed number used for optimization and sampling --
    #     popedInput$dSeed=-1
    ## -- Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov --
    #popedInput$EACriteria=1
    ## -- Filename and path for a run file that is used instead of the regular PopED call --
    #popedInput$strRunFile=''

    ## --------------------------
    ## ---- PopED specific
    ## --------------------------

    ## -- The current PopED version --
    #popedInput$strPopEDVersion='2.13'

    #     ## -- Use graph output during search --
    #     popedInput$bShowGraphs=0
    #     ## -- If a log file should be used (0=FALSE, 1=TRUE) --
    #     popedInput$use_logfile=0
    ## -- Filname and path of the output file during search --
    #popedInput$output_file=paste(match.call()[[1]],'_summary',sep='')
    ## -- Filname suffix of the result function file --
    #popedInput$output_function_file=paste(match.call()[[1]],'_output_',sep='')
    ## -- Filename and path for storage of current optimal design --
    #popedInput$strIterationFileName=paste(match.call()[[1]],'_current.R',sep='')
    #     ## -- Method used to calculate M1
    #     ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    #     popedInput$m1_switch=1
    #     ## -- Method used to calculate M2
    #     ## (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    #     popedInput$m2_switch=1
    #     ## -- Method used to calculate linearization of residual error
    #     ## (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
    #     popedInput$hle_switch=1
    #     ## -- Method used to calculate the gradient of the model
    #     ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    #     popedInput$gradff_switch=1
    #     ## -- Method used to calculate the gradient of the parameter vector g
    #     ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
    #     popedInput$gradfg_switch=1
    #     ## -- Number of iterations in random search between screen output --
    #     popedInput$rsit_output=5
    #     ## -- Number of iterations in stochastic gradient search between screen output --
    #     popedInput$sgit_output=1
    #     ## -- Step length of derivative of linearized model w.r.t. typical values --
    #     popedInput$hm1=0.0001
    #     ## -- Step length of derivative of model w.r.t. g --
    #     popedInput$hlf=0.0001
    #     ## -- Step length of derivative of g w.r.t. b --
    #     popedInput$hlg=0.0001
    #     ## -- Step length of derivative of variance w.r.t. typical values --
    #     popedInput$hm2=0.0001
    #     ## -- Step length of derivative of OFV w.r.t. time --
    #     popedInput$hgd=0.0001
    #     ## -- Step length of derivative of model w.r.t. sigma --
    #     popedInput$hle=0.0001
    #     ## -- The absolute tolerance for the diff equation solver --
    #     popedInput$AbsTol=1e-05
    #     ## -- The relative tolerance for the diff equation solver --
    #     popedInput$RelTol=1e-05
    #     ## -- The diff equation solver method, 0=ode45, 1=ode15s --
    #     popedInput$iDiffSolverMethod=0
    #     ## -- If the differential equation results should be stored in memory (1) or not (0) --
    #     popedInput$bUseMemorySolver=0
    #     ## -- Number of Random search iterations --
    #     popedInput$rsit=300
    #     ## -- Number of Stochastic gradient search iterations --
    #     popedInput$sgit=150
    #     ## -- Number of Random search iterations with discrete optimization --
    #     popedInput$intrsit=250
    #     ## -- Number of Stochastic Gradient search iterations with discrete optimization --
    #     popedInput$intsgit=50
    #     ## -- Iterations until adaptive narrowing in random search --
    #     popedInput$maxrsnullit=250
    #     ## -- Stoachstic Gradient convergence value,
    #     ## (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
    #     popedInput$convergence_eps=1e-08
    #     ## -- Random search locality factor for sample times --
    #     popedInput$rslxt=10
    #     ## -- Random search locality factor for covariates --
    #     popedInput$rsla=10
    #     ## -- Stochastic Gradient search first step factor for sample times --
    #     popedInput$cfaxt=0.001
    #     ## -- Stochastic Gradient search first step factor for covariates --
    #     popedInput$cfaa=0.001
    #     ## -- Use greedy algorithm for group assignment optimization --
    #     popedInput$bGreedyGroupOpt=0
    #     ## -- Exchange Algorithm StepSize --
    #     popedInput$EAStepSize=0.01
    #     ## -- Exchange Algorithm NumPoints --
    #     popedInput$EANumPoints=0
    #     ## -- Exchange Algorithm Convergence Limit/Criteria --
    #     popedInput$EAConvergenceCriteria=1e-20
    #     ## -- Avoid replicate samples when using Exchange Algorithm --
    #     popedInput$bEANoReplicates=0
    #     ## -- BFGS Minimizer Convergence Criteria Minimum Step --
    #     popedInput$BFGSConvergenceCriteriaMinStep=1e-08
    #     ## -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
    #     popedInput$BFGSProjectedGradientTol=0.0001
    #     ## -- BFGS Minimizer Line Search Tolerance f --
    #     popedInput$BFGSTolerancef=0.001
    #     ## -- BFGS Minimizer Line Search Tolerance g --
    #     popedInput$BFGSToleranceg=0.9
    #     ## -- BFGS Minimizer Line Search Tolerance x --
    #     popedInput$BFGSTolerancex=0.1
    #     ## -- Number of iterations in ED-optimal design to calculate convergence criteria --
    #     popedInput$ED_diff_it=30
    #     ## -- ED-optimal design convergence criteria in percent --
    #     popedInput$ED_diff_percent=10
    #     ## -- Number of grid points in the line search --
    #     popedInput$line_search_it=50
    #     ## -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
    #     popedInput$iNumSearchIterationsIfNotLineSearch=1
    
 

    #     ## -- Parallel options for PopED -- --
    #     ## -- Compile option for PopED
    #     ## -1 = No compilation,
    #     ## 0 or 3 = Full compilation,
    #     ## 1 or 4 = Only using MCC (shared lib),
    #     ## 2 or 5 = Only MPI,
    #     ## Option 0,1,2 runs PopED and option 3,4,5 stops after compilation --
    #     popedInput$parallelSettings$iCompileOption=-1
    #     ## -- Parallel method to use (0 = Matlab PCT, 1 = MPI) --
    #     popedInput$parallelSettings$iUseParallelMethod=1
    #     ## -- Additional dependencies used in MCC compilation (mat-files), if several space separated --
    #     popedInput$parallelSettings$strAdditionalMCCCompilerDependencies=''
    #     ## -- Compilation output executable name --
    #     popedInput$parallelSettings$strExecuteName='calc_fim.exe'
    #     ## -- Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) --
    #     popedInput$parallelSettings$iNumProcesses=2
    #     ## -- Number of design evaluations that should be evaluated in each process before getting new work from job manager --
    #     popedInput$parallelSettings$iNumChunkDesignEvals=-2
    #     ## -- The prefix of the input mat file to communicate with the excutable --
    #     popedInput$parallelSettings$strMatFileInputPrefix='parallel_input'
    #     ## -- The prefix of the output mat file to communicate with the excutable --
    #     popedInput$parallelSettings$strMatFileOutputPrefix='parallel_output'
    #     ## -- Extra options send to e$g. the MPI exectuable or a batch script, see execute_parallel$m for more information and options --
    #     popedInput$parallelSettings$strExtraRunOptions=''
    #     ## -- Polling time to check if the parallel execution is finished --
    #     popedInput$parallelSettings$dPollResultTime=1.000000e-01
    #     ## -- The file containing the popedInput structure that should be used to evaluate the designs --
    #     popedInput$parallelSettings$strFunctionInputName='function_input_fo_reduced'
    #     ## -- If the random search is going to be executed in parallel --
    #     #popedInput$parallelSettings$bParallelRS=0
    #     ## -- If the stochastic gradient search is going to be executed in parallel --
    #     #popedInput$parallelSettings$bParallelSG=0
    #     ## -- If the line search is going to be executed in parallel --
    #     #popedInput$parallelSettings$bParallelLS=0
    #     ## -- If the modified exchange algorithm is going to be executed in parallel --
    #     #popedInput$parallelSettings$bParallelMFEA=0

    return( popedInput ) 
}
