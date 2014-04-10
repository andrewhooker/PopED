poped.db.2 <- create.poped.database(list(),
                                    groupsize=20,
                                    sigma=diag(c( 0.15,0.15)),
                                    bpop=c( 0.5,0.2,1,1,1),  # should be in model definition and be named values
                                    d=c( 0.01,0.01,0.01,0.01,0.01), # should be in model definition and be named values
                                    xt=c( 0.33,0.66,0.25,0.5,0.75),
                                    maxxt=c( 1,1,1,1,1), # should associated with each xt more tightly
                                    minxt=c( 0,0,0,0,0), # should associated with each xt more tightly
                                    a=2.75,
                                    maxa=5,# should associated with each xt more tightly
                                    mina=0.5,# should associated with each xt more tightly
                                    model_switch=c( 1,1,2,2,2),# should associated with each xt more tightly
                                    iFIMCalculationType=0,
                                    modtit='One comp IV bolus with direct Emax effect' ## should be in model definition here there should be a design title
)  




one.comp.emax.1 <- function(){
  
  popedInput <- list()
  
  # -- Auto generated input file for the run One comp IV bolus with direct Emax effect --
  # -- Created with PopED version 2.13, 31-Jan-2014 10:11:55
  
  # -- The current PopED version --
  popedInput$strPopEDVersion='2.13'
  
  # ---- Design variables ----
  
  # -- The length of the g parameter vector --
  popedInput$ng=6
  # -- Number of typical values --
  popedInput$nbpop=5
  # -- Number of IIV parameters --
  popedInput$nb=5
  # -- Number of IOV variance parameters --
  popedInput$ndocc=0
  # -- Number of discrete variables --
  popedInput$nx=0
  # -- Number of covariates --
  popedInput$na=1
  # -- Number of occassions --
  popedInput$NumOcc=0
  
  # -- Number of individuals/groups --
  popedInput$m=1
  # -- Max number of samples per group/individual --
  popedInput$maxni=5
  # -- Min number of samples per group/individual --
  popedInput$minni=1
  
  # -- D-optimal design (1) or ED-optimal design (0) --
  popedInput$d_switch=1
  # -- Approximation method, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI --
  popedInput$iApproximationMethod=0
  # -- Num indivduals in each conditional step --
  popedInput$iFOCENumInd=1000
  # -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
  popedInput$iEDCalculationType=0
  # -- Use random search (1=TRUE, 0=FALSE) --
  popedInput$bUseRandomSearch=1
  # -- Use Stochastic Gradient search (1=TRUE, 0=FALSE) --
  popedInput$bUseStochasticGradient=1
  # -- Use Line search (1=TRUE, 0=FALSE) --
  popedInput$bUseLineSearch=1
  
  # -- Use Exchange algorithm (1=TRUE, 0=FALSE) --
  popedInput$bUseExchangeAlgorithm=0
  
  # -- Use BFGS Minimizer (1=TRUE, 0=FALSE) --
  popedInput$bUseBFGSMinimizer=0
  
  # -- OFV calculation for FIM (1=Determinant of FIM,4=log determinant of FIM,6=determinant of interesting part of FIM (Ds)) --
  popedInput$ofv_calc_type=1
  # -- The prior FIM --
  popedInput$prior_fim=t(zeros(1,0))
  
  
  # -- Vector of optimization tasks (1=TRUE,0=FALSE) matrix(c(Samples per subject,Sampling schedule,Discrete design var,Covariates,Num ind per group),nrow=1,byrow=T) --
  # -- All elements set to zero => only calculate the FIM with current design --
  popedInput$optsw=matrix(c( 0,1,0,1,0),nrow=1,byrow=T)
  
  # -- Vector for line search on covariates (1=TRUE,0=FALSE) --
  popedInput$line_opta=1
  
  # -- Vector for line search on discrete design variables (1=TRUE,0=FALSE) --
  popedInput$line_optx=t(zeros(1,0))
  
  
  # -- The seed number used for optimization and sampling --
  popedInput$dSeed=-1
  
  # -- Vector defining the size of the different groups (num individuals in each group) --
  popedInput$design$groupsize=20
  
  # -- Vector defining the max size of the different groups (max num individuals in each group) --
  popedInput$design$maxgroupsize=20
  
  # -- Vector defining the min size of the different groups (min num individuals in each group) --
  popedInput$design$mingroupsize=20
  
  
  # -- The total maximal groupsize --
  popedInput$design$maxtotgroupsize=20
  
  # -- The total minimal groupsize --
  popedInput$design$mintotgroupsize=20
  
  # -- Vector defining the variances of the errors --
  popedInput$design$sigma=matrix(c( 0.15,0,
                                    0,0.15),nrow=2,byrow=T)
  
  
  # -- Matrix defining the typical values, the typical value variances and the distribution --
  popedInput$design$bpop=matrix(c( 0,0.5,0,
                                   0,0.2,0,
                                   0,1,0,
                                   0,1,0,
                                   0,1,0),nrow=5,byrow=T)
  
  # -- Matrix defining the IIV, the IIV variances and the IIV distribution --
  popedInput$design$d=matrix(c( 0,0.01,0,
                                0,0.01,0,
                                0,0.01,0,
                                0,0.01,0,
                                0,0.01,0),nrow=5,byrow=T)
  
  # -- Matrix defining the covariances of the IIV variances --
  popedInput$design$covd=matrix(c( 0,0,0,0,0,0,0,0,0,0),nrow=1,byrow=T)
  
  # -- Matrix defining the IOV, the IOV variances and the IOV distribution --
  popedInput$design$docc=t(zeros(3,0))
  
  # -- Matrix defining the covariance of the IOV --
  popedInput$design$covdocc=t(zeros(0,1))
  
  
  # -- Vector defining the number of sample for each group --
  popedInput$design$ni=5
  
  
  # -- Matrix defining the initial sampling schedule --
  popedInput$design$xt=matrix(c( 0.33,0.66,0.25,0.5,0.75),nrow=1,byrow=T)
  
  # -- Matrix defining the max value for each sample --
  popedInput$design$maxxt=matrix(c( 1,1,1,1,1),nrow=1,byrow=T)
  
  # -- Matrix defining the min value for each sample --
  popedInput$design$minxt=matrix(c( 0,0,0,0,0),nrow=1,byrow=T)
  
  
  # -- Matrix defining the initial discrete values --
  popedInput$design$x=t(zeros(0,1))
  
  # -- Cell defining the discrete variables --
  popedInput$design$discrete_x=t(cell(0,1))
  
  
  # -- Vector defining the initial covariate values --
  popedInput$design$a=2.75
  
  # -- Vector defining the max value for each covariate --
  popedInput$design$maxa=5
  
  # -- Vector defining the min value for each covariate --
  popedInput$design$mina=0.5
  
  
  # -- Vector defining which response a certain sampling time belongs to --
  popedInput$design$model_switch=matrix(c( 1,1,2,2,2),nrow=1,byrow=T)
  
  
  # -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_bpop=matrix(c( 1,1,1,1,1),nrow=1,byrow=T)
  
  # -- Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_d=matrix(c( 1,1,1,1,1),nrow=1,byrow=T)
  
  # -- Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_covd=matrix(c( 0,0,0,0,0,0,0,0,0,0),nrow=1,byrow=T)
  
  # -- Vector defining if a IOV variance is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_docc=t(zeros(0,1))
  
  # -- Vector defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_covdocc=t(zeros(0,1))
  
  # -- Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_sigma=matrix(c( 1,1),nrow=1,byrow=T)
  
  # -- Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed) --
  popedInput$notfixed_covsigma=0
  
  
  # -- Use grouped time points (1=TRUE, 0=FALSE) --
  popedInput$bUseGrouped_xt=0
  # -- Matrix defining the grouping of sample points --
  popedInput$design$G=matrix(c( 1,2,3,4,5),nrow=1,byrow=T)
  
  # -- Use grouped covariates (1=TRUE, 0=FALSE) --
  popedInput$bUseGrouped_a=0
  # -- Matrix defining the grouping of covariates --
  popedInput$design$Ga=1
  
  # -- Use grouped discrete design variables (1=TRUE, 0=FALSE) --
  popedInput$bUseGrouped_x=0
  # -- Matrix defining the grouping of discrete design variables --
  popedInput$design$Gx=t(zeros(0,1))
  
  
  # -- Filname and path of the model file --
  popedInput$ff_file='ff'
  # -- Filname and path of the g parameter file --
  popedInput$fg_file='sfg'
  # -- Filname and path of the error model file --
  popedInput$fError_file='feps'
  # -- Filname and path for the user defined distributions --
  popedInput$strUserDistributionFile=''
  # -- Filname and path for the ED penalty function, empty string means no penalty --
  popedInput$strEDPenaltyFile=''
  # -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
  popedInput$strAutoCorrelationFile=''
  # -- The model title --
  popedInput$modtit='One comp IV bolus with direct Emax effect'
  # -- Use graph output during search --
  popedInput$bShowGraphs=1
  # -- If a log file should be used (0=FALSE, 1=TRUE) --
  popedInput$use_logfile=0
  # -- Filname and path of the output file during search --
  popedInput$output_file='Output$txt'
  # -- Filname suffix of the result function file --
  popedInput$output_function_file='function_output'
  # -- Filename and path for storage of current optimal design --
  popedInput$strIterationFileName=''
  # -- Filename and path for a run file that is used instead of the regular PopED call --
  popedInput$strRunFile=''
  
  ## -- Filname and path of the output file during search --
  popedInput$output_file=paste(match.call()[[1]],'_summary',sep='')
  ## -- Filname suffix of the result function file --
  popedInput$output_function_file=paste(match.call()[[1]],'_output_',sep='')
  ## -- Filename and path for storage of current optimal design --
  popedInput$strIterationFileName=paste(match.call()[[1]],'_current.R',sep='')
  
  
  # -- Method used to calculate M1 (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  popedInput$m1_switch=1
  
  # -- Method used to calculate M2 (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  popedInput$m2_switch=1
  
  # -- Method used to calculate linearization of residual error (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
  popedInput$hle_switch=1
  
  # -- Method used to calculate the gradient of the model (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  popedInput$gradff_switch=1
  
  # -- Method used to calculate the gradient of the parameter vector g (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
  popedInput$gradfg_switch=1
  
  
  # -- Use normal random sampling or Latin Hypercube sampling, 0=Random Sampling, 1=LatinHyperCube --
  popedInput$bLHS=0
  
  # -- Value to interpret as zero in design --
  popedInput$ourzero=1e-05
  # -- Number of iterations in random search between screen output --
  popedInput$rsit_output=5
  # -- Number of iterations in stochastic gradient search between screen output --
  popedInput$sgit_output=1
  
  
  # -- Step length of derivative of linearized model w$r.t. typical values --
  popedInput$hm1=0.0001
  # -- Step length of derivative of model w$r.t. g --
  popedInput$hlf=0.0001
  # -- Step length of derivative of g w$r.t. b --
  popedInput$hlg=0.0001
  # -- Step length of derivative of variance w$r.t. typical values --
  popedInput$hm2=0.0001
  # -- Step length of derivative of OFV w$r.t. time --
  popedInput$hgd=0.0001
  # -- Step length of derivative of model w$r.t. sigma --
  popedInput$hle=0.0001
  # -- The absolute tolerance for the diff equation solver --
  popedInput$AbsTol=1e-05
  # -- The relative tolerance for the diff equation solver --
  popedInput$RelTol=1e-05
  # -- The diff equation solver method, 0=ode45, 1=ode15s --
  popedInput$iDiffSolverMethod=0
  # -- If the differential equation results should be stored in memory or not --
  popedInput$bUseMemorySolver=0
  # -- Fisher Information Matrix type (0=Full FIM, 1=Reduced FIM) --
  popedInput$iFIMCalculationType=0
  # -- Number of Random search iterations --
  popedInput$rsit=300
  # -- Number of Stochastic gradient search iterations --
  popedInput$sgit=150
  # -- Number of Random search iterations when discrete optimization --
  popedInput$intrsit=100
  # -- Number of Stochastic Gradient search iterations when discrete optimization --
  popedInput$intsgit=50
  # -- Iterations until adaptive narrowing in random search --
  popedInput$maxrsnullit=50
  # -- Stoachstic Gradient convergence value, (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
  popedInput$convergence_eps=1e-08
  # -- Random search locality factor for sample times --
  popedInput$rslxt=4
  # -- Random search locality factor for covariates --
  popedInput$rsla=4
  # -- Stochastic Gradient search first step factor for sample times --
  popedInput$cfaxt=0.001
  # -- Stochastic Gradient search first step factor for covariates --
  popedInput$cfaa=0.001
  # -- Use greedy algorithm for group assignment optimization --
  popedInput$bGreedyGroupOpt=0
  # -- Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov --
  popedInput$EACriteria=1
  # -- Exchange Algorithm StepSize --
  popedInput$EAStepSize=0.01
  # -- Exchange Algorithm NumPoints --
  popedInput$EANumPoints=0
  # -- Exchange Algorithm Convergence Limit/Criteria --
  popedInput$EAConvergenceCriteria=1e-20
  # -- Avoid replicate samples when using Exchange Algorithm --
  popedInput$bEANoReplicates=0
  # -- BFGS Minimizer Convergence Criteria Minimum Step --
  popedInput$BFGSConvergenceCriteriaMinStep=1e-08
  # -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
  popedInput$BFGSProjectedGradientTol=0.0001
  # -- BFGS Minimizer Line Search Tolerance f --
  popedInput$BFGSTolerancef=0.001
  # -- BFGS Minimizer Line Search Tolerance g --
  popedInput$BFGSToleranceg=0.9
  # -- BFGS Minimizer Line Search Tolerance x --
  popedInput$BFGSTolerancex=0.1
  # -- Sample size for ED-optimal design distribution --
  popedInput$ED_samp_size=45
  # -- Number of iterations in ED-optimal design to calculate convergence criteria --
  popedInput$ED_diff_it=30
  # -- ED-optimal design convergence criteria in percent --
  popedInput$ED_diff_percent=10
  # -- Number of grid points in the line search --
  popedInput$line_search_it=50
  # -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
  popedInput$iNumSearchIterationsIfNotLineSearch=1
  
  # -- Ds_index, set index to 1 if a parameter is uninteresting, otherwise 0. size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma --
  popedInput$CriterionOptions$ds_index=matrix(c( 0,0,0,0,0,0,0,0,0,0,0,0),nrow=1,byrow=T)
  
  # -- Parallel options for PopED -- --
  
  # -- Compile option for PopED (-1 = No compilation, 0 or 3 = Full compilation, 1 or 4 = Only using MCC (shared lib), 2 or 5 = Only MPI, Option 0,1,2 runs PopED and option 3,4,5 }s after compilation --
  popedInput$parallelSettings$iCompileOption=-1
  # -- Parallel method to use (0 = Matlab PCT, 1 = MPI) --
  popedInput$parallelSettings$iUseParallelMethod=1
  # -- Additional dependencies used in MCC compilation (mat-files), if several space separated --
  popedInput$parallelSettings$strAdditionalMCCCompilerDependencies=''
  # -- Compilation output executable name --
  popedInput$parallelSettings$strExecuteName='calc_fim.exe'
  # -- Number of processes to use when running in parallel (E$g. 3 = 2 workers, 1 job manager) --
  popedInput$parallelSettings$iNumProcesses=2
  # -- Number of design evaluations that should be evaluated in each process before getting new work from job manager --
  popedInput$parallelSettings$iNumChunkDesignEvals=-2
  # -- The prefix of the input mat file to communicate with the excutable --
  popedInput$parallelSettings$strMatFileInputPrefix='parallel_input'
  # -- The prefix of the output mat file to communicate with the excutable --
  popedInput$parallelSettings$strMatFileOutputPrefix='parallel_output'
  # -- Extra options send to e$g. the MPI exectuable or a batch script, see execute_parallel$m for more information and options --
  popedInput$parallelSettings$strExtraRunOptions=''
  # -- Polling time to check if the parallel execution is finished --
  popedInput$parallelSettings$dPollResultTime=1.000000e-01
  # -- The file containing the popedInput structure that should be used to evaluate the designs --
  popedInput$parallelSettings$strFunctionInputName='function_input'
  # -- If the random search is going to be executed in parallel --
  popedInput$parallelSettings$bParallelRS=0
  # -- If the stochastic gradient search is going to be executed in parallel --
  popedInput$parallelSettings$bParallelSG=0
  # -- If the line search is going to be executed in parallel --
  popedInput$parallelSettings$bParallelLS=0
  # -- If the modified exchange algorithm is going to be executed in parallel --
  popedInput$parallelSettings$bParallelMFEA=0
  
  # -- User defined data structure that e$g. could be used to send in data to the model --
  popedInput$user_data=cell(0,0)
  
  
  
  return( popedInput ) 
}
