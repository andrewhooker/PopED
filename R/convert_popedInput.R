#' Convert the poped input file to a database that can be used by PopED.
#' 
#' The function takes the defined input from poped and creates a database for use in other PopED functions.
#' Typically this function will not be accessed by the user.  Instead, use the function \code{\link{create.poped.database}}
#' 
#' @param popedInput An input file to PopED 
#' @param iCompileOption \bold{******START OF PARALLEL OPTIONS**********} Compile options for PopED
#' \itemize{
#' \item -1 = No compilation,
#' \item 0 or 3 = Full compilation,
#' \item 1 or 4 = Only using MCC (shared lib),
#' \item 2 or 5 = Only MPI,
#' \item Option 0,1,2 runs PopED and option 3,4,5 stops after compilation
#' }
#' @return A PopED database
#' @family poped_input

## Function translated automatically using 'matlab.to.r()'
## then manually adjusted
## Author: Andrew Hooker


convert_popedInput <- 
  function(popedInput,
           
          
           ## --------------------------
           ## ---- Model definition
           ## --------------------------
           
           # -- Filname and path of the model file --
           ff_file=poped.choose(popedInput[["ff_file"]],"ff"),
           # -- Filname and path of the g parameter file --
           fg_file=poped.choose(popedInput$fg_file,'sfg'),
           # -- Filname and path of the error model file --
           fError_file=poped.choose(popedInput$fError_file,'feps'),
                      
           ## --------------------------
           ## ---- What to optimize
           ## --------------------------
           
           ## -- Vector of optimization tasks (1=TRUE,0=FALSE)
           ## (Samples per subject, Sampling schedule, Discrete design variable, Continuous design variable, Number of id per group)
           ## -- All elements set to zero => only calculate the FIM with current design --
           optsw=poped.choose(popedInput$optsw,cbind(0,0,0,0,0)),           
           
           ## --------------------------
           ## ---- Initial Design 
           ## --------------------------

           ## -- Matrix defining the initial sampling schedule --
           xt=poped.choose(popedInput$design[["xt"]],stop("'xt' needs to be defined")),
           ## -- Number of groups/individuals --
           m=poped.choose(popedInput[["m"]],size(xt,1)),                     
           ## -- Matrix defining the initial discrete values --
           x=poped.choose(popedInput$design[["x"]],zeros(m,0)),               
           ## -- Number of discrete design variables --
           nx=poped.choose(popedInput$nx,size(x,2)),      
           ## -- Vector defining the initial covariate values --
           a=poped.choose(popedInput$design[["a"]],zeros(m,0)),    
           ## number of continuous design variables that are not time (e.g. continuous covariates)
           na=poped.choose(popedInput$na,size(a,2)),      
           ## -- Vector defining the size of the different groups (num individuals in each group) --
           groupsize=poped.choose(popedInput$design$groupsize,stop("'groupsize' needs to be defined")),      
           ## -- Vector defining the number of samples for each group --
           ni=poped.choose(popedInput$design$ni,matrix(size(xt,2),m,1)),    
           ## -- Vector defining which response a certain sampling time belongs to --
           model_switch=poped.choose(popedInput$design$model_switch,ones(size(xt,1),size(xt,2))),
   
           ## --------------------------
           ## ---- Design space
           ## --------------------------
           
           ## -- Max number of samples per group/individual --
           maxni=poped.choose(popedInput$maxni,size(xt,2)),                     
           ## -- Min number of samples per group/individual --
           minni=poped.choose(popedInput$minni,1),    
           ## -- Vector defining the max size of the different groups (max num individuals in each group) --
           maxgroupsize=poped.choose(popedInput$design$maxgroupsize,groupsize),       
           ## -- Vector defining the min size of the different groups (min num individuals in each group) --
           mingroupsize=poped.choose(popedInput$design$mingroupsize,ones(m,1)),   
           ## -- The total maximal groupsize over all groups--
           maxtotgroupsize=poped.choose(popedInput$design$maxtotgroupsize,sum(groupsize)),   
           ## -- The total minimal groupsize over all groups--
           mintotgroupsize=poped.choose(popedInput$design$mintotgroupsize,sum(mingroupsize)),   
           ## -- Matrix defining the max value for each sample --
           maxxt=poped.choose(popedInput$design$maxxt,xt),   
           ## -- Matrix defining the min value for each sample --
           minxt=poped.choose(popedInput$design$minxt,xt),   
           ## -- Cell defining the discrete variables --
           discrete_x=poped.choose(popedInput$design$discrete_x,cell(m,nx)),     
           ## -- Vector defining the max value for each covariate --
           maxa=poped.choose(popedInput$design$maxa,a),   
           ## -- Vector defining the min value for each covariate --
           mina=poped.choose(popedInput$design$mina,a),   
           ## -- Use grouped time points (1=TRUE, 0=FALSE) --
           bUseGrouped_xt=poped.choose(popedInput$bUseGrouped_xt,FALSE),               
           ## -- Matrix defining the grouping of sample points --
           G_xt=poped.choose(popedInput$design$G,matrix(seq(1,length(xt),len=length(xt)),size(xt,1),size(xt,2),byrow=T)),               
           ## -- Use grouped covariates (1=TRUE, 0=FALSE) --
           bUseGrouped_a=poped.choose(popedInput$bUseGrouped_a,FALSE),               
           ## -- Matrix defining the grouping of covariates --
           G_a=poped.choose(popedInput$design$Ga,matrix(seq(1,length(a),len=length(a)),size(a,1),size(a,2),byrow=T)),               
           ## -- Use grouped discrete design variables (1=TRUE, 0=FALSE) --
           bUseGrouped_x=poped.choose(popedInput$bUseGrouped_x,FALSE),               
           ## -- Matrix defining the grouping of discrete design variables --
           G_x=poped.choose(popedInput$design$Gx,matrix(seq(1,length(x),len=length(x)),size(x,1),size(x,2),byrow=T)),               
           
           
           ## --------------------------
           ## ---- FIM calculation 
           ## --------------------------
           
           ## -- Fisher Information Matrix type
           ## (0=Full FIM,
           ##  1=Reduced FIM,
           ##  2=weighted models,
           ##  3=Loc models,
           ##  4=reduced FIM with derivative of SD of sigma as pfim,
           ##  5=FULL FIM parameterized with A,B,C matrices & derivative of variance,
           ##  6=Calculate one model switch at a time, good for large matrices,
           ##  7=Reduced FIM parameterized with A,B,C matrices & derivative of variance) --
           iFIMCalculationType=poped.choose(popedInput$iFIMCalculationType,1),
           ## -- Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI --
           iApproximationMethod=poped.choose(popedInput$iApproximationMethod,0),
           ## -- Num indivduals in each step of FOCE --
           iFOCENumInd=poped.choose(popedInput$iFOCENumInd,1000),
           ## -- The prior FIM (added to calculated FIM) --
           prior_fim=poped.choose(popedInput$prior_fim,matrix(0,0,1)),
           ## -- Filname and path for the Autocorrelation function, empty string means no autocorrelation --
           strAutoCorrelationFile=poped.choose(popedInput$strAutoCorrelationFile,''),
           
           ## --------------------------
           ## ---- Criterion specification
           ## --------------------------
           
           ## -- D-family design (1) or ED-familty design (0) (with or without parameter uncertainty) --
           d_switch=poped.choose(popedInput$d_switch,1),
           ## -- OFV calculation type for FIM (1=Determinant of FIM,4=log determinant of FIM,6=determinant of interesting part of FIM (Ds)) --
           ofv_calc_type=poped.choose(popedInput$ofv_calc_type,1),
           ## -- Ds_index, set index to 1 if a parameter is uninteresting, otherwise 0.
           ## size=(1,num unfixed parameters). First unfixed bpop, then unfixed d, then unfixed docc and last unfixed sigma --
           ## default is the fixed effects being important
           ds_index=popedInput$CriterionOptions$ds_index,   
           ## -- Penalty function, empty string means no penalty.  User defined criterion --
           strEDPenaltyFile=poped.choose(popedInput$strEDPenaltyFile,''),
           
           
           
           ## --------------------------
           ## ---- E-family Criterion options
           ## --------------------------
           ## -- ED Integral Calculation, 0=Monte-Carlo-Integration, 1=Laplace Approximation, 2=BFGS Laplace Approximation  -- --
           iEDCalculationType=poped.choose(popedInput$iEDCalculationType,0),     
           ## -- Sample size for E-family sampling --
           ED_samp_size=poped.choose(popedInput$ED_samp_size,45),     
           ## -- How to sample from distributions in E-family calculations. 0=Random Sampling, 1=LatinHyperCube --
           bLHS=poped.choose(popedInput$bLHS,1),     
           ## -- Filname and path for user defined distributions for E-family designs --
           strUserDistributionFile=poped.choose(popedInput$strUserDistributionFile,''), 
           
           ## --------------------------
           ## ---- Model parameters 
           ## --------------------------
           
           ## -- Number of typical values --
           nbpop=popedInput$nbpop,
           ## -- Number of IIV parameters --
           NumRanEff=popedInput$nb,
           ## -- Number of IOV variance parameters --
           NumDocc=popedInput$ndocc,
           ## -- Number of occassions --
           NumOcc= popedInput$NumOcc,
           ## -- The length of the g parameter vector --
           ng=popedInput$ng,
           
           ## -- Matrix defining the fixed effects, per row (row number = parameter_number),
           ## the type of the distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform,
           ## 3 = User Defined Distribution, 4 = lognormal and 5 = truncated normal).
           ## The second column defines the mean.
           ## The third column defines the variance of the distribution.
           # can also just supply the parameter values as a c()
           bpop=poped.choose(popedInput$design$bpop,stop('bpop must be defined')),
           ## -- Matrix defining the diagnonals of the IIV (same logic as for the fixed efects) --
           # can also just supply the parameter values as a c()
           d=poped.choose(popedInput$design$d,stop('d must be defined')),
           ## -- Matrix defining the covariances of the IIV variances --
           # set to zero if not defined
           covd=popedInput$design$covd,
           ## -- Matrix defining the variances of the residual variability terms --
           ## REQUIRED! No defaults given.
           sigma=popedInput$design$sigma,
           ## -- Matrix defining the IOV, the IOV variances and the IOV distribution --
           docc=poped.choose(popedInput$design$docc,matrix(0,0,3)),
           ## -- Matrix defining the covariance of the IOV --
           covdocc=poped.choose(popedInput$design$covdocc,zeros(1,length(docc[,2,drop=F])*(length(docc[,2,drop=F])-1)/2)),
           
           ## --------------------------
           ## ---- Model parameters fixed or not
           ## --------------------------
           ## -- Vector defining if a typical value is fixed or not (1=not fixed, 0=fixed) --
           notfixed_bpop=popedInput$notfixed_bpop,
           ## -- Vector defining if a IIV is fixed or not (1=not fixed, 0=fixed) --
           notfixed_d=popedInput$notfixed_d,
           ## -- Vector defining if a covariance IIV is fixed or not (1=not fixed, 0=fixed) --
           notfixed_covd=popedInput$notfixed_covd,           
           ## -- Vector defining if an IOV variance is fixed or not (1=not fixed, 0=fixed) --
           notfixed_docc=popedInput$notfixed_docc,
           ## -- Vector row major order for lower triangular matrix defining if a covariance IOV is fixed or not (1=not fixed, 0=fixed) --
           notfixed_covdocc=poped.choose(popedInput$notfixed_covdocc,zeros(1,length(covdocc))),
           ## -- Vector defining if a residual error parameter is fixed or not (1=not fixed, 0=fixed) --
           notfixed_sigma=poped.choose(popedInput$notfixed_sigma,t(rep(1,length(diag(sigma))))),   
           ## -- Vector defining if a covariance residual error parameter is fixed or not (1=not fixed, 0=fixed) --
           ## default is fixed
           notfixed_covsigma=poped.choose(popedInput$notfixed_covsigma,zeros(1,length(notfixed_sigma)*(length(notfixed_sigma)-1)/2)),  
           
           
           ## --------------------------
           ## ---- Optimization algorithm choices
           ## --------------------------
           
           ## -- Use random search (1=TRUE, 0=FALSE) --
           bUseRandomSearch=poped.choose(popedInput$bUseRandomSearch,TRUE),           
           ## -- Use Stochastic Gradient search (1=TRUE, 0=FALSE) --
           bUseStochasticGradient=poped.choose(popedInput$bUseStochasticGradient,TRUE),
           ## -- Use Line search (1=TRUE, 0=FALSE) --
           bUseLineSearch=poped.choose(popedInput$bUseLineSearch,TRUE),
           ## -- Use Exchange algorithm (1=TRUE, 0=FALSE) --
           bUseExchangeAlgorithm=poped.choose(popedInput$bUseExchangeAlgorithm,FALSE),       
           ## -- Use BFGS Minimizer (1=TRUE, 0=FALSE) --
           bUseBFGSMinimizer=poped.choose(popedInput$bUseBFGSMinimizer,FALSE),
           ## -- Exchange Algorithm Criteria, 1 = Modified, 2 = Fedorov --
           EACriteria=poped.choose(popedInput$EACriteria,1),
           ## -- Filename and path for a run file that is used instead of the regular PopED call --
           strRunFile=poped.choose(popedInput$strRunFile,''),
           
           ## --------------------------
           ## ---- Labeling and file names
           ## --------------------------
           
           ## -- The current PopED version --
           strPopEDVersion=poped.choose(popedInput$strPopEDVersion, packageVersion("PopED")),  
           ## -- The model title --
           modtit=poped.choose(popedInput$modtit,'PopED model'),
           ## -- Filname and path of the output file during search --
           output_file=poped.choose(popedInput$output_file,paste("PopED_output",'_summary',sep='')),
           ## -- Filname suffix of the result function file --
           output_function_file=poped.choose(popedInput$output_function_file,paste("PopED",'_output_',sep='')),
           ## -- Filename and path for storage of current optimal design --
           strIterationFileName=poped.choose(popedInput$strIterationFileName,paste("PopED",'_current.R',sep='')),
           
           
           ## --------------------------
           ## ---- Misc options
           ## --------------------------
           ## -- User defined data structure that, for example could be used to send in data to the model --
           user_data=poped.choose(popedInput$user_data,cell(0,0)),
           ## -- Value to interpret as zero in design --
           ourzero=poped.choose(popedInput$ourzero,1e-5),                    
           ## -- The seed number used for optimization and sampling -- integer or -1 which creates a random seed
           dSeed=poped.choose(popedInput$dSeed,-1),
           ## -- Vector for line search on continuous design variables (1=TRUE,0=FALSE) --
           line_opta=poped.choose(popedInput$line_opta,TRUE),
           ## -- Vector for line search on discrete design variables (1=TRUE,0=FALSE) --
           line_optx=poped.choose(popedInput$line_optx,TRUE), #matrix(0,0,1)           
           ## -- Use graph output during search --
           bShowGraphs=poped.choose(popedInput$bShowGraphs,FALSE),
           ## -- If a log file should be used (0=FALSE, 1=TRUE) --
           use_logfile=poped.choose(popedInput$use_logfile,FALSE),          
           ## -- Method used to calculate M1
           ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           m1_switch=poped.choose(popedInput$m1_switch,1),
           ## -- Method used to calculate M2
           ## (0=Central difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           m2_switch=poped.choose(popedInput$m2_switch,1),
           ## -- Method used to calculate linearization of residual error
           ## (0=Complex difference, 1=Central difference, 30=Automatic differentiation) --
           hle_switch=poped.choose(popedInput$hle_switch,1),
           ## -- Method used to calculate the gradient of the model
           ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           gradff_switch=poped.choose(popedInput$gradff_switch,1),
           ## -- Method used to calculate the gradient of the parameter vector g
           ## (0=Complex difference, 1=Central difference, 20=Analytic derivative, 30=Automatic differentiation) --
           gradfg_switch=poped.choose(popedInput$gradfg_switch,1),
           ## -- Number of iterations in random search between screen output --
           rsit_output=poped.choose(popedInput$rsit_output,5),          
           ## -- Number of iterations in stochastic gradient search between screen output --
           sgit_output=poped.choose(popedInput$sgit_output,1),
           ## -- Step length of derivative of linearized model w.r.t. typical values --
           hm1=poped.choose(popedInput$hm1,0.00001),
           ## -- Step length of derivative of model w.r.t. g --
           hlf=poped.choose(popedInput$hlf,0.00001),
           ## -- Step length of derivative of g w.r.t. b --
           hlg=poped.choose(popedInput$hlg,0.00001),
           ## -- Step length of derivative of variance w.r.t. typical values --
           hm2=poped.choose(popedInput$hm2,0.00001),
           ## -- Step length of derivative of OFV w.r.t. time --
           hgd=poped.choose(popedInput$hgd,0.00001),
           ## -- Step length of derivative of model w.r.t. sigma --
           hle=poped.choose(popedInput$hle,0.00001),
           ## -- The absolute tolerance for the diff equation solver --
           AbsTol=poped.choose(popedInput$AbsTol,0.00001),
           ## -- The relative tolerance for the diff equation solver --
           RelTol=poped.choose(popedInput$RelTol,0.00001),
           ## -- The diff equation solver method, 0, no other option --
           iDiffSolverMethod=poped.choose(popedInput$iDiffSolverMethod,0),
           ## -- If the differential equation results should be stored in memory (1) or not (0) --
           bUseMemorySolver=poped.choose(popedInput$bUseMemorySolver,FALSE),
           ## -- Number of Random search iterations --
           rsit=poped.choose(popedInput$rsit,300),
           ## -- Number of Stochastic gradient search iterations --
           sgit=poped.choose(popedInput$sgit,150),
           ## -- Number of Random search iterations with discrete optimization --
           intrsit=poped.choose(popedInput$intrsit,250),
           ## -- Number of Stochastic Gradient search iterations with discrete optimization --
           intsgit=poped.choose(popedInput$intsgit,50),
           ## -- Iterations until adaptive narrowing in random search --
           maxrsnullit=poped.choose(popedInput$maxrsnullit,50),
           ## -- Stoachstic Gradient convergence value,
           ## (difference in OFV for D-optimal, difference in gradient for ED-optimal) --
           convergence_eps=poped.choose(popedInput$convergence_eps,1e-08),
           ## -- Random search locality factor for sample times --
           rslxt=poped.choose(popedInput$rslxt,10),
           ## -- Random search locality factor for covariates --
           rsla=poped.choose(popedInput$rsla,10),
           ## -- Stochastic Gradient search first step factor for sample times --
           cfaxt=poped.choose(popedInput$cfaxt,0.001),
           ## -- Stochastic Gradient search first step factor for covariates --
           cfaa=poped.choose(popedInput$cfaa,0.001),
           ## -- Use greedy algorithm for group assignment optimization --
           bGreedyGroupOpt=poped.choose(popedInput$bGreedyGroupOpt,FALSE),
           ## -- Exchange Algorithm StepSize --
           EAStepSize=poped.choose(popedInput$EAStepSize,0.01),
           ## -- Exchange Algorithm NumPoints --
           EANumPoints=poped.choose(popedInput$EANumPoints,FALSE),
           ## -- Exchange Algorithm Convergence Limit/Criteria --
           EAConvergenceCriteria=poped.choose(popedInput$EAConvergenceCriteria,1e-20),
           ## -- Avoid replicate samples when using Exchange Algorithm --
           bEANoReplicates=poped.choose(popedInput$bEANoReplicates,FALSE),
           ## -- BFGS Minimizer Convergence Criteria Minimum Step --
           BFGSConvergenceCriteriaMinStep=poped.choose(popedInput$BFGSConvergenceCriteriaMinStep,1e-08),
           ## -- BFGS Minimizer Convergence Criteria Normalized Projected Gradient Tolerance --
           BFGSProjectedGradientTol=poped.choose(popedInput$BFGSProjectedGradientTol,0.0001),
           ## -- BFGS Minimizer Line Search Tolerance f --
           BFGSTolerancef=poped.choose(popedInput$BFGSTolerancef,0.001),
           ## -- BFGS Minimizer Line Search Tolerance g --
           BFGSToleranceg=poped.choose(popedInput$BFGSToleranceg,0.9),
           ## -- BFGS Minimizer Line Search Tolerance x --
           BFGSTolerancex=poped.choose(popedInput$BFGSTolerancex,0.1),
           ## -- Number of iterations in ED-optimal design to calculate convergence criteria --
           ED_diff_it=poped.choose(popedInput$ED_diff_it,30),
           ## -- ED-optimal design convergence criteria in percent --
           ED_diff_percent=poped.choose(popedInput$ED_diff_percent,10),
           ## -- Number of grid points in the line search --
           line_search_it=poped.choose(popedInput$line_search_it,50),
           ## -- Number of iterations of full Random search and full Stochastic Gradient if line search is not used --
           iNumSearchIterationsIfNotLineSearch=poped.choose(popedInput$iNumSearchIterationsIfNotLineSearch,1),
           
           ## --------------------------
           ## -- Parallel options for PopED -- --
           ## --------------------------
           #     ## -- Compile option for PopED
           #     ## -1 = No compilation,
           #     ## 0 or 3 = Full compilation,
           #     ## 1 or 4 = Only using MCC (shared lib),
           #     ## 2 or 5 = Only MPI,
           #     ## Option 0,1,2 runs PopED and option 3,4,5 stops after compilation --
           iCompileOption=poped.choose(popedInput$parallelSettings$iCompileOption,-1),
           ## -- Parallel method to use (0 = Matlab PCT, 1 = MPI) --
           iUseParallelMethod=poped.choose(popedInput$parallelSettings$iUseParallelMethod,1),
           ## -- Additional dependencies used in MCC compilation (mat-files), if several space separated --
           strAdditionalMCCCompilerDependencies=poped.choose(popedInput$parallelSettings$strAdditionalMCCCompilerDependencies,''),
           ## -- Compilation output executable name --
           strExecuteName=poped.choose(popedInput$parallelSettings$strExecuteName,'calc_fim.exe'),
           ## -- Number of processes to use when running in parallel (e.g. 3 = 2 workers, 1 job manager) --
           iNumProcesses=poped.choose(popedInput$parallelSettings$iNumProcesses,2),
           ## -- Number of design evaluations that should be evaluated in each process before getting new work from job manager --
           iNumChunkDesignEvals=poped.choose(popedInput$parallelSettings$iNumChunkDesignEvals,-2),
           ## -- The prefix of the input mat file to communicate with the excutable --
           strMatFileInputPrefix=poped.choose(popedInput$parallelSettings$strMatFileInputPrefix,'parallel_input'),
           ## -- The prefix of the output mat file to communicate with the excutable --
           strMatFileOutputPrefix=poped.choose(popedInput$parallelSettings$strMatFileOutputPrefix,'parallel_output'),
           ## -- Extra options send to e$g. the MPI exectuable or a batch script, see execute_parallel$m for more information and options --
           strExtraRunOptions=poped.choose(popedInput$parallelSettings$strExtraRunOptions,''),
           ## -- Polling time to check if the parallel execution is finished --
           dPollResultTime=poped.choose(popedInput$parallelSettings$dPollResultTime,0.1),
           ## -- The file containing the popedInput structure that should be used to evaluate the designs --
           strFunctionInputName=poped.choose(popedInput$parallelSettings$strFunctionInputName,'function_input'),
           ## -- If the random search is going to be executed in parallel --
           bParallelRS=poped.choose(popedInput$parallelSettings$bParallelRS,FALSE),
           ## -- If the stochastic gradient search is going to be executed in parallel --
           bParallelSG=poped.choose(popedInput$parallelSettings$bParallelSG,FALSE),
           ## -- If the modified exchange algorithm is going to be executed in parallel --
           bParallelMFEA=poped.choose(popedInput$parallelSettings$bParallelMFEA,FALSE),
           ## -- If the line search is going to be executed in parallel --
           bParallelLS=poped.choose(popedInput$parallelSettings$bParallelLS,FALSE)
           
  ){
    globalStructure <- list()
    
    #     # update popedInput with options supplied in function
    #     called_args <- match.call()
    #     default_args <- formals()
    #     for(i in names(called_args)[-1]){
    #       if(length(grep("^popedInput$",capture.output(default_args[[i]])))==1) {
    #         eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
    #       }
    #     }
    
    #modifyList(settings, list(…)$settings)
    
    ## compare to a default input function.
    #   ## -- Filname and path of the model file --
    #   popedInput$ff_file='ff'
    #   ## -- Filname and path of the parameter file --
    #   popedInput$fg_file='sfg'
    #   ## -- Filname and path of the residual error model file --
    #   popedInput$fError_file='feps.add.prop'
    #   ## -- The model title --
    #   popedInput$modtit='Sigmoidal Emax model'
    
    Engine = list(Type=1,Version=version$version.string)
    
    globalStructure$strPopEDVersion = strPopEDVersion
    globalStructure$user_distribution_pointer=''
    
    
    
    globalStructure$nx = nx
    globalStructure$na = na
    
    globalStructure$bLHS = bLHS
    
    globalStructure$discrete_x = popedInput$design$discrete_x
    
    globalStructure$m = m
    globalStructure$maxni=maxni
    globalStructure$minni=minni
    globalStructure$maxm=globalStructure$m # max number of individuals or groups    
    globalStructure$maxmaxni=globalStructure$maxni # max number of samples per individual or group
    
    globalStructure$bUseGrouped_xt = bUseGrouped_xt
    globalStructure$bUseGrouped_a  = bUseGrouped_a
    globalStructure$bUseGrouped_x  = bUseGrouped_x
    
    globalStructure$d_switch = d_switch
    globalStructure$iApproximationMethod = iApproximationMethod
    globalStructure$iFOCENumInd = iFOCENumInd
    globalStructure$bUseRandomSearch = bUseRandomSearch
    globalStructure$bUseStochasticGradient = bUseStochasticGradient
    globalStructure$bUseLineSearch = bUseLineSearch
    globalStructure$bUseExchangeAlgorithm = bUseExchangeAlgorithm
    
    globalStructure$bUseBFGSMinimizer=bUseBFGSMinimizer
    
    
    globalStructure$iEDCalculationType=iEDCalculationType
    
    
    
    globalStructure$BFGSConvergenceCriteriaMinStep=BFGSConvergenceCriteriaMinStep
    
    
    globalStructure$BFGSProjectedGradientTol=BFGSProjectedGradientTol
    
    
    
    globalStructure$BFGSTolerancef=BFGSTolerancef
    
    
    globalStructure$BFGSToleranceg=BFGSToleranceg
    
    
    globalStructure$BFGSTolerancex=BFGSTolerancex
    
    
    
    
    
    globalStructure$covdocc=covdocc
    
    
    globalStructure$notfixed_covdocc=notfixed_covdocc
    
    
    
    
    globalStructure$notfixed_covsigma=notfixed_covsigma
    
    
    
    
    
    globalStructure$parallelSettings$iCompileOption = iCompileOption
    globalStructure$parallelSettings$strAdditionalMCCCompilerDependencies = strAdditionalMCCCompilerDependencies
    globalStructure$parallelSettings$iUseParallelMethod = iUseParallelMethod
    globalStructure$parallelSettings$strExecuteName = strExecuteName
    globalStructure$parallelSettings$iNumProcesses = iNumProcesses
    globalStructure$parallelSettings$iNumChunkDesignEvals = iNumChunkDesignEvals
    globalStructure$parallelSettings$strMatFileInputPrefix = strMatFileInputPrefix
    globalStructure$parallelSettings$strMatFileOutputPrefix = strMatFileOutputPrefix
    globalStructure$parallelSettings$strExtraRunOptions = strExtraRunOptions
    globalStructure$parallelSettings$dPollResultTime = dPollResultTime
    globalStructure$parallelSettings$strFunctionInputName = strFunctionInputName
    globalStructure$parallelSettings$bParallelRS = bParallelRS  
    globalStructure$parallelSettings$bParallelSG = bParallelSG
    globalStructure$parallelSettings$bParallelLS = bParallelLS
    globalStructure$parallelSettings$bParallelMFEA = bParallelMFEA
    
    
    globalStructure$hm1=hm1
    globalStructure$hlf=hlf
    globalStructure$hlg=hlg
    globalStructure$hm2=hm2
    globalStructure$hgd=hgd
    globalStructure$hle=hle
    globalStructure$numericalSettings$AbsTol = AbsTol
    globalStructure$numericalSettings$RelTol = RelTol
    globalStructure$numericalSettings$iDiffSolverMethod = iDiffSolverMethod
    
    #Temp thing for memory solvers
    globalStructure$numericalSettings$bUseMemorySolver = bUseMemorySolver
    globalStructure$numericalSettings$solved_solutions = cell(0,0)
    globalStructure$numericalSettings$maxtime = max(max(maxxt))+hgd
    
    globalStructure$iFIMCalculationType = iFIMCalculationType
    globalStructure$rsit=rsit
    globalStructure$sgit=sgit
    globalStructure$intrsit=intrsit
    globalStructure$intsgit=intsgit
    globalStructure$maxrsnullit=maxrsnullit
    globalStructure$convergence_eps=convergence_eps
    globalStructure$rslxt=rslxt
    globalStructure$rsla=rsla
    globalStructure$cfaxt=cfaxt
    globalStructure$cfaa=cfaa
    
    globalStructure$EACriteria = EACriteria
    globalStructure$EAStepSize = EAStepSize
    globalStructure$EANumPoints = EANumPoints
    globalStructure$EAConvergenceCriteria = EAConvergenceCriteria
    
    
    globalStructure$ED_samp_size=ED_samp_size
    globalStructure$ED_diff_it = ED_diff_it
    globalStructure$ED_diff_percent = ED_diff_percent
    globalStructure$ls_step_size=line_search_it
    
    globalStructure$ofv_calc_type = ofv_calc_type
    
    globalStructure$iNumSearchIterationsIfNotLineSearch = iNumSearchIterationsIfNotLineSearch
    
    globalStructure$ourzero=ourzero
    globalStructure$rsit_output=rsit_output
    globalStructure$sgit_output=sgit_output
    
    
    if(exists(fg_file)){
      globalStructure$fg_pointer = fg_file
    } else {
      source(fg_file)
      returnArgs <-  fileparts(fg_file) 
      strfgModelFilePath <- returnArgs[[1]]
      strfgModelFilename  <- returnArgs[[2]]
      ## if (~strcmp(strfgModelFilePath,''))
      ##    cd(strfgModelFilePath);
      ## end
      globalStructure$fg_pointer = strfgModelFilename
    }
    
    
    if((!as.character(strUserDistributionFile)=='')){
      if(exists(strUserDistributionFile)){
        globalStructure$user_distribution_pointer = strUserDistributionFile
      } else {
        source(strUserDistributionFile) 
        returnArgs <-  fileparts(strUserDistributionFile) 
        strUserDistFilePath <- returnArgs[[1]]
        strUserDistFilename  <- returnArgs[[2]]
        ##     if (~strcmp(strUserDistFilePath,''))
        ##        cd(strUserDistFilePath);
        ##     end
        globalStructure$user_distribution_pointer = strUserDistFilename
      }
    }
    
    globalStructure$ed_penalty_pointer=zeros(1,0)
    
    if((!as.character(strEDPenaltyFile)=='')){
      if(exists(strEDPenaltyFile)){
        globalStructure$ed_penalty_pointer = strEDPenaltyFile
      } else {
        source(popedInput$strEDPenaltyFile) 
        returnArgs <-  fileparts(popedInput$strEDPenaltyFile) 
        strEDPenaltyFilePath <- returnArgs[[1]]
        strEDPenaltyFilename  <- returnArgs[[2]]
        ##     if (~strcmp(strEDPenaltyFilePath,''))
        ##        cd(strEDPenaltyFilePath);
        ##     end
        globalStructure$ed_penalty_pointer = strEDPenaltyFilename
      }
    }    
    globalStructure$auto_pointer=zeros(1,0)
    
    if((!as.character(strAutoCorrelationFile)=='')){
      if(exists(strAutoCorrelationFile)){
        globalStructure$auto_pointer = strAutoCorrelationFile
      } else {
        source(popedInput$strAutoCorrelationFile) 
        returnArgs <-  fileparts(popedInput$strAutoCorrelationFile) 
        strAutoCorrelationFilePath <- returnArgs[[1]]
        strAutoCorrelationFilename  <- returnArgs[[2]]
        ##     if (~strcmp(strAutoCorrelationFilePath,''))
        ##        cd(strAutoCorrelationFilePath);
        ##     end
        globalStructure$auto_pointer = strAutoCorrelationFilename
      }
    }
  
    if(exists(ff_file)){
      globalStructure$ff_pointer = ff_file 
    } else {
      source(ff_file)
      returnArgs <-  fileparts(ff_file) 
      strffModelFilePath <- returnArgs[[1]]
      strffModelFilename  <- returnArgs[[2]]
      ## if (~strcmp(strffModelFilePath,''))
      ##    cd(strffModelFilePath);
      ## end
      globalStructure$ff_pointer = strffModelFilename
    }
    
    #Check if there is any sub models defined
    if((isfield(popedInput,'SubModels'))){
      i=1
      while(isfield(popedInput$SubModels,sprintf('ff_file%d',i))){
        source(eval(sprintf('popedInput$SubModels$ff_file%d',i))) ##ok<NASGU> 
        returnArgs <-  fileparts(eval(sprintf('popedInput$SubModels$ff_file%d',i))) ##ok<NASGU> 
        strffModelFilePath <- returnArgs[[1]]
        strffModelFilename  <- returnArgs[[2]]
        ##         if (~strcmp(strffModelFilePath,''))
        ##             cd(strffModelFilePath);
        ##         end
        globalStructure$subffPointers[paste('ff_pointer',i,sep='')] = strffModelFilename
        i=i+1
      }
    }
    
    if(exists(fError_file)){
      globalStructure$ferror_pointer = fError_file
    } else {
      source(fError_file)
      returnArgs <-  fileparts(fError_file) 
      strErrorModelFilePath <- returnArgs[[1]]
      strErrorModelFilename  <- returnArgs[[2]]
      ## if (~strcmp(strErrorModelFilePath,''))
      ##    cd(strErrorModelFilePath);
      ## end
      globalStructure$ferror_pointer = strErrorModelFilename
    }
    
    if((Engine$Type==1) ){#Matlab
      ##  %Set the model file string path
      globalStructure$model_file = ff_file
      ##   model_file = eval('functions(globalStructure.ff_pointer)');
      ##   if (~strcmp(model_file.file,''))
      ##       globalStructure.model_file = eval('char(model_file.file)');
      ##   else
      ##       globalStructure.model_file = eval('char(model_file.function)');
      ##   end
    } else {   #FreeMat
      globalStructure$model_file = ff_file
    }
    
    #==================================
    # Initialize the randomization    
    #==================================
    if((dSeed == -1)){
      globalStructure$dSeed = as.integer(Sys.time())
    } else {
      globalStructure$dSeed = dSeed
    }
    if((Engine$Type==1) ){#Matlab
      set.seed(globalStructure$dSeed)
      ##    randn('state', globalStructure.dSeed);
    } else {                #FreeMat
      ##seed(round(globalStructure$dSeed),round(globalStructure$dSeed))
      ##tmp = round(rand(1)*2^30)
      ##seed(globalStructure$dSeed,round(tmp))
    }
    
    globalStructure$nbpop = poped.choose(nbpop,find.largest.index(globalStructure$fg_pointer,"bpop"))
    globalStructure$NumRanEff = poped.choose(NumRanEff,find.largest.index(globalStructure$fg_pointer,"b"))
    globalStructure$NumDocc = poped.choose(NumDocc,find.largest.index(globalStructure$fg_pointer,"bocc",mat=T,mat.row=T))
    globalStructure$NumOcc = poped.choose(NumOcc,find.largest.index(globalStructure$fg_pointer,"bocc",mat=T,mat.row=F))
    globalStructure$ng = poped.choose(ng,length(sfg(0,0,0,0,0)))    
    
    globalStructure$notfixed_docc = poped.choose(notfixed_docc,matrix(1,nrow=1,ncol=globalStructure$NumDocc))
    globalStructure$notfixed_d = poped.choose(notfixed_d,matrix(1,nrow=1,ncol=globalStructure$NumRanEff))
    globalStructure$notfixed_bpop = poped.choose(notfixed_bpop,matrix(1,nrow=1,ncol=globalStructure$nbpop))
    
    if(size(d,1)==1 && size(d,2)==globalStructure$NumRanEff){ # we have just the parameter values not the uncertainty
      d_descr <- zeros(globalStructure$NumRanEff,3)
      d_descr[,2] <- d
      d_descr[,1] <- 0 # point values
      d_descr[,3] <- 0 # variance
      d <- d_descr
    }
    
    if(size(bpop,1)==1 && size(bpop,2)==globalStructure$nbpop){ # we have just the parameter values not the uncertainty
      bpop_descr <- zeros(globalStructure$nbpop,3)
      bpop_descr[,2] <- bpop
      bpop_descr[,1] <- 0 # point values
      bpop_descr[,3] <- 0 # variance
      bpop <- bpop_descr
    }    
    
    covd = poped.choose(covd,zeros(1,globalStructure$NumRanEff)*(globalStructure$NumRanEff-1)/2)
    globalStructure$covd = covd
    
    tmp <- ones(1,length(covd))
    tmp[covd==0] <- 0
    globalStructure$notfixed_covd=poped.choose(notfixed_covd,tmp)
    
    
    #==================================
    # Sample the individual eta's for FOCE and FOCEI
    #==================================
    if((globalStructure$iApproximationMethod!=0 && globalStructure$iApproximationMethod!=3)){
      
      iMaxCorrIndNeeded = 100
      
      bzeros=zeros(globalStructure$NumRanEff,1)
      bones = matrix(1,globalStructure$NumRanEff,1)
      bocczeros=zeros(globalStructure$NumDocc,1)
      boccones=matrix(1,globalStructure$NumDocc,1)
      
      globalStructure$b_global=zeros(globalStructure$NumRanEff,max(globalStructure$iFOCENumInd,iMaxCorrIndNeeded))
      
      fulld = getfulld(d[,2],globalStructure$covd)
      fulldocc = getfulld(docc[,2,drop=F],globalStructure$covdocc)
      
      globalStructure$bocc_global = cell(globalStructure$iFOCENumInd,1)
      
      
      if((globalStructure$d_switch)){
        globalStructure$b_global = 
          pargen(
            matrix(
              c(
                bones,
                bzeros,
                d[,2]),nrow=1,byrow=T),
            0,
            max(globalStructure$iFOCENumInd,iMaxCorrIndNeeded),
            globalStructure$bLHS,
            t(zeros(1,0)),
            globalStructure)
        
        for(i in 1:globalStructure$iFOCENumInd){
          globalStructure$bocc_global[[i]]=zeros(size(docc,1),globalStructure$NumOcc)
          globalStructure$bocc_global[[i]]=pargen(matrix(c(boccones,bocczeros,docc[,2,drop=F]),nrow=1,byrow=T),0,globalStructure$NumOcc,globalStructure$bLHS,t(zeros(1,0)),globalStructure)
        }
      } else {
        d_dist=pargen(d,globalStructure$user_distribution_pointer,max(globalStructure$iFOCENumInd,iMaxCorrIndNeeded),globalStructure$bLHS,zeros(1,0),globalStructure)
        docc_dist=pargen(docc,globalStructure$user_distribution_pointer,globalStructure$iFOCENumInd,globalStructure$bLHS,zeros(1,0),globalStructure)
        
        if((!isempty(d_dist))){
          for(i in 1:max(globalStructure$iFOCENumInd,iMaxCorrIndNeeded)){
            tmp_d_dist = 
              matrix(
                c(matrix(1,length(d[,2]),1),
                  t(zeros(length(d[,2]),1)),d_dist(i,)),nrow=1,byrow=T)
            globalStructure$b_global[,i] = pargen(tmp_d_dist,globalStructure$user_distribution_pointer,1,globalStructure$bLHS,i,globalStructure)
          }
        }
        
        if((!isempty(docc_dist))){
          for(i in 1:globalStructure$iFOCENumInd){
            tmp_docc_dist = 
              matrix(
                c(
                  matrix(1,length(docc[,2,drop=F]),1),
                  t(zeros(length(docc[,2,drop=F]),1)),docc_dist(i,)),nrow=1,byrow=T)
            globalStructure$bocc_global[[i]]=t(pargen(tmp_docc_dist,globalStructure$user_distribution_pointer,globalStructure$NumOcc,globalStructure$bLHS,i,globalStructure))
          }
        }
      }
      
      #We can use the same correlateSamples for ED because the correlation
      #matrix for the different Omegas is the same
      globalStructure$b_global = t(correlateSamples(globalStructure$b_global,fulld))
      #Add correlation for bocc
      for(i in 1:length(globalStructure$bocc_global)){
        globalStructure$bocc_global[[i]]= t(correlateSamples(globalStructure$bocc_global[[i]],fulldocc))
      }
      
    } else {
      globalStructure$b_global=zeros(globalStructure$NumRanEff,1)
      globalStructure$bocc_global = cell(1,1)
      globalStructure$bocc_global[[1]]=zeros(size(docc,1),globalStructure$NumOcc)
      globalStructure$iFOCENumInd = 1
    }
    
    globalStructure$modtit=modtit
    globalStructure$exptit= sprintf('%s_exp$mat',modtit) #experiment settings title
    globalStructure$opttit= sprintf('%s_opt$mat',modtit) #optimization settings title
    
    globalStructure$bShowGraphs=bShowGraphs
    globalStructure$use_logfile=use_logfile
    globalStructure$output_file=output_file
    globalStructure$output_function_file=output_function_file
    
    globalStructure$optsw=optsw
    
    globalStructure$line_opta=line_opta
    globalStructure$line_optx=line_optx
    globalStructure$design = popedInput$design
    globalStructure$design$bpop = bpop
    globalStructure$design$d = d
    globalStructure$design$covd = covd
    globalStructure$design$sigma = sigma
    globalStructure$design$docc = docc
    globalStructure$design$covdocc = covdocc
    globalStructure$design$x = x
    globalStructure$design$xt = xt
    globalStructure$design$ni = ni
    
    globalStructure$design$G = G_xt
    globalStructure$design$Ga = G_a
    globalStructure$design$Gx = G_x
    
    globalStructure$design$a = a
    globalStructure$design$groupsize = groupsize
    globalStructure$design$maxgroupsize = maxgroupsize
    globalStructure$design$mingroupsize = mingroupsize
    globalStructure$design$maxtotgroupsize = maxtotgroupsize
    globalStructure$design$mintotgroupsize = mintotgroupsize
    if(length(maxxt)==1) maxxt=ones(size(xt,1),size(xt,2))*maxxt
    globalStructure$design$maxxt = maxxt
    if(length(minxt)==1) minxt=ones(size(xt,1),size(xt,2))*minxt
    globalStructure$design$minxt = minxt
    globalStructure$design$discrete_x = discrete_x
    globalStructure$design$maxa = maxa
    globalStructure$design$mina = mina
    globalStructure$design$model_switch = model_switch
    
    
    globalStructure$m1_switch = m1_switch
    globalStructure$m2_switch = m2_switch
    globalStructure$hle_switch = hle_switch
    globalStructure$gradff_switch=gradff_switch
    globalStructure$gradfg_switch = gradfg_switch
    
    globalStructure$prior_fim = prior_fim
    
    globalStructure$notfixed_sigma =  notfixed_sigma
    globalStructure$sigma = sigma
    globalStructure$docc = docc
    
    
    globalStructure$ds_index = ds_index
    
    globalStructure$strIterationFileName = strIterationFileName
    
    globalStructure$user_data = user_data
    
    globalStructure$bUseSecondOrder = FALSE
    globalStructure$bCalculateEBE = FALSE
    
    globalStructure$bGreedyGroupOpt = bGreedyGroupOpt
    
    
    globalStructure$bEANoReplicates = bEANoReplicates
    
    
    if(strRunFile==''){
      globalStructure$run_file_pointer=zeros(1,0)
    } else {
      
      if(exists(strRunFile)){
        globalStructure$strRunFile = strRunFile
      } else {
        source(popedInput$strRunFile)
        returnArgs <-  fileparts(popedInput$strRunFile) 
        strRunFilePath <- returnArgs[[1]]
        strRunFilename  <- returnArgs[[2]]
        ## if (~strcmp(strErrorModelFilePath,''))
        ##    cd(strErrorModelFilePath);
        ## end
        globalStructure$run_file_pointer = strRunFileFilename
      }
    }
    
    
    
    
    return( globalStructure) 
  }

poped.choose <- function(arg1,arg2){
  #ifelse(!is.null(arg1), arg1, arg2)
  if(!is.null(arg1)){
    return(arg1)
  } else {
    return(arg2)
  }
}



find.largest.index <- function (func.str="sfg",lab="bpop",mat=F,mat.row=T) {
  txt <- capture.output(eval(parse(text=func.str)))
  txt <- grep(paste("^[^\\#]*",lab,"\\[",sep=""),txt,value=T)
  ind <- 0
  if(length(txt)!=0 && !mat)  ind <- gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
  if(length(txt)!=0 && mat && mat.row)  ind <- gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*,.*?\\].*",sep=""),"\\1",txt)
  if(length(txt)!=0 && mat && !mat.row)  ind <- gsub(paste("^[^\\#]*",lab,"\\[.*?,\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
  
  max(as.numeric(ind))
}
# 
#  find.largest.index("sfg","bpop")
#  find.largest.index("sfg","b")
#find.largest.index("sfg","bocc",mat=T,mat.row=T)
#  find.largest.index("sfg","x")
#  find.largest.index("sfg","a")
# 
# txt <- capture.output(eval(parse(text="sfg")))
# txt <- grep(paste("^[^\\#]*bpop","\\[",sep=""),txt,value=T)
# txt
# ind <- gsub(paste("^[^\\#]*","bpop","\\[","(\\d+)\\].*",sep=""),"\\1",txt)
# max(as.numeric(ind))
