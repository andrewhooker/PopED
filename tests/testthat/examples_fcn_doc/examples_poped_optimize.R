# warfarin optimize model

\dontrun{
  
  ##############
  # Optimization
  ##############
  
  # below are a number of ways to optimize the problem
  
  
  # RS+SG+LS optimization of sample times
  # optimization with just a few iterations
  # only to check that things are working
  output <- poped_optimize(poped.db,opt_xt=T,
                           rsit=5,sgit=5,ls_step_size=5)
  
  # RS+SG+LS optimization of sample times 
  # (longer run time than above but more likely to reach a maximum)
  output <- poped_optimize(poped.db,opt_xt=T)
  get_rse(output$fmf,output$poped.db)
  plot_model_prediction(output$poped.db)
  
  # MFEA optimization with only integer times allowed
  mfea.output <- poped_optimize(poped.db,opt_xt=1,
                                bUseExchangeAlgorithm=1,
                                EAStepSize=1)
  get_rse(mfea.output$fmf,mfea.output$poped.db)
  plot_model_prediction(mfea.output$poped.db)
  
  # Examine efficiency of sampling windows
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=0.5)
  plot_efficiency_of_windows(mfea.output$poped.db,xt_windows=1)
  
  # Random search (just a few samples here)
  rs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,rsit=20,
                              bUseRandomSearch= 1,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0)
  get_rse(rs.output$fmf,rs.output$poped.db)
  
  # line search, DOSE and sample time optimization
  ls.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 0,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 1,
                              ls_step_size=10)
  
  # Stochastic gradient search, DOSE and sample time optimization
  sg.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1, 
                              bUseRandomSearch= 0,
                              bUseStochasticGradient = 1,
                              bUseBFGSMinimizer = 0,
                              bUseLineSearch = 0,
                              sgit=20)
  
  # BFGS search, DOSE and sample time optimization
  bfgs.output <- poped_optimize(poped.db,opt_xt=1,opt_a=1,
                                bUseRandomSearch= 0,
                                bUseStochasticGradient = 0,
                                bUseBFGSMinimizer = 1,
                                bUseLineSearch = 0)
  
}
