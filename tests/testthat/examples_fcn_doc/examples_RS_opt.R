
# Just a few iterations, optimize on DOSE and sample times using the full FIM
out_1 <- RS_opt(poped.db,opt_xt=1,opt_a=1,rsit=3,fim.calc.type=0)

\dontrun{
  
  # More iterations
  rs.output <- RS_opt(poped.db)
  
}  
