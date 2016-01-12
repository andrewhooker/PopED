
# very sparse grid to evaluate (4 points for each design valiable)
output <- a_line_search(poped.db, opt_xt=TRUE, opt_a=TRUE, ls_step_size=4)

\dontrun{  
  
  # longer run time
  output <- a_line_search(poped.db,opt_xt=TRUE)
  
  # output to a text file
  output <- a_line_search(poped.db,opt_xt=TRUE,out_file="tmp.txt")
  
}

