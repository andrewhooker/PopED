


# Examine efficiency of sampling windows at plus/minus 0.5 hours from 
# sample points in the design
plot_efficiency_of_windows(poped.db,xt_windows=0.5)


if(interactive()){  

  plot_efficiency_of_windows(poped.db,
                             xt_plus=c( 0.5,1,2,1,2,3,7,1),
                             xt_minus=c( 0.1,2,5,4,2,3,6,2))
  
  plot_efficiency_of_windows(poped.db,xt_windows=c( 0.5,1,2,1,2,3,7,1))
  
  
  plot_efficiency_of_windows(poped.db,
                             xt_plus=c( 0.5,1,2,1,2,3,7,1),
                             xt_minus=c( 0.1,2,5,4,2,3,6,2),
                             y_rse=FALSE)
  
  plot_efficiency_of_windows(poped.db,
                             xt_plus=c( 0.5,1,2,1,2,3,7,1),
                             xt_minus=c( 0.1,2,5,4,2,3,6,2),
                             y_eff=FALSE)
}

