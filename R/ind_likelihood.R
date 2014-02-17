#Calculates the individual ln likelihood
#Written for PopED by JN
ind_likelihood <- function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,globalStructure,lC,det_res_var,b_ind){
  
  if((!bUDDLike)){
    fg0=feval(globalStructure$fg_pointer,x,a,bpop,b_ind,bocc_ind)
    ipred = feval(globalStructure$ff_pointer,model_switch,xt_ind,fg0,globalStructure)
    res = data-ipred #Individual residuals
    
    if((bInter)){
      #For Cases WITH interaction, linearize around eta = eta^
      #eps = zeros(size(tdata),1),size(sigma,2))
      eps = zeros(size(t(data),1),size(sigma,2))
      h = LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,globalStructure) #The covariance for this individual
      res_var = diag_matlab(diag_matlab(h%*%sigma%*%t(h)))
      lC=solve(chol(res_var),diag_matlab(length(res_var)))
      det_res_var = det(res_var)
    } else {
      #Cases WITHOUT interaction, linearize around eta = 0
      #  h = LinMatrixH(tdata,cdata,theta,zeros(size(eta)),eps) #The covariance for this individual
    }
    
    R=(t(res)%*%lC)
    li = -1/2*log(det_res_var)-1/2*R%*%t(R) # + const
    
  } else {
    #%UDD likelihood
    #li=sum(model(tdata,cdata,theta,eta))
    stop("User defined likelihood not implemented for PopED in R")  
  }
  return( li ) 
}
