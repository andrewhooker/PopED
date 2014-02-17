
hesskalpha2 <- function(alpha, model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,globalStructure,ha,Engine){
  #D2KALPHA2 calculates the hessian of k with respect to alpha
  #   Detailed explanation goes here
  returnArgs <- log_prior_pdf(alpha,bpopdescr,ddescr) 
  p <- returnArgs[[1]]
  gradp <- returnArgs[[2]]
  hessp <- returnArgs[[3]]
  #get dF/dAlpha and fim
  returnArgs <- dfimdalpha(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,globalStructure,ha) 
  d_fim <- returnArgs[[1]]
  fim <- returnArgs[[2]]
  try({
    L=chol(fim)
    # calc inverse
    iL=solve(L,diag_matlab(length(L)))
    ifim=t(iL)%*%iL
    # calc trace of iF*dF/dAlpha(i)*iF*dF/dAlpha(j)
    for(i in 1:size(d_fim,3)){
      for(j in 1:i){
        tigi[i,j]=trace_matrix(ifim*d_fim[,,i]*ifim*d_fim[,,j])
        tigi[j,i]=tigi[i,j]
      }
    }
    d2=d2fimdalpha2(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,globalStructure,1e-4)
    d2Fim=reshape_matlab(d2,length(fim)^2,length(hessp)^2)
    d2logdfim=t(d2Fim)%*%ifim
    hess=-(hessp+reshape_matlab(d2logdfim,length(hessp),length(hessp))-tigi)
  })
#   catch
#   if((Engine$Type==1)){
#     exception = lasterror
#     if(exception$identifier=='MATLAB:posdef'){
#       hess=zeros(length(alpha))
#     } else {
#       rethrow(exception)
#     }
#   } else {
#     exception = lasterr
#     if(exception$identifier==''){
#       hess=zeros(length(alpha))
#     } else {
#       stop(sprintf(exception))
#     }
#   }

return(  hess  ) 
}

