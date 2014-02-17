line_search_uc <- function(x0, f0, g0, d, f_handle,f_options,exp_index){
  dim=length(x0)
  ls_options <- c()
  ls_options$ls_stepmax=1 #max step length for line search
  ls_options$ls_delta_alpha=1e-4 #convergence criterion on alpha
  ls_options$ls_fdecreas=1e-4 #sufficent decrease for line search
  if((norm(d)>ls_options$ls_stepmax)){
    d=d*ls_options$ls_stepmax/norm(d)
  }
  #compute lambda_min
  lambda_min=ls_options$ls_delta_alpha/max(abs(d)/max(matrix(c(abs(x0),matrix(1,dim,1)),nrow=1,byrow=T)))
  lambda=1 #newton step
  slope=t(g0)%*%d #slope at beginning of line search
  
  #trans  <- function(x) matrix(c(x[bpop_index],exp(x[d_index])),ncol=1,byrow=T)
  
  trans  <- function(x){ 
    x[exp_index] <- exp(x[exp_index])
    return(x)
  }
  
  while(TRUE){
    x1=x0+lambda*d
    #evaluate function
    f_options[[1]]=trans(x1)
    fval1=do.call(f_handle,f_options)
    fval1 <- fval1$k
    if(lambda<lambda_min){
      x_min=x0
      f_min=f0
      break
    } else if (fval1<f0+ls_options$ls_fdecreas*lambda*slope){
      x_min=x1
      f_min=fval1
      break
    } else {
      if(lambda==1){
        lambda_tmp=-slope/(2*(fval1-f0-slope))
      } else {
        rhs1=fval1-f0-lambda*slope
        rhs2=fval2-f0-lambda2*slope
        a=(rhs1/lambda^2-rhs2/lambda2^2)/(lambda-lambda2)
        b=(-lambda2*rhs1/lambda^2+lambda*rhs2/lambda2^2)/(lambda-lambda2)
        if(a==0){
          lambda_tmp=-slope/(2*b)
        } else {
          disc=b^2-3*a*slope
          if(disc<0){
            lambda_tmp=0.5*lambda
          } else if (b<0){
            lambda_tmp=(-b*sqrt(disc))/(3*a)
          } else {
            lambda_tmp=-slope/(b+sqrt(disc))
          }
          lambda_tmp=min(matrix(c(0.5*lambda, lambda_tmp),nrow=1,byrow=T))
        }
      }
      lambda2=lambda
      lambda=max(lambda_tmp,0.1*lambda)
      fval2=fval1
    }
  }
  return(list( x_min= x_min, f_min= f_min)) 
}

