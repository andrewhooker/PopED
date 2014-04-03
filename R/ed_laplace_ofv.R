#' Evaluate the expectation of determinant the Fisher Information Matrix (FIM) using the Laplace approximation.
#' 
#' Compute the expectation of the \code{det(FIM)} using the Laplace approximation to the expectation.  
#' Computations are made based on the model, parameters, distributions of parameter uncertainty, 
#' design and methods defined in the 
#' PopED database or as arguments to the funciton. 
#' 
#' @inheritParams evaluate.fim
#' @inheritParams create.poped.database
#' @inheritParams Doptim
#' @param xtopto the xtopto value
#' @param xopto the xopto value
#' @param optxt If sampling times are optimized
#' @param opta If continuous design variables are optimized
#' @param aopto the aopto value
#' 
#' @param return_gradient Should the gradient be returned.
#' 
#' @return The FIM and the hessian of the FIM.
#' 
#' @family FIM
#' @family E-family
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker
## right now function only works for normal and log-normal priors

ed_laplace_ofv <- function(x,optxt, opta, model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc, poped.db,return_gradient=FALSE){
  if(any(ddescr[,1,drop=F]!=0&ddescr[,1,drop=F]!=4)){
    stop(sprintf('Only lognormal prior is supported for random effects!')) 
  }
  if(any(bpopdescr[,1,drop=F]!=0 & bpopdescr[,1,drop=F]!=4 & bpopdescr[,1,drop=F]!=1)){
    stop(sprintf('Only Normal and lognormal priors are supported for fixed effects!')) 
  }
  Engine = list(Type=1,Version=version$version.string)
  
  x2=x
  if(!isempty(x)){
    if(optxt){
      notfixed=poped.db$design$minxt!=poped.db$design$maxxt
      if(poped.db$bUseGrouped_xt){
        xtopto[notfixed]=x[poped.db$design$G[notfixed]]
        x[1:numel(unique(poped.db$design$G[notfixed]))]=matrix(0,0,0)
      } else {
        xtopto[notfixed]=x[1:numel(xtopto[notfixed])]
        x=x[-c(1:numel(xtopto[notfixed]))]
      }
    }
    if(opta){
      notfixed=poped.db$design$mina!=poped.db$design$maxa
      if(poped.db$bUseGrouped_a){
        aopto[notfixed]=x[poped.db$design$Ga[notfixed]]
      } else {
        aopto[notfixed]=x
      }
    }
    x=x2    
  }
  
  alpha_k=matrix(c(bpopdescr[bpopdescr[,1]!=0,2], ddescr[ddescr[,1]!=0,2]),ncol=1,byrow=T)
  
  
  ## do log transformation of parameters for unconstrained optimization 
  ## (silly, and should be changed, distribution should define allowed values)
  d_index=(sum(bpopdescr[,1]!=0)+1):length(alpha_k)
  bpop_index=1:sum(bpopdescr[,1]!=0)
  unfixed_bpop <- bpopdescr[bpopdescr[,1]!=0,]
  exp_index=c(unfixed_bpop[,1]==4,d_index==d_index)
  alpha_k_log=alpha_k
  alpha_k_log[exp_index]=log(alpha_k[exp_index])
  #alpha_k_log[d_index]=log(alpha_k[d_index])
  trans  <- function(x){ 
    x[exp_index] <- exp(x[exp_index])
    return(x)
  }
  #trans  <- function(x) matrix(c(x[bpop_index],exp(x[d_index])),ncol=1,byrow=T)
  
  #calc initial function value and gradient
  returnArgs <- calc_k(alpha_k,model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                       return_gradient=T) 
  f_k <- returnArgs[[1]]
  gf_k <- returnArgs[[2]]
  
  #transform gradient for ds (log(alpha))
  #   if(!isempty(d_index)){
  #     gf_k[d_index]=gf_k[d_index]*exp(alpha_k_log[d_index])
  #   }
  if(!isempty(exp_index)){
    gf_k[exp_index]=gf_k[exp_index]*exp(alpha_k_log[exp_index])
  }
  if(isnan(f_k)){
    f=0
    gf=zeros(size(x))
    return
  }
  #initialize optimization variables
  dim=length(alpha_k)
  
  H_k=diag_matlab(dim)
  B_k=H_k
  niter=0
  while(norm(gf_k,type="2")>0.001){ 	# while inner conv. krit. not met{
    #determine search direction for line search
    p_k=-H_k%*%gf_k
    f_name  <- "calc_k"
    #f_options <- list(trans(alpha),model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine)
    f_options <- list("replace",model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine)
    returnArgs <- line_search_uc(alpha_k_log,f_k,gf_k,p_k,f_name,f_options,exp_index)
    alpha_k1_log <- returnArgs[[1]]
    f_k1 <- returnArgs[[2]]
    s_k=alpha_k1_log-alpha_k_log
    if(max(abs(t(s_k))/max(matrix(c(t(alpha_k1_log), matrix(1,1,length(t(alpha_k1_log)))),nrow=2,byrow=T)))<1e-3){ 
      # check this that it is the same as in matlab
      break
    }
    returnArgs <- calc_k(trans(alpha_k1_log),model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                         return_gradient=T) 
    f_k1 <- returnArgs[[1]]
    gf_k1 <- returnArgs[[2]]
    #transform gradient for ds (log(alpha))
    #     if(!isempty(d_index)){
    #       gf_k1[d_index]=gf_k1[d_index]*exp(alpha_k1_log(d_index))
    #     }
    if(!isempty(exp_index)){
      gf_k1[exp_index]=gf_k1[exp_index]*exp(alpha_k1_log[exp_index])
    }
    y_k=gf_k1-gf_k
    rho_k=1/(t(y_k)%*%s_k)
    rho_k <- rho_k[,,drop=T]
    if((t(y_k)%*%s_k)/(-t(gf_k)%*%s_k) > .Machine$double.eps){
      H_k=(diag_matlab(dim)-rho_k*s_k%*%t(y_k))%*%H_k%*%(diag_matlab(dim)-rho_k*y_k%*%t(s_k))+rho_k*s_k%*%t(s_k)
    }
    alpha_k_log=alpha_k1_log
    gf_k=gf_k1
    f_k=f_k1
    niter=niter+1
  }
  alpha_k=trans(alpha_k_log)
  #if the number of iterations is smaller than the dimension of the problem
  #we have to calculate the hessian explicitly
  if((niter<length(B_k)||poped.db$iEDCalculationType==1)){
    hess=hesskalpha2(alpha_k, model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,1e-6,Engine)
    detHessPi=det(hess)*(2*pi)^(-length(hess))
  } else {
    temp=matrix(1,size(gf_k))
    #temp[d_index]=1/alpha_k[d_index]
    temp[exp_index]=1/alpha_k[exp_index]
    
    iH_k=inv(H_k)
    hess=iH_k*(temp%*%t(temp))
    detHessPi=det(hess)*(2*pi)^(-length(hess))
  }
  f=Re(-exp(-f_k)/sqrt(detHessPi))
  
  if(return_gradient){
    bpop=bpopdescr[,2,drop=F]
    bpop[bpopdescr[,1,drop=F]!=0]=alpha_k[1:sum(bpopdescr[,1,drop=F]!=0),drop=F]
    d=ddescr[,2,drop=F]
    d[ddescr[,1]!=0]=alpha_k[sum(bpopdescr[,1,drop=F]!=0)+1:end,drop=F]
    d=getfulld(d,covd)
    
    gradxt=matrix(0,0,0)
    grada=matrix(0,0,0)
    if((optxt==TRUE)){
      notfixed=poped.db$design$minxt!=poped.db$design$maxxt
      gradxt=-gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
      gradxt=gradxt(notfixed)
      if(poped.db$bUseGrouped_xt){
        index=unique(poped.db$design$G)
        gradxt=gradxt(index)
      }
    }
    if((opta==TRUE)){
      notfixed=poped.db$design$mina!=poped.db$design$maxa
      grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)                
      grada=grada(notfixed)
      if(poped.db$bUseGrouped_a){
        index=unique(poped.db$design$Ga)
        grada=grada(index)
      }
    }
    dkdxt=matrix(c(gradxt,grada),nrow=1,byrow=T)
    
    h_alpha=1e-4
    tensor=array(0,dim=c(length(alpha_k), length(alpha_k), length(dkdxt)))
    for(i in 1:length(alpha_k)){
      for(j in 1:i){
        alpha_plus_plus=alpha_k
        alpha_plus_plus[i]=alpha_plus_plus[i]+h_alpha
        alpha_plus_plus[j]=alpha_plus_plus[j]+h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_plus_plus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design$minxt!=poped.db$design$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$bUseGrouped_xt){
            index=unique(poped.db$design$G)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design$mina!=poped.db$design$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$bUseGrouped_a){
            index=unique(poped.db$design$Ga)
            grada=grada(index)
          }
        }
        dkdxt_plus_plus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        alpha_minus_plus=alpha_k
        alpha_minus_plus[i]=alpha_minus_plus[i]-h_alpha
        alpha_minus_plus[j]=alpha_minus_plus[j]+h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_minus_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_minus_plus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design$minxt!=poped.db$design$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$bUseGrouped_xt){
            index=unique(poped.db$design$G)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design$mina!=poped.db$design$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$bUseGrouped_a){
            index=unique(poped.db$design$Ga)
            grada=grada(index)
          }
        }
        dkdxt_minus_plus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        alpha_plus_minus=alpha_k
        alpha_plus_minus[i]=alpha_plus_minus[i]+h_alpha
        alpha_plus_minus[j]=alpha_plus_minus[j]-h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus_minus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_plus_minus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design$minxt!=poped.db$design$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$bUseGrouped_xt){
            index=unique(poped.db$design$G)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design$mina!=poped.db$design$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$bUseGrouped_a){
            index=unique(poped.db$design$Ga)
            grada=grada(index)
          }
        }
        dkdxt_plus_minus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        alpha_minus_minus=alpha_k
        alpha_minus_minus[i]=alpha_minus_minus[i]-h_alpha
        alpha_minus_minus[j]=alpha_minus_minus[j]-h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_minus_minus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_minus_minus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design$minxt!=poped.db$design$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$bUseGrouped_xt){
            index=unique(poped.db$design$G)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design$mina!=poped.db$design$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$bUseGrouped_a){
            index=unique(poped.db$design$Ga)
            grada=grada(index)
          }
        }
        dkdxt_minus_minus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        tensor[i,j,]=((dkdxt_plus_plus-dkdxt_plus_minus-dkdxt_minus_plus+dkdxt_minus_minus))/(4*h_alpha^2)
        tensor[j,i,]=tensor[i,j,]
      }
    }
    ddetHessdxt=zeros(length(dkdxt),1)
    for(i in 1:length(dkdxt)){
      ddetHessdxt[i]=detHessPi*trace_matrix(inv(hess)*(-tensor[,,i]))
    }
    
    gf=Re(exp(-f_k)*(2*detHessPi*dkdxt+ddetHessdxt)/(2*detHessPi^(3/2)))
  }
  return(list( f= f,gf=gf)) 
}

calc_k <- function(alpha, model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,return_gradient=F){
  bpop=bpopdescr[,2,drop=F]
  bpop[bpopdescr[,1,drop=F]!=0]=alpha[1:sum(bpopdescr[,1,drop=F]!=0),drop=F]
  d=ddescr[,2,drop=F]
  d[ddescr[,1]==4]=alpha[(sum(bpopdescr[,1,drop=F]!=0)+1):length(alpha),drop=F]
  d=getfulld(d,covd)
  retargs=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,poped.db)
  fim <- retargs$ret
  if((!return_gradient)){
    k=-log_prior_pdf(alpha, bpopdescr, ddescr)-log(det(fim))
    grad_k=matrix(0,0,0)
  } else {
    returnArgs <- log_prior_pdf(alpha, bpopdescr, ddescr,return_gradient=T) 
    logp <- returnArgs[[1]]
    grad_p <- returnArgs[[2]]
    returnArgs <- dfimdalpha(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped.db,1e-6) 
    d_fim <- returnArgs[[1]]
    fim <- returnArgs[[2]]
    ifim <- inv(fim)
    dim(ifim) <- c(length(ifim),1)
    gradlogdfim=t(reshape_matlab(d_fim,length(fim),length(grad_p)))%*%ifim
    grad_k=-(gradlogdfim+grad_p)
    ## if not positive definite set grad_k=zeros(length(alpha),1)
    
    k=-logp-log(det(fim))
  }
  return(list( k= k, grad_k= grad_k)) 
}
