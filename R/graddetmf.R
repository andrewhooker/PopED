## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

graddetmf <- function(model_switch,aX,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db,lndet=FALSE,gradxt=FALSE){
  
  n = get_fim_size(poped.db)
  m=size(ni,1)
  if (gradxt == FALSE) {
    gdmf=matrix(1,m,size(a,2))
  } else {
    gdmf=matrix(1,m,size(xt,2))
  }
  
  iParallelN = (poped.db$settings$parallel$bParallelSG==1) + 1 #1 if no parallel, 2 if parallel
  
  if((iParallelN == 2)){
    designsin = cell(1,0)
    it = 1
  }
  for(p in 1:iParallelN){
    if((p==2)){
      #Execute parallel designs
      stop("Parallel execution not yet implemented in PopED for R")
      designout = designsin
      #designout = execute_parallel(designsin,poped.db)
    }
    if((iParallelN==1)){
      returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
      mft <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
    } else {
      if((p==1)){
        designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
      } else {
        mft = designout[[it]]$FIM
        it = it+1
      }
    }
    
    if((iParallelN==1 || p==2)){
      if(all(size(poped.db$settings$prior_fim) == size(mft))){
        mft = mft + poped.db$settings$prior_fim
      }
      imft=inv(mft)
      if((is.infinite(imft[1,1]))){
        imft = zeros(size(mft))
      }
    }
    
    for(i in 1:m){
      if(groupsize[i]==0){
        gdmf[i,1:ni[i]]=zeros(1,ni(i))
      } else {
        if((!isempty(x))){
          x_i = t(x[i,,drop=F])
        } else {
          x_i =  zeros(0,1)
        }
        if((!isempty(a))){
          a_i = t(a[i,,drop=F])      
        } else {
          a_i =  zeros(0,1)
        }
        if((iParallelN ==1)){
          returnArgs <- mf_all(t(model_switch[i,1:ni[i,drop=F],drop=F]),t(xt[i,1:ni[i,drop=F],drop=F]),x_i,a_i,bpop,d,sigma,docc,poped.db) 
          mf_tmp <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
        } else {
          if((p==1)){
            designsin = update_designinlist(designsin,1,ni,xt,x,a,-1,i)
          } else {
            mf_tmp = designout[[it]]$FIM
            it = it+1
          }
        }
        mfb=groupsize[i]*mf_tmp
        
	a0 = a
	xt0 = xt
	nCtl = ifelse(gradxt==FALSE, size(poped.db$design$a,2), ni[i])
        for(ct1 in 1:nCtl){
          if(aX[i,ct1]!=0){
	    if (gradxt==FALSE) {
              a=a0
              a[i,ct1]=a[i,ct1]+poped.db$settings$hgd
            } else {
              xt=xt0
              xt[i,ct1]=xt[i,ct1]+poped.db$settings$hgd
            }
            if((iParallelN ==1)){
              returnArgs <- mf_all(t(model_switch[i,1:ni[i,drop=F],drop=F]),t(xt[i,1:ni[i,drop=F],drop=F]),x_i,t(a[i,,drop=F]),bpop,d,sigma,docc,poped.db) 
              mf_tmp <- returnArgs[[1]]
              poped.db <- returnArgs[[2]]
            } else {
              if((p==1)){
                designsin = update_designinlist(designsin,1,ni,xt,x,a,-1,i)
              } else {
                mf_tmp = designout[[it]]$FIM
                it = it+1
              }
            }
            if((iParallelN ==1 || p==2)){
              mf_plus = groupsize[i]*mf_tmp
              ir=(mf_plus-mfb)/poped.db$settings$hgd
              s=0 #Calc the tr(A^-1 * dA/dX) for some X
              for(ct2 in 1:n){
                s=s+imft[ct2,,drop=F]%*%ir[,ct2,drop=F]
              }
              
              if((s==0) ){#The model doesn't depend on design variable (a or t), e.g. PD is only dependent on time and not dose, fix the a-gradient to a small value or PD with Placebo dose, fix the xt-gradient to a small value.
                s = 1e-12
              }
              gdmf[i,ct1]=s
            }
          }
        }
      }
    }
  }
  if (lndet == FALSE) {
    ret=gdmf*det(mft)
  } else {
    ret=gdmf
  }
  return(list(ret = ret, poped.db = poped.db)) 
}

