## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

graddetmf_ext <- function(model_switch,aX,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db,lndet=FALSE,gradxt=FALSE){
  
  n = get_fim_size(poped.db)
  m=size(ni,1)
  if (gradxt==FALSE) {
    gdmf = matrix(1,m,size(a,2))
    G_X = poped.db$design_space$G_a
  } else {
    gdmf = matrix(1,m,size(xt,2))
    G_X = poped.db$design_space$G_xt
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
      #If we have a prior
      if(all(size(poped.db$settings$prior_fim) == size(mft))){
        mft = mft + poped.db$settings$prior_fim
      }
      imft=inv(mft)
      if((is.infinite(imft[1,1]))){
        imft = zeros(size(mft))
      }
    }
    
    a0 = a
    xt0 = xt
    for(k in 1:max(max(max(G_X)),0)){
      inters = (G_X == k)
      if((sum(sum(inters))!=0) ){#If we have a covariate or time-point defined here (accord. to G_X)
	if (gradxt==FALSE) {
          a = a0+poped.db$settings$hgd*inters
        } else {
          xt = xt0+poped.db$settings$hgd*inters
        }
        if((iParallelN==1)){
          returnArgs <-  mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
          mf_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
        } else {
          if((p==1)){
            designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
          } else {
            mf_plus = designout[[it]]$FIM
            it = it+1
          }
        }
        
        if((iParallelN==1 || p==2)){
          #If we have a prior
          if(all(size(poped.db$settings$prior_fim)==size(mft))){
            mf_plus = mf_plus + poped.db$settings$prior_fim
          }
          ir=(mf_plus-mft)/poped.db$settings$hgd
          
          s=0 #Calc the tr(A^-1 * dA/dX) for some X
          for(ct2 in 1:n){
            s=s+imft[ct2,,drop=F]%*%ir[,ct2,drop=F]
          }
          
          if((s==0) ){#The model doesn't depend on a or xt, e$g. PD is only dependent on time and not dose, fix the a-gradient to a small value or PD with Placebo dose, fix the xt-gradient to a small value
            s = 1e-12
          }
          gdmf[inters==1 & aX!=0] = s
          # for(i in 1:size(a,1)){
          #   for(j in 1:size(a,2)){
          #     if((inters[i,j]==1 && aX[i,j]!=0)){
          #       gdmf[i,j]=s
          #     }
          #   }
          # }
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


