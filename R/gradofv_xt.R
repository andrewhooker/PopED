## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradofv_xt <- function(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure){
  
  #Input: the prior FIM or (empty) and all the other things to calculate the
  #grad with
  #Return a vector that is the gradient and the global structure
  
  if((globalStructure$ofv_calc_type==1) ){#Determinant Design
    if((!globalStructure$bUseGrouped_xt)){
      returnArgs <- graddetmfxt(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
      ofv_grad <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
    } else {
      returnArgs <- graddetmfxt_ext(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
      ofv_grad <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
    }
    return
  }
  
  if((globalStructure$ofv_calc_type==2) ){#Trace of inverse 
    if((!globalStructure$bUseGrouped_xt)){
      returnArgs <-  gradtrmfxt(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
      ofv_grad <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
      return
    } else {
      fprintf('Warning: grad mf for grouped xt on A-optimal design is not implemented analytically\nNumerical difference is used instead\n')
    }
  }
  
  if((globalStructure$ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('grad mf for xt on S-optimal design not implemented yet'))
  }
  
  if((globalStructure$ofv_calc_type==4) ){#Log determinant
    if((!globalStructure$bUseGrouped_xt)){
      returnArgs <- gradlndetmfxt(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
      ofv_grad <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
    } else {
      returnArgs <- gradlndetmfxt_ext(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
      ofv_grad <- returnArgs[[1]]
      globalStructure <- returnArgs[[2]]
    }
    return
  }
  
  #All other types of criterions, i$e. Ds-optimal design, CV-optimal etc.
  iParallelN = (globalStructure$parallelSettings$bParallelSG==1) + 1 #1 if no parallel, 2 if parallel
  
  if((iParallelN == 2)){
    designsin = cell(1,0)
    it=1
  }
  
  if((!globalStructure$bUseGrouped_xt) ){#All other OFV with ungrouped xt
    m=size(ni,1)
    gdmf=matrix(1,m,size(xt,2))
    
    for(p in 1:iParallelN){
      if((p==2)){
        #Execute parallel designs
        stop("Parallel execution not yet implemented in PopED for R")
        designout = designsin
        #designout = execute_parallel(designsin,globalStructure)
      }
      if((iParallelN==1)){
        returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
        mf_tmp <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
      } else {
        if((p==1)){
          designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
        } else {
          mf_tmp = designout[[it]]$FIM
          it = it+1
        }
      }
      
      if((iParallelN==1 || p==2)){
        mft = ofv_fim(mf_tmp,globalStructure)
      }
      
      for(i in 1:m){
        if((groupsize[i]==0)){
          gdmf[i,1:ni[i]]=zeros(1,ni(i))
        } else {
          for(ct1 in 1:ni[i]){
            if((axt[i,ct1]!=0)){
              xt_plus=xt
              xt_plus[i,ct1]=xt_plus[i,ct1]+globalStructure$hgd
              
              if((iParallelN ==1)){
                returnArgs <-  mftot(model_switch,groupsize,ni,xt_plus,x,a,bpop,d,sigma,docc,globalStructure) 
                fmf_tmp <- returnArgs[[1]]
                globalStructure <- returnArgs[[2]]
              } else {
                if((p==1)){
                  designsin = update_designinlist(designsin,groupsize,ni,xt_plus,x,a,-1,0)
                } else {
                  fmf_tmp = designout[[it]]$FIM
                  it = it+1
                }
              }
              
              if((iParallelN ==1 || p==2)){
                mft_plus = ofv_fim(fmf_tmp,globalStructure)
                tmp=(mft_plus-mft)/globalStructure$hgd
                if((tmp==0)){
                  tmp=1e-12
                }
                gdmf[i,ct1]=tmp
              }
            }
          }
        }
      }
    }
    ofv_grad = gdmf
    return
  }
  if((globalStructure$bUseGrouped_xt) ){#All other OFV with grouped xt
    m=size(ni,1)
    gdmf=matrix(1,m,size(xt,2))
    
    for(p in 1:iParallelN){
      if((p==2)){
        #Execute parallel designs
        stop("Parallel execution not yet implemented in PopED for R")
        #designout = execute_parallel(designsin,globalStructure)
      }
      if((iParallelN==1)){
        returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
        mf_tmp <- returnArgs[[1]]
        globalStructure <- returnArgs[[2]]
      } else {
        if((p==1)){
          designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
        } else {
          mf_tmp = designout[[it]]$FIM
          it = it+1
        }
      }
      
      if((iParallelN==1 || p==2)){
        mft = ofv_fim(mf_tmp,globalStructure)
      }
      
      for(k in 1:max(max(max(globalStructure$G)),0)){
        tmp = matrix(1,size(xt,1),size(xt,2))*k
        inters = (globalStructure$G==tmp)
        if((sum(sum(inters))!=0) ){#If we have a time-point defined here (accord. to G)
          xt_plus = xt+globalStructure$hgd*inters
          if((iParallelN ==1)){
            returnArgs <-  mftot(model_switch,groupsize,ni,xt_plus,x,a,bpop,d,sigma,docc,globalStructure) 
            mf_plus <- returnArgs[[1]]
            globalStructure <- returnArgs[[2]]
          } else {
            if((p==1)){
              designsin = update_designinlist(designsin,groupsize,ni,xt_plus,x,a,-1,0)
            } else {
              mf_plus = designout[[it]]$FIM
              it = it+1
            }
          }
          
          if((iParallelN ==1 || p==2)){
            mft_plus = ofv_fim(mf_plus,globalStructure)
            s=(mft_plus-mft)/globalStructure$hgd
            
            if((s==0) ){#The model doesn't depend on t e$g. PD with Placebo dose, fix the xt-gradient to a small value
              s = 1e-12
            }
            
            for(i in 1:size(xt,1)){
              for(j in 1:size(xt,2)){
                if((inters[i,j]==1 && axt[i,j]!=0)){
                  gdmf[i,j]=s
                }
              }
            }
          }
        }
      }
      ofv_grad = gdmf
    }
  }
  return(list( ofv_grad= ofv_grad,globalStructure =globalStructure )) 
}
