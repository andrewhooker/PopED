graddetmfa_ext <- function(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure){

n = get_fim_size(globalStructure)

m=size(ni,1)
gdmf=matrix(1,m,size(a,2))

iParallelN = (globalStructure$parallelSettings$bParallelSG==1) + 1 #1 if no parallel, 2 if parallel

if((iParallelN == 2)){
    designsin = cell(1,0)
}


it = 1
for(p in 1:iParallelN){
    if((p==2)){
        #Execute parallel designs
      stop("Parallel execution not yet implemented in PopED for R")
      designout = designsin
      #designout = execute_parallel(designsin,globalStructure)
    }
    
    if((iParallelN==1)){
         returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
mft <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    } else {
        if((p==1)){
            designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
        } else {
            mft = designout[[it]]$FIM
            it = it+1
        }
    }
    
    if((iParallelN==1 || p==2)){
        if(all(size(globalStructure$prior_fim)==size(mft))){
            mft = mft + globalStructure$prior_fim
        }
        imft=inv(mft)
        if((isinf(imft[1,1]))){
            imft = zeros(size(mft))
        }
    }
    
    
    for(k in 1:max(max(max(globalStructure$Ga)),0)){
        tmp = matrix(1,size(a,1),size(a,2))*k
        inters = (globalStructure$Ga==tmp)
        if((sum(sum(inters))!=0) ){#If we have a covariate defined here (accord. to Ga)
            a_plus = a+globalStructure$hgd*inters
            if((iParallelN==1)){
                 returnArgs <-  mftot(model_switch,groupsize,ni,xt,x,a_plus,bpop,d,sigma,docc,globalStructure) 
mf_plus <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
            } else {
                if((p==1)){
                    designsin = update_designinlist(designsin,groupsize,ni,xt,x,a_plus,-1,0)
                } else {
                    mf_plus = designout[[it]]$FIM
                    it = it+1
                }
            }
            #If we have a prior
            
            if((iParallelN ==1 || p==2)){
                if(all(size(globalStructure$prior_fim)==size(mft))){
                    mf_plus = mf_plus + globalStructure$prior_fim
                }
                ir=(mf_plus-mft)/globalStructure$hgd

                s=0 #Calc the tr(A^-1 * dA/dX) for some X
                for(ct2 in 1:n){
                    s=s+imft[ct2,,drop=F]%*%ir[,ct2,drop=F]
                }

                if((s==0) ){#The model doesn't depend on a, e$g. PD is only dependent on time and not dose, fix the a-gradient to a small value
                    s = 1e-12
                }

                for(i in 1:size(a,1)){
                    for(j in 1:size(a,2)){
                        if((inters[i,j]==1 && aa[i,j]!=0)){
                            gdmf[i,j]=s
                        }
                    }
                }
            }
        }
    }
}
ret=gdmf*det(mft)
return(list( ret= ret,globalStructure=globalStructure)) 
}




