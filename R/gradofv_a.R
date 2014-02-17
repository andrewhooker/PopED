## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradofv_a <- function(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure){

#Input: the prior FIM or (empty) and all the other things to calculate the
#grad with for a
#Return a vector that is the gradient

if((!isempty(globalStructure$ed_penalty_pointer)) ){#If a penalty function is used
    #[ni, xt, model_switch, x, a, bpop, n, d, maxxt, minxt, maxa,mina]=downsizing_general_design
    na = size(a,2)
    m=size(ni,1)
    gdmf=zeros(m,na)
     returnArgs <- ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,globalStructure$covd,sigma,docc,globalStructure) 
fmf <- returnArgs[[1]]
mft <- returnArgs[[2]]
globalStructure <- returnArgs[[3]]
    for(i in 1:m){
        if((groupsize[i]==0)){
            gdmf[i,1:ni[i]]=zeros(1,ni(i))
        } else {
            for(ct1 in 1:na){
                if((aa[i,ct1]!=0)){
                    a_plus=a
                    a_plus[i,ct1]=a_plus[i,ct1]+globalStructure$hgd
                     returnArgs <-  ed_mftot(model_switch,groupsize,ni,xt,x,a_plus,bpop,d,globalStructure$covd,sigma,docc,globalStructure) 
fmf_tmp <- returnArgs[[1]]
mft_plus <- returnArgs[[2]]
globalStructure <- returnArgs[[3]]
                    tmp=(mft_plus-mft)/globalStructure$hgd
                    if((tmp==0)){
                        gdmf[i,ct1]=1e-12
                    } else {
                        gdmf[i,ct1]=tmp
                    }
                }
            }
        }
    }
    ofv_grad = gdmf
    return
}

if((globalStructure$ofv_calc_type==1) ){#determinant
    if((!globalStructure$bUseGrouped_a)){
         returnArgs <- graddetmfa(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
ofv_grad <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    } else {
         returnArgs <- graddetmfa_ext(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
ofv_grad <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    }
    return
}

if((globalStructure$ofv_calc_type==2) ){#A-Optimal Design
    if((!globalStructure$bUseGrouped_a)){
         returnArgs <-  gradtrmfa(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
ofv_grad <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
        return
    } else {
        fprintf('Warning: grad mf for grouped a on A-optimal design is not implemented analytically\nNumerical difference is used instead\n')
    }
}

if((globalStructure$ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('S-optimal design for a is not implemented yet!'))
}

if((globalStructure$ofv_calc_type==4) ){#Log determinant
    if((!globalStructure$bUseGrouped_a)){
         returnArgs <- gradlndetmfa(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
ofv_grad <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
    } else {
         returnArgs <- gradlndetmfa_ext(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,globalStructure) 
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

if((!globalStructure$bUseGrouped_a) ){#If not a grouped criterion
    na = size(a,2)
    m=size(ni,1)
    gdmf=zeros(m,na)

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
                for(ct1 in 1:na){
                    if((aa[i,ct1]!=0)){
                        a_plus=a
                        a_plus[i,ct1]=a_plus[i,ct1]+globalStructure$hgd

                        if((iParallelN ==1)){
                             returnArgs <-  mftot(model_switch,groupsize,ni,xt,x,a_plus,bpop,d,sigma,docc,globalStructure) 
fmf_tmp <- returnArgs[[1]]
globalStructure <- returnArgs[[2]]
                        } else {
                            if((p==1)){
                                designsin = update_designinlist(designsin,groupsize,ni,xt,x,a_plus,-1,0)
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

if((globalStructure$bUseGrouped_a) ){#If grouped a
    m=size(ni,1)
    gdmf=matrix(1,m,size(a,2))

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

        for(k in 1:max(max(max(globalStructure$Ga)),0)){
            tmp = matrix(1,size(a,1),size(a,2))*k
            inters = (globalStructure$Ga==tmp)
            if((sum(sum(inters))!=0) ){#If we have a time-point defined here (accord. to G)
                a_plus = a+globalStructure$hgd*inters
                if((iParallelN ==1)){
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

                if((iParallelN ==1 || p==2)){
                    mft_plus = ofv_fim(mf_plus,globalStructure)
                    s=(mft_plus-mft)/globalStructure$hgd

                    if((s==0) ){#The model doesn't depend on a e$g.
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
        ofv_grad = gdmf
    }
}
return(list( ofv_grad= ofv_grad,globalStructure =globalStructure )) 
}
