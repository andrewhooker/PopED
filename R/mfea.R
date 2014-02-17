## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

mfea <- function(globalStructure,model_switch,ni,xt,x,a,bpopdescr,ddescr,maxxt,minxt,maxa,mina,fmf,dmf){
  
  if((globalStructure$EACriteria!=1)){
    stop(sprintf('The criteria that can be used is Modified Fedorov Exchange Algorithm (EACriteria = 1)'))
  }
  globalStructure$parallelSettings$bParallelMFEA
  if((globalStructure$parallelSettings$bParallelMFEA == 1)){
    designsin = cell(1,0) # temporary solution
    iParallelN = 2
  } else {
    iParallelN = 1
  }
  
  m=size(ni,1)
  
  rho=globalStructure$EAConvergenceCriteria
  xt_current = xt
  a_current = a
  x_current = x
  bpop=bpopdescr[,2,drop=F]
  d=getfulld(ddescr[,2,drop=F],globalStructure$covd)
  doccfull = getfulld(globalStructure$docc[,2,drop=F],globalStructure$covdocc)
  
  stepsize = globalStructure$EAStepSize
  numpoints = globalStructure$EANumPoints
  
  delta_max = rho+0.1
  it = 0
  
  previous_index_i = 0
  previous_index_ct = 0
  previous_type = 0
  
  # if((globalStructure$bShowGraphs)){
  #     cfg = figure(1)
  #     if((globalStructure$Engine$Type==1) ){#Matlab
  #         set(cfg,'Name','Exchange Algorithm')
  #     }
  # }
  
  while(delta_max>rho){
    dmfvector <- c()
    best_it = 0
    if((it==0)){
      if((globalStructure$d_switch)){
        if((iParallelN==1)){ #IF NOT PARALLEL MFEA
          returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_current,a_current,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
          nfmf <- returnArgs[[1]]
          globalStructure <- returnArgs[[2]]
          optofv=ofv_fim(nfmf,globalStructure)
        } else { ##STORE DESIGN FOR PARALLEL EXEC.
          designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_current,x_current,a_current,-1,0)
        }
      } else {
        returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_current,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
        nfmf <- returnArgs[[1]]
        optofv <- returnArgs[[2]]
        globalStructure <- returnArgs[[3]]
      }
    }
    if((iParallelN==1)){
      if((it==0 && (optofv<=0 || isnan(optofv)))){
        fprintf('The OFV is negative or zero, new initial values might be needed. OFV : %g\n',optofv)
      }
    }
    delta_max = 0
    t_max = 0
    
    bFirstInGroup_a =  TRUE
    bFirstInGroup_x =  TRUE
    bFirstInGroup_xt = TRUE
    
    for(i in 1:m ){#For all groups
      if((globalStructure$optsw[2]==TRUE) ){#Optimize over samples
        for(p  in 1:iParallelN){
          if((p==2)){
            stop("Parallel execution not yet implemented in PopED for R")
            #designout = execute_parallel(designsin,globalStructure)
            designout = designsin
            des_num = 1
            if((it==0)){
              nfmf = designout[[des_num]]$FIM
              optofv = ofv_fim(nfmf,globalStructure)
              des_num = des_num +1
              if((it==0 && (optofv<=0 || isnan(optofv)))){
                fprintf('The OFV is negative or zero, new initial values might be needed. OFV : %g\n',optofv)
              }
            }
          }
          
          if((globalStructure$bUseGrouped_xt)){
            if((bFirstInGroup_xt)){
              bFirstInGroup_xt = FALSE
              for(k in 1:max(max(globalStructure$G))){
                tmp = matrix(1,size(xt,1),size(xt,2))*k
                inters = (globalStructure$G==tmp)
                if((sum(sum(inters))!=0) ){#If we have a sample defined here (accord. to G)
                  returnArgs <-  find(inters) 
                  col <- returnArgs[[1]]
                  row <- returnArgs[[2]]
                  min_val = minxt[col(1),row(1)]
                  max_val = maxxt[col(1),row(1)]
                  if((min_val!=max_val)){
                    if((stepsize==0)){
                      tt=linspace(min_val,max_val,numpoints)
                    } else {
                      tt=linspace(min_val,max_val,floor((max_val-min_val)/stepsize)+1)
                    }
                    
                    if(globalStructure$bEANoReplicates){
                      #remove already taken points
                      tt=setdiff(tt,unique(xt_current[unique(col),,drop=F]))
                    }
                    
                    ntt = length(tt)
                    for(j in 1:ntt){
                      ti = tt[j]
                      xt_tmp = xt_current*(inters==0)+ti*(inters!=0)
                      
                      if((globalStructure$d_switch)){
                        if((iParallelN==1)){
                          returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
                          nfmf <- returnArgs[[1]]
                          globalStructure <- returnArgs[[2]]
                          nofv=ofv_fim(nfmf,globalStructure)
                        } else {
                          if((p==1)){
                            designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,-1,0)
                          } else {
                            nfmf = designout[[des_num]]$FIM
                            nofv = designout[[des_num]]$ofv
                            des_num = des_num +1
                          }
                        }
                      } else {
                        returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
                        nfmf <- returnArgs[[1]]
                        nofv <- returnArgs[[2]]
                        globalStructure <- returnArgs[[3]]
                      }
                      if((iParallelN==1 || (p==2 && des_num > 1))){
                        delta = nofv/optofv -1 
                        
                        if((delta>delta_max)){
                          i_index = k
                          t_max = ti
                          type = 1
                          delta_max = delta
                          dmf = nofv
                          best_it = best_it + 1
                          dmfvector[best_it]=dmf
                          fmf = nfmf
                          graph_det(globalStructure$bShowGraphs,globalStructure$optsw,xt_tmp,a_current,x_current,ni,m,dmfvector,best_it,globalStructure$Engine)
                          if((!isempty(globalStructure$strIterationFileName))){
                            write_iterationfile('Modified Exchange Algorithm',best_it,xt_tmp,a_current,x_current,ni,globalStructure$groupsize,fmf,dmf,globalStructure)
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            for(ct in 1:ni[i] ){#For all samples in group
              if((previous_type!=1 || previous_index_i!=i || previous_index_ct!=ct)){
                if((minxt[i,ct]!=maxxt[i,ct])){
                  if((stepsize==0)){
                    tt=linspace(minxt[i,ct],maxxt[i,ct],numpoints)
                  } else {
                    tt=linspace(minxt[i,ct],maxxt[i,ct],floor((maxxt[i,ct]-minxt[i,ct])/stepsize))
                  }
                  ntt = length(tt)
                  xt_tmp = xt_current
                  for(j in 1:ntt){
                    ti = tt[j]
                    xt_tmp[i,ct] = ti
                    
                    if((globalStructure$d_switch)){
                      if((iParallelN==1)){
                        returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
                        nfmf <- returnArgs[[1]]
                        globalStructure <- returnArgs[[2]]
                        nofv=ofv_fim(nfmf,globalStructure)
                      } else {
                        if((p==1)){
                          designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,-1,0)
                        } else {
                          
                          nfmf = designout[[des_num]]$FIM
                          nofv = designout[[des_num]]$ofv
                          des_num = des_num +1
                        }
                      }
                    } else {
                      returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
                      nfmf <- returnArgs[[1]]
                      nofv <- returnArgs[[2]]
                      globalStructure <- returnArgs[[3]]
                    }
                    if((iParallelN==1 || (p==2 && des_num > 1))){
                      delta = nofv/optofv -1 
                      
                      if((delta>delta_max)){
                        i_index = i
                        ct_index = ct
                        t_max = ti
                        type = 1
                        delta_max = delta
                        dmf = nofv
                        fmf = nfmf
                        best_it = best_it+1
                        dmfvector[best_it]=dmf
                        graph_det(globalStructure$bShowGraphs,globalStructure$optsw,xt_tmp,a_current,x_current,ni,m,dmfvector,best_it,globalStructure$Engine)
                        if((!isempty(globalStructure$strIterationFileName))){
                          write_iterationfile('Modified Exchange Algorithm',best_it,xt_tmp,a_current,x_current,ni,globalStructure$groupsize,fmf,dmf,globalStructure)
                        }
                      }
                      xt_tmp[i,ct] = xt_current[i,ct]
                    }
                  }
                }
              }
            }
          }
        }
        if((iParallelN == 2)){
          designsin = cell(1,0)
        }
        
      }
      
      if((globalStructure$optsw[3])){
        for(p  in 1:iParallelN){
          if((p==2)){
            stop("Parallel execution not yet implemented in PopED for R")
            #designout = execute_parallel(designsin,globalStructure)
            des_num = 1
            if((it==0)){
              nfmf = designout[[des_num]]$FIM
              optofv = ofv_fim(nfmf,globalStructure)
              des_num = des_num +1
              if((it==0 && (optofv<=0 || isnan(optofv)))){
                fprintf('The OFV is negative or zero, new initial values might be needed. OFV : %g\n',optofv)
              }
            }
          }
          
          if((globalStructure$bUseGrouped_x)){
            if((bFirstInGroup_x)){
              bFirstInGroup_x = FALSE
              for(k in 1:max(max(globalStructure$Gx))){
                tmp = matrix(1,size(x,1),size(x,2))*k
                inters = (globalStructure$Gx==tmp)
                if((sum(sum(inters))!=0) ){#If we have a discrete variabel defined here (accord. to Gx)
                  tt=get_discrete_val(inters,globalStructure$discrete_x)
                  if((length(tt)>1)){
                    ntt = length(tt)
                    for(j in 1:ntt){
                      ti = tt[j]
                      x_tmp = x_current*(inters==0)+ti*(inters!=0)
                      
                      if((globalStructure$d_switch)){
                        if((iParallelN==1)){
                          returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_tmp,a_current,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
                          nfmf <- returnArgs[[1]]
                          globalStructure <- returnArgs[[2]]
                          nofv=ofv_fim(nfmf,globalStructure)
                        } else { #PARALLEL DESIGN STORAGE
                          if((p==1)){
                            designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_current,x_tmp,a_current,-1,0)
                          } else {
                            nfmf = designout[[des_num]]$FIM
                            nofv = designout[[des_num]]$ofv
                            des_num = des_num +1
                          }
                        }
                      } else {
                        returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_tmp,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
                        nfmf <- returnArgs[[1]]
                        nofv <- returnArgs[[2]]
                        globalStructure <- returnArgs[[3]]
                      }
                      if((iParallelN==1 || (p==2 && des_num > 1))){
                        delta = nofv/optofv -1 
                        
                        if((delta>delta_max)){
                          i_index = k
                          t_max = ti
                          type = 3
                          delta_max = delta
                          dmf = nofv
                          best_it = best_it+1
                          dmfvector[best_it]=dmf
                          fmf = nfmf
                          graph_det(globalStructure$bShowGraphs,globalStructure$optsw,xt_current,a_current,x_tmp,ni,m,dmfvector,best_it,globalStructure$Engine)
                          if((!isempty(globalStructure$strIterationFileName))){
                            write_iterationfile('Modified Exchange Algorithm',best_it,xt_current,a_current,x_tmp,ni,globalStructure$groupsize,fmf,dmf,globalStructure)
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            for(ct in 1:length(x[i,]) ){#For all discrete variables in group
              if((previous_type!=3 || previous_index_i!=i || previous_index_ct!=ct)){
                discrete_val = globalStructure$discrete_x[[i,ct]]
                if((length(discrete_val)>1)){
                  ntt = length(discrete_val)
                  x_tmp = x_current
                  for(j in 1:ntt){
                    ti = discrete_val[j]
                    x_tmp[i,ct] = ti
                    
                    if((globalStructure$d_switch)){
                      if((iParallelN==1)){
                        returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_tmp,a_current,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
                        nfmf <- returnArgs[[1]]
                        globalStructure <- returnArgs[[2]]
                        nofv=ofv_fim(nfmf,globalStructure)
                      } else { #PARALLEL DESIGN STORAGE
                        if((p==1)){
                          designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_current,x_tmp,a_current,-1,0)
                        } else {
                          nfmf = designout[[des_num]]$FIM
                          nofv = designout[[des_num]]$ofv
                          des_num = des_num +1
                        }
                      }
                    } else {
                      returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_tmp,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
                      nfmf <- returnArgs[[1]]
                      nofv <- returnArgs[[2]]
                      globalStructure <- returnArgs[[3]]
                    }
                    if((iParallelN==1 || (p==2 && des_num > 1))){
                      delta = nofv/optofv -1 
                      
                      if((delta>delta_max)){
                        i_index = i
                        ct_index = ct
                        t_max = ti
                        type = 3
                        delta_max = delta
                        dmf = nofv
                        best_it = best_it+1
                        dmfvector[best_it]=dmf
                        fmf = nfmf
                        graph_det(globalStructure$bShowGraphs,globalStructure$optsw,xt_current,a_current,x_tmp,ni,m,dmfvector,best_it,globalStructure$Engine)
                        if((!isempty(globalStructure$strIterationFileName))){
                          write_iterationfile('Modified Exchange Algorithm',best_it,xt_current,a_current,x_tmp,ni,globalStructure$groupsize,fmf,dmf,globalStructure)
                        }
                      }
                    }
                    x_tmp[i,ct] = x_current[i,ct]
                  }
                }
              }
            }
          }
        }
        if((iParallelN == 2)){
          designsin = cell(1,0)
        }
      }
      
      if((globalStructure$optsw[4]) ){#Optimize over covariates
        for(p  in 1:iParallelN){
          if((p==2)){
            stop("Parallel execution not yet implemented in PopED for R")
            #designout = execute_parallel(designsin,globalStructure)
            des_num = 1
            if((it==0)){
              nfmf = designout[[des_num]]$FIM
              optofv = ofv_fim(nfmf,globalStructure)
              des_num = des_num +1
              if((it==0 && (optofv<=0 || isnan(optofv)))){
                fprintf('The OFV is negative or zero, new initial values might be needed. OFV : %g\n',optofv)
              }
            }
          }
          if((globalStructure$bUseGrouped_a)){
            if((bFirstInGroup_a)){
              bFirstInGroup_a = FALSE
              for(k in 1:max(max(globalStructure$Ga))){
                tmp = matrix(1,size(a,1),size(a,2))*k
                inters = (globalStructure$Ga==tmp)
                if((sum(sum(inters))!=0) ){#If we have a covariate defined here (accord. to Ga)
                  returnArgs <-  find(inters) 
                  col <- returnArgs[[1]]
                  row <- returnArgs[[2]]
                  min_val = mina[col(1),row(1)]
                  max_val = maxa[col(1),row(1)]
                  if((min_val!=max_val)){
                    if((stepsize==0)){
                      tt=linspace(min_val,max_val,numpoints)
                    } else {
                      tt=linspace(min_val,max_val,floor((max_val-min_val)/stepsize))
                    }
                    ntt = length(tt)
                    for(j in 1:ntt){
                      ti = tt[j]
                      a_tmp = a_current*(inters==0)+ti*(inters!=0)
                      
                      if((globalStructure$d_switch)){
                        if((iParallelN==1)){
                          returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_current,a_tmp,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
                          nfmf <- returnArgs[[1]]
                          globalStructure <- returnArgs[[2]]
                          nofv=ofv_fim(nfmf,globalStructure)
                        } else { #PARALLEL DESIGN STORAGE
                          if((p==1)){
                            designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_current,x_current,a_tmp,-1,0)
                          } else {
                            nfmf = designout[[des_num]]$FIM
                            nofv = designout[[des_num]]$ofv
                            des_num = des_num +1
                          }
                        }
                      } else {
                        returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
                        nfmf <- returnArgs[[1]]
                        nofv <- returnArgs[[2]]
                        globalStructure <- returnArgs[[3]]
                      }
                      if((iParallelN==1 || (p==2 && des_num > 1))){
                        delta = nofv/optofv -1
                        if((delta>delta_max)){
                          i_index = k
                          t_max = ti
                          type = 2
                          delta_max = delta
                          dmf = nofv
                          best_it = best_it+1
                          dmfvector[best_it]=dmf
                          fmf = nfmf
                          graph_det(globalStructure$bShowGraphs,globalStructure$optsw,xt_current,a_tmp,x_current,ni,m,dmfvector,best_it,globalStructure$Engine)
                          if((!isempty(globalStructure$strIterationFileName))){
                            write_iterationfile('Modified Exchange Algorithm',best_it,xt_current,a_tmp,x_current,ni,globalStructure$groupsize,fmf,dmf,globalStructure)
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            for(ct in 1:length(a[i,]) ){#For all covariates in group
              if((previous_type!=2 || previous_index_i!=i || previous_index_ct!=ct)){
                if((mina[i,ct]!=maxa[i,ct])){
                  if((stepsize==0)){
                    tt=linspace(mina[i,ct],maxa[i,ct],numpoints)
                  } else {
                    tt=linspace(mina[i,ct],maxa[i,ct],floor((maxa[i,ct]-mina[i,ct])/stepsize))
                  }
                  ntt = length(tt)
                  a_tmp = a_current
                  for(j in 1:ntt){
                    ti = tt[j]
                    a_tmp[i,ct] = ti
                    
                    if((globalStructure$d_switch)){
                      if((iParallelN==1)){
                        returnArgs <-  mftot(model_switch,globalStructure$groupsize,ni,xt_current,x_current,a_tmp,bpop,d,globalStructure$sigma,doccfull,globalStructure) 
                        nfmf <- returnArgs[[1]]
                        globalStructure <- returnArgs[[2]]
                        nofv=ofv_fim(nfmf,globalStructure)
                      } else { #PARALLEL DESIGN STORAGE
                        if((p==1)){
                          designsin = update_designinlist(designsin,globalStructure$groupsize,ni,xt_current,x_current,a_tmp,-1,0)
                        } else {
                          nfmf = designout[[des_num]]$FIM
                          nofv = designout[[des_num]]$ofv
                          des_num = des_num +1
                        }
                      }
                    } else {
                      returnArgs <-  ed_mftot(model_switch,globalStructure$groupsize,ni,xt_tmp,x_current,a_current,bpopdescr,ddescr,globalStructure$covd,globalStructure$sigma,globalStructure$docc,globalStructure) 
                      nfmf <- returnArgs[[1]]
                      nofv <- returnArgs[[2]]
                      globalStructure <- returnArgs[[3]]
                    }
                    if((iParallelN==1 || (p==2 && des_num > 1))){
                      delta = nofv/optofv -1 
                      if((delta>delta_max)){
                        i_index = i
                        ct_index = ct
                        t_max = ti
                        type = 2
                        delta_max = delta
                        dmf = nofv
                        best_it = best_it+1
                        dmfvector[best_it]=dmf
                        fmf = nfmf
                        graph_det(globalStructure$bShowGraphs,globalStructure$optsw,xt_current,a_tmp,x_current,ni,m,dmfvector,best_it,globalStructure$Engine)
                        if((!isempty(globalStructure$strIterationFileName))){
                          write_iterationfile('Modified Exchange Algorithm',best_it,xt_current,a_tmp,x_current,ni,globalStructure$groupsize,fmf,dmf,globalStructure)
                        }
                      }
                    }
                    a_tmp[i,ct] = a_current[i,ct]
                  }
                }
              }
            }
          }
        }
        if((iParallelN == 2)){
          designsin = cell(1,0)
        }
      }
    }
    it = it+1
    fprintf('FE - It. : %d\n',it)
    
    if((delta_max>rho)){
      previous_index_i = i_index
      previous_type = type
      
      if((type==1)){
        if((!globalStructure$bUseGrouped_xt)){
          fprintf('Exchanged sample %d in group/ind %d from %g to %g\n',ct_index,i_index,xt_current[i_index,ct_index],t_max)
          previous_index_ct = ct_index
          xt_current[i_index,ct_index] = t_max
        } else {
          fprintf('Exchanged a grouped sample (grouped as %d) to %g\n',i_index,t_max)
          tmp = matrix(1,size(xt,1),size(xt,2))*i_index
          inters = (globalStructure$G==tmp)
          xt_current=xt_current*(inters==0)+t_max*(inters!=0)
        }
      } else {
        if((type==2)){
          if((!globalStructure$bUseGrouped_a)){
            fprintf('Exchanged covariate %d in group/ind %d from %g to %g\n',ct_index,i_index,a_current[i_index,ct_index],t_max)
            previous_index_ct = ct_index
            a_current[i_index,ct_index] = t_max
          } else {
            fprintf('Exchanged a grouped covariate (grouped as %d) to %g\n',i_index,t_max)
            tmp = matrix(1,size(a,1),size(a,2))*i_index
            inters = (globalStructure$Ga==tmp)
            a_current=a_current*(inters==0)+t_max*(inters!=0)
          }
        } else {
          if((type==3)){
            if((!globalStructure$bUseGrouped_x)){
              fprintf('Exchanged discrete variable %d in group/ind %d from %g to %g\n',ct_index,i_index,x_current[i_index,ct_index],t_max)
              previous_index_ct = ct_index
              x_current[i_index,ct_index] = t_max
            } else {
              fprintf('Exchanged a grouped discrete variable (grouped as %d) to %g\n',i_index,t_max)
              tmp = matrix(1,size(x,1),size(x,2))*i_index
              inters = (globalStructure$Gx==tmp)
              x_current=x_current*(inters==0)+t_max*(inters!=0)
            }
          }
        }
      }
      optofv=dmf
    }
    fprintf('Delta : %g   OFV. : %g\n',delta_max,dmf)
  }
  xt = xt_current
  a = a_current
  x = x_current
  
  return(list( xt= xt,x=x,a=a,fmf=fmf,dmf=dmf,globalStructure=globalStructure)) 
}

graph_det <- function(bShowGraphs,optsw,xtopt,aopt,xopt,gni,m,dmfvector,it,Engine){
  ## screen output
  # if((bShowGraphs)){
  #     if((Engine$Type==1)){
  #         set(0,'CurrentFigure',1)
  #     } else {
  #         figure(1)
  #     }
  #     clf
  #     optSum = optsw[2]+optsw[3]+optsw[4]
  #     numRows = 1
  #     if((optSum>1)){
  #         numRows = 2
  #     }
  # 
  #     if((optsw[2])){
  #         subplot(numRows,2,1)
  #         title('The current sampling times for each group')
  #         xlabel('Sampling time')
  #         ylabel('Group nr')
  #         ##hold on
  #         for(i in 1:m){
  #             plot(xtopt[i,1:gni(i)],matrix(1,1,gni(i))*i,'b*')
  #         }
  #         ##hold off
  #     }
  #     if(optsw[3]){
  #         subplot(numRows,2,1+optsw[2])
  #         title('The current discrete var. for each group')
  #         xlabel('Design var.-value')
  #         ylabel('Group nr')
  #         ##hold on
  #         for(i in 1:m){
  #             plot(xopt(i,),matrix(1,1,size(xopt,2))*i,'b*')
  #         }
  #         ##hold off
  #     }
  #     if(optsw[4]){
  #         subplot(numRows,2,1+optsw[2]+optsw[3])
  #         title('The current covariates for each group')
  #         xlabel('Covariate-value')
  #         ylabel('Group nr')
  #         ##hold on
  #         for(i in 1:m){
  #             plot(aopt[i,],matrix(1,1,size(aopt,2))*i,'b*')
  #         }
  #         ##hold off
  #     }
  # 
  #     subplot(numRows,2,1+optSum)
  #     ##hold on
  #     title('OFV(FIM)')
  #     xlabel('Iteration')
  #     ylabel('Value')
  #     plot(1:it,dmfvector,'-b')
  #     ##hold off
  #     drawnow
  # }
  
}

#Get the discrete values for an intersection matrix inters
get_discrete_val <- function(inters,discrete_x){
  for(i in 1:size(inters,1)){
    for(j in 1:size(inters,2)){
      if((inters[i,j]==1)){
        discrete_val=discrete_x[[i,j]]
        return
      }
    }
  }
  return( discrete_val ) 
}
