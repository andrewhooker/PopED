## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

ed_mftot <- function(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,globalStructure){
  #+++++++++++++++++++++++ ED OFV(MF) VALUE
  s=0
  s1=0
  
  fim_list=cell(1,globalStructure$ED_samp_size)
  d_gen_list=cell(1,globalStructure$ED_samp_size)
  docc_gen_list=cell(1,globalStructure$ED_samp_size)
  
  
  bpop_gen = pargen(bpopdescr,globalStructure$user_distribution_pointer,globalStructure$ED_samp_size,globalStructure$bLHS,zeros(1,0),globalStructure)
  
  for(ct in 1:globalStructure$ED_samp_size){
    d_gen = getfulld(pargen(ddescr,globalStructure$user_distribution_pointer,1,globalStructure$bLHS,ct,globalStructure),covd)
    docc_gen = getfulld(pargen(docc,globalStructure$user_distribution_pointer,1,globalStructure$bLHS,ct,globalStructure),globalStructure$covdocc)
    returnArgs <- mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop_gen[ct,],d_gen,sigma,docc_gen,globalStructure) 
    mftmp <- returnArgs[[1]]
    globalStructure <- returnArgs[[2]]
    s=s+ofv_fim(mftmp,globalStructure)
    s1=s1+mftmp
    fim_list[[ct]]=mftmp
    d_gen_list[[ct]]=d_gen
    docc_gen_list[[ct]]=docc_gen
  }
  if((!isempty(globalStructure$ed_penalty_pointer))){
    returnArgs <- feval(globalStructure$ed_penalty_pointer,fim_list,bpop_gen,d_gen_list,docc_gen_list,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,globalStructure) 
    ED_fim <- returnArgs[[1]]
    ED_ofv <- returnArgs[[2]]
    globalStructure <- returnArgs[[3]]
  } else {
    ED_ofv=s/globalStructure$ED_samp_size
    ED_fim=s1/globalStructure$ED_samp_size
  }
  return(list( ED_fim= ED_fim,ED_ofv=ED_ofv,globalStructure=globalStructure)) 
}

