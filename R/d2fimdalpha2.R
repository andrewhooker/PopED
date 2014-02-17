d2fimdalpha2 <- function(alpha, model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,globalStructure,ha){

bpop=bpopdescr[,2,drop=F]
bpop[bpopdescr[,1,drop=F]!=0]=alpha[1:sum(bpopdescr[,1,drop=F]!=0)]
d=ddescr[,2,drop=F]
d[ddescr[,1]!=0]=alpha[sum(bpopdescr[,1,drop=F]!=0)+1:end]
d=getfulld(d,covd)

fim=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,globalStructure)

i=1
hess=array(0,dim=c(length(fim),length(fim),length(alpha),length(alpha)))
for(i in 1:length(alpha)){
    alpha_plus=alpha
    alpha_plus[i]=alpha[i]+ha
    bpop=bpopdescr[,2,drop=F]
    bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
    d=ddescr[,2,drop=F]
    d[ddescr[,1]!=0]=alpha_plus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
    d=getfulld(d,covd)
    fim_plus=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,globalStructure)
    
    for(j in 1:i){
        alpha_plus2=alpha
        alpha_plus2[j]=alpha[j]+ha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus2[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_plus2[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        fim_plus2=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,globalStructure)
        
        alpha_plus_plus=alpha
        alpha_plus_plus[i]=alpha_plus_plus[i]+ha
        alpha_plus_plus[j]=alpha_plus_plus[j]+ha
        
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_plus_plus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        fim_plus_plus=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,globalStructure)
        
        hess[,,i,j]=(fim_plus_plus-fim_plus-fim_plus2+fim)/ha^2
        hess[,,j,i]=hess[,,i,j]
    } 
}
return(list( hess= hess,fim =fim )) 
}
