library(PopED)


# typical WT prediction from Sumpter & Holford 2011
# PMA in years
get_typical_weight <- function(PMA,SEX=2){
  WTmax1 <- 2.76
  TM50wt1 <- 38.5
  HILLwt1 <- 12.9
  HILL2wt1 <- 2.74
  PMA1 <- PMA*52
  
  WTmax2 <- 16.4
  TM50wt2 <- 2.1
  HILLwt2 <- 2.04
  HILL2wt2 <- 1
  
  WTmax3 <- 40.2
  TM50wt3 <- 12.4
  HILLwt3 <- 2.87
  HILL2wt3 <- 0
  
  TLAGwt4 <- 12.4
  WTmax4 <- 33.6
  THALFwt4 <- 3.61
  
  FFEM <- 1
  if(SEX==1) FFEM <- 0.884
  
  wt1 <- FFEM*WTmax1/(1+(TM50wt1/PMA1)^((PMA1<TM50wt1)*HILLwt1+(PMA1>=TM50wt1)*HILL2wt1))
  wt2 <- WTmax2/(1+(TM50wt2/PMA)^((PMA<TM50wt2)*HILLwt2+(PMA>=TM50wt2)*HILL2wt2))
  wt3 <- WTmax3/(1+(TM50wt3/PMA)^((PMA<TM50wt3)*HILLwt3+(PMA>=TM50wt3)*HILL2wt3))
  wt4 <- (PMA > TLAGwt4)*(FFEM*WTmax4*(1-exp(-log(2)/THALFwt4*(PMA-TLAGwt4))))
  
  return(wt1 + wt2 + wt3 + wt4)
}

# plot typical weight function
PMA = seq(0,80,by=0.1)
df1 <- data.frame(PMA = PMA, WT=get_typical_weight(PMA,SEX=2),SEX=as.factor(2))
df2 <- data.frame(PMA = PMA,WT=get_typical_weight(PMA,SEX=1),SEX=as.factor(1))
df <- rbind.data.frame(df1,df2)
levels(df$SEX) <- c("Male","Female")
library(ggplot2)
(plot1 <- ggplot(data=df,aes(x = PMA, y=WT, group=SEX)) + geom_line(aes(color=SEX)) + 
  labs(y="Weight (kg)",x="Post Menstrual Age (years)", title = "Prediction of typical weight from PMA"))


PK.1.comp.maturation.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    
    WT <- get_typical_weight(PMA,SEX=SEX)
    CL=CL*(WT/70)^(3/4)*(PMA^HILL)/(TM50^HILL+PMA^HILL)
    V=V*(WT/70)
    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*xt) 
    
    return(list( y= y,poped.db=poped.db))
  })
}


PK.1.comp.maturation.fg <- function(x,a,bpop,b,bocc){
  parameters=c( CL=bpop[1]*exp(b[1]),
                V=bpop[2]*exp(b[2]),
                TM50=bpop[3]*exp(b[3]),
                HILL=bpop[4],
                PMA=a[1],
                SEX=a[2])
  return( parameters ) 
}

poped.db <- create.poped.database(ff_file="PK.1.comp.maturation.ff",
                                  fError_file="feps.add.prop",
                                  fg_file="PK.1.comp.maturation.fg",
                                  groupsize=rbind(50,20,20,20),
                                  m=4,
                                  sigma=c(0.015,0.0015),
                                  notfixed_sigma = c(1,0),
                                  bpop=c(CL=3.8,V=20,TM50=60/52,HILL=3), 
                                  d=c(CL=0.05,V=0.05,TM50=0.05), 
                                  xt=c( 1,2,4,6,8,24),
                                  minxt=0,
                                  maxxt=24,
                                  bUseGrouped_xt=1,
                                  a=list(c(PMA=25,SEX=2),
                                         c(PMA=15,SEX=2),
                                         c(PMA=10,SEX=2),
                                         c(PMA=5,SEX=2)),
                                  maxa=list(c(PMA=30,SEX=2)),
                                  mina=list(c(PMA=1,SEX=2)))



##  create plot of model 
plot_model_prediction(poped.db)
plot_model_prediction(poped.db,IPRED=T,DV=T,separate.groups=T)


## evaluate initial design
evaluate_design(poped.db)
shrinkage(poped.db)


# Optimization of sample times and WT
output <- poped_optim(poped.db,opt_xt=T,opt_a=T)

summary(output)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db)




