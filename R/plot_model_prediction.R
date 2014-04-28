#' Model predictions 
#' 
#' Function generates model predictions for the typical value in the population,
#' individual predictions and data predictions.
#' 
#' @inheritParams RS_opt
#' @param models_to_use Which model number should we use?
#' @param model_num_points How many points should be plotted.  If not a number then the design in poped.db is used.
#' @param model_minxt The minimum of the sample times for the predictions.
#' @param model_maxxt The maximum of the sample times for the predictions.
#' @param include_sample_times Should the sample times from poped.db be included in the predictions?
#' @param IPRED Should we simulate individual predictions?
#' @param DV should we simulate observations?
#' @param num_ids The number of individuals to simulate if using IPRED or DV.
#' @param groups_to_use Which groups should we use for predictions from the poped.db.
#' @return A dataframe of simulated data, either with some dense grid of samples or based on the design in the poped
#' database.
#' 
#' @family evaluate_design
#' @family Simulation
#' 
#' @example tests/testthat/examples_fcn_doc/examples_model_prediction.R
#' 

#' 
## allow for input not from poped.db
model_prediction <- function(poped.db,
                             models_to_use="all",
                             model_num_points=NULL,
                             model_minxt=NULL,model_maxxt=NULL,
                             include_sample_times=T,
                             groups_to_use="all",
                             IPRED=FALSE,
                             DV=FALSE,
                             num_ids=100){
  
  if(DV) IPRED=T
  
  design = poped.db$design
  design$xt <- poped.db$gxt
  design$x <- poped.db$gx
  design$a <- poped.db$ga
  
  docc_size = 0
  if((!isempty(design$docc[,2]))){
    docc_size = size(design$docc[,2,drop=F],1)
  }
  d_size = 0
  if((!isempty(design$d[,2]))){
    d_size = size(design$d[,2,drop=F],1)
  }
  
  if(IPRED){
    #d_dist = cbind(rep(1,d_size),rep(0,d_size), design$d[,2])
    #docc_dist = cbind(rep(1,docc_size),rep(0,docc_size), design$docc[,2])
    fulld = getfulld(design$d[,2],design$covd)
    fulldocc = getfulld(design$docc[,2,drop=F],design$covdocc)
    
    #Sample all the b's and bocc's (same samples for all sub models)
    b_sim_matrix = zeros(num_ids,length(design$d[,2]))
    bocc_sim_matrix = zeros(num_ids*poped.db$NumOcc,length(design$docc[,2,drop=F]))
    
    #     if((!globalStructure$d_switch) ){#ED-optimal
    #       d_dist=pargen(design$d,globalStructure$user_distribution_pointer,max(max(model_num_simulations),iMaxCorrIndNeeded),globalStructure$bLHS,zeros(0,1),globalStructure)
    #       for(j in 1:max(max(model_num_simulations),iMaxCorrIndNeeded)){
    #         tmp_d_dist = matrix(c(matrix(1,length(design$d(,2)),1),t(zeros(length(design$d(,2)),1),d_dist(j,))),nrow=1,byrow=T)
    #         b_sim_matrix(j,) = pargen(tmp_d_dist,globalStructure$user_distribution_pointer,1,globalStructure$bLHS,j,globalStructure)
    #       }
    #     } else { #D-optimal design
    b_sim_matrix = rmvnorm(num_ids,sigma=fulld)
    #b_sim_matrix = pargen(d_dist,poped.db$user_distribution_pointer,num_ids,poped.db$bLHS,zeros(0,1),poped.db)
    #var(b_sim_matrix)
    #d_dist
    #}
    #b_sim_matrix = correlateSamples(t(b_sim_matrix),fulld)
    #}
    
    if(nrow(fulldocc)!=0){
      #       if((!globalStructure$d_switch) ){#ED-optimal
      #         docc_dist=pargen(design$docc,globalStructure$user_distribution_pointer,max(max(model_num_simulations)*popedInput$NumOcc,iMaxCorrIndNeeded),globalStructure$bLHS,zeros(0,1),globalStructure)
      #         for(j in 1:max(max(model_num_simulations)*popedInput$NumOcc,iMaxCorrIndNeeded)){
      #           tmp_docc_dist = matrix(c(matrix(1,length(design$docc[,2,drop=F]),1),t(zeros(length(design$docc[,2,drop=F]),1),docc_dist(j,))),nrow=1,byrow=T)
      #           bocc_sim_matrix(j,) = pargen(tmp_docc_dist,globalStructure$user_distribution_pointer,1,globalStructure$bLHS,j,globalStructure)
      #         }
      #       } else { #D-optimal design
      bocc_sim_matrix = rmvnorm(num_ids*poped.db$NumOcc,sigma=fulldocc)
      #bocc_sim_matrix = pargen(docc_dist,globalStructure$user_distribution_pointer,max(max(model_num_simulations)*popedInput$NumOcc,iMaxCorrIndNeeded),globalStructure$bLHS,zeros(0,1),globalStructure)
      #}
      #Add correlation of bocc, might not work
      #bocc_sim_matrix = correlateSamples(t(bocc_sim_matrix),fulldocc)
    }
    
  }
  
  if(all(groups_to_use=="all")){
    groups_to_use = 1:size(design$xt,1)
  }
  if(all(models_to_use=="all")){
    models_to_use = unique(as.vector(design$model_switch))
  }
  
  df <- data.frame()
  
  for(i in 1:length(groups_to_use)){
    if((isempty(design$a))){
      a_i = zeros(0,1)
    } else {
      a_i = design$a[groups_to_use[i],,drop=F]
    }
    if((isempty(design$x))){
      x_i = zeros(0,1)
    } else {
      x_i = design$x[groups_to_use[i],,drop=F]
    }
    
    if(all(is.null(model_num_points))){
      xt_i = design$xt[groups_to_use[i],1:design$ni[groups_to_use[i]]]
      model_switch_i = design$model_switch[groups_to_use[i],1:design$ni[groups_to_use[i]]]
      if(!all(models_to_use == unique(as.vector(design$model_switch)))){ ## needs testing
        xt_i = xt_i[model_switch_i %in% models_to_use]
        model_switch_i = model_switch_i[model_switch_i %in% models_to_use]
      }
    } else {
      xt_i <- c()
      model_switch_i <- c()
      if(length(models_to_use)>1 && length(model_num_points)==1) model_num_points <- rep(model_num_points,length(models_to_use))
      for(j in models_to_use){
        if(is.null(model_minxt)){
          minv <- min(as.vector(design$minxt[design$model_switch==j])) 
        } else {                    
          minv = model_minxt[j]
        }
        if(is.null(model_maxxt)){
          maxv <- max(as.vector(design$maxxt[design$model_switch==j])) 
        } else {
          maxv = model_maxxt[j]
        }                #xt = t(seq(minv,maxv,length.out=model_num_points[i]))
        
        xt_i= c(xt_i,seq(minv,maxv,length.out=model_num_points[j]))
        
        #model.pred <- rbind(xt)
        #model.pred <- data.frame(Time=xt)
        #model.pred <- c(model.pred,foo=xt)
        #browser()
        model_switch_i = c(model_switch_i,j*matrix(1,1,model_num_points[j]))
      }
      if(include_sample_times){
        xt_i_extra = design$xt[groups_to_use[i],1:design$ni[groups_to_use[i]]]
        model_switch_i_extra = design$model_switch[groups_to_use[i],1:design$ni[groups_to_use[i]]]
        if(!all(models_to_use == unique(as.vector(design$model_switch)))){ ## needs testing
          xt_i_extra = xt_i_extra[model_switch_i_extra %in% models_to_use]
          model_switch_i_extra = model_switch_i_extra[model_switch_i_extra %in% models_to_use]
        }
        tmp.include <- !(xt_i_extra %in% xt_i)
        xt_i <- c(xt_i,xt_i_extra[tmp.include])
        model_switch_i <- c(model_switch_i,model_switch_i_extra[tmp.include])
        tmp.order <- order(xt_i)
        xt_i <- xt_i[tmp.order]
        model_switch_i <- model_switch_i[tmp.order]
      }
    }
    g0 = feval(poped.db$fg_pointer,x_i,a_i,design$bpop[,2,drop=F],zeros(1,d_size),zeros(docc_size,poped.db$NumOcc))
    
    pred <- feval(poped.db$ff_pointer,model_switch_i,xt_i,g0,poped.db)
    pred <- drop(pred[[1]])
        
    group.df <- data.frame(Time=xt_i,PRED=pred,Group=groups_to_use[i],Model=model_switch_i)
    #     group.df <- data.frame(Time=xt_i,PRED=drop(pred[[1]]),Group=groups_to_use[i],
    #                            ##paste("Group",i,sep="_"),
    #                            Model=model_switch_i)
    #     
    
    if(IPRED){
      group.df.ipred <- data.frame()
      bocc_start= 1
      for(j in 1:num_ids){
        tmp.df <- group.df
        bocc_stop=bocc_start + poped.db$NumOcc - 1
        if(nrow(fulldocc)==0){ 
          bocc_start=0
          bocc_stop=0
        }
        fg_sim = feval(poped.db$fg_pointer,x_i,a_i,design$bpop[,2,drop=F],b_sim_matrix[j,],t(bocc_sim_matrix[bocc_start:bocc_stop,]))
        bocc_start = bocc_stop + 1
        ipred <- feval(poped.db$ff_pointer,model_switch_i,xt_i,fg_sim,poped.db)
        ipred <- drop(ipred[[1]])
        ID <- (i-1)*num_ids+j
        tmp.df["ID"] <- ID
        tmp.df["IPRED"] <- ipred
        
        if(DV){
          eps_sim = rmvnorm(length(xt_i),sigma=design$sigma)
          dv <- feval(poped.db$ferror_pointer,model_switch_i,xt_i,fg_sim,eps_sim,poped.db) 
          dv <- drop(dv[[1]])
          tmp.df["DV"] <- dv
        }
        
        group.df.ipred <- rbind(group.df.ipred,tmp.df)
      }
      group.df <- group.df.ipred
    }
    
    df <- rbind(df,group.df)
    #model.pred  <- rbind(model.pred,i=returnArgs[[1]])
    #model.pred[paste("Group",i,sep="_")]  <- returnArgs[[1]]
  }
  df$Group <- as.factor(df$Group)
  df$Model <- as.factor(df$Model)
  if(IPRED) df$ID <- as.factor(df$ID)
  return( df ) 
}

#' Plot model predictions 
#' 
#' Function plots model predictions for the typical value in the population,
#' individual predictions and data predictions.
#' 
#' @inheritParams RS_opt
#' @inheritParams model_prediction
#' @param separate.groups Should there be separate plots for each group.
#' @param sample.times Should sample times be shown on the plots.
#' @param sample.times.IPRED Should sample times be shown based on the IPRED y-values.
#' @param sample.times.DV Should sample times be shown based on the DV y-values.
#' @param PRED Should a PRED line be drawn.
#' @param IPRED.lines Should IPRED lines be drawn?
#' @param alpha.IPRED.lines What should the transparency for the IPRED.lines be?
#' @param alpha.IPRED What should the tranparency of the IPRED CI?
#' @param sample.times.size What should the size of the sample.times be?
#' @param alpha.DV What should the tranparency of the DV CI?
#' @param DV.lines Should DV lines be drawn?
#' @param DV.points Should DV points be drawn?
#' @param alpha.DV.lines What should the transparency for the DV.lines be?
#' @param alpha.DV.points What should the transparency for the DV.points be?
#' @param sample.times.DV.points TRUE or FALSE.
#' @param sample.times.DV.lines TRUE or FALSE.
#' @param alpha.sample.times.DV.points What should the transparency for the sample.times.DV.points be?
#' @param alpha.sample.times.DV.lines What should the transparency for the sample.times.DV.lines be?
#' @param y_lab The label of the y-axis.
#' @param facet_scales Can be "free", "fixed", "free_x" or "free_y"
#' @param facet_label_names TRUE or FALSE
#' 
#' 
#' @return A \link[ggplot2]{ggplot2} object.
#' 
#' @family evaluate_design
#' @family Simulation
#' @family Graphics
#' 
#' @example tests/testthat/examples_fcn_doc/examples_plot_model_prediction.R
#' 

plot_model_prediction <- function(poped.db,
                                  ##models_to_use="all",
                                  model_num_points=100,
                                  ##model_minxt,
                                  ##model_maxxt,
                                  ##groups_to_plot,
                                  ##bShowGraph,bPlotInSame,
                                  separate.groups=F, 
                                  sample.times=T, 
                                  sample.times.IPRED=F,
                                  sample.times.DV=F,
                                  PRED=T,
                                  IPRED=F,
                                  IPRED.lines=F,
                                  alpha.IPRED.lines=0.1,
                                  alpha.IPRED=0.3,
                                  sample.times.size=4,
                                  DV=F,
                                  alpha.DV=0.3,
                                  DV.lines=F,
                                  DV.points=F,
                                  alpha.DV.lines=0.3,
                                  alpha.DV.points=0.3,
                                  sample.times.DV.points=F,
                                  sample.times.DV.lines=F,
                                  alpha.sample.times.DV.points=0.3,
                                  alpha.sample.times.DV.lines=0.3,
                                  y_lab="Model Predictions",
                                  facet_scales="fixed", # could be "free", "fixed", "free_x" or "free_y"
                                  facet_label_names = T, 
                                  ...){
  
  df <-  model_prediction(poped.db,
                          #models_to_use,
                          model_num_points=model_num_points,
                          ##model_minxt,
                          ##model_maxxt,
                          ##groups_to_plot,
                          ...)
  
  if(sample.times){
    df.2 <-  model_prediction(poped.db,
                              ##models_to_use,
                              ##model_num_points=NULL,
                              ##model_minxt,model_maxxt,
                              ##groups_to_plot,
                              ...)
  }
  
  if(IPRED || IPRED.lines || DV){
    df.ipred <-  model_prediction(poped.db,
                                  #models_to_use,
                                  model_num_points=model_num_points,
                                  ##model_minxt,
                                  ##model_maxxt,
                                  ##groups_to_plot,
                                  IPRED=T,
                                  DV=DV,
                                  ...)
    
    if(sample.times.IPRED || sample.times.DV || sample.times.DV.points || sample.times.DV.lines){
      df.ipred.samples <- df.ipred[df.ipred$Time %in% poped.db$design$xt,]  
    }
  }
  
  ## if((bPlotModel(i))){
  ##     if((bShowGraph)){
  ##         if((!bPlotInSame)){
  ##             fi = figure
  ##         }
  ##         ##hold on
  ##         title(sprintf('Model #d',i))
  
  ##         xlabel('time')
  ##         ylabel('DV')
  ##     }
  
  ## if((bShowGraph)){
  ##     ##hold on
  ##     plot(xt,model_values(i,))
  ##     xlim(matrix(c(minv, maxv),nrow=1,byrow=T))
  ##     #ylim(matrix(c(0,max(model_values(i,))*1.1),nrow=1,byrow=T))
  ##     ##hold off
  ## }
  
  ##df <- df[[1]]
  
  ##library(reshape2)
  ##dat <- melt(df,id=c("Time"))
  ##names(dat) <- c("Time","variable","DV")
  
  many.models <- F
  if(length(unique(df$Model))>1) many.models <- T
  
  many.groups <- F
  if(length(unique(df$Group))>1) many.groups <- T
  
  #library(ggplot2)
  Group <- c()
  Time <- c()
  ID <- c()
   
  if(facet_label_names){
    if(exists("df")) levels(df$Model) <- paste("Model:",levels(df$Model))
    if(exists("df")) levels(df$Group) <- paste("Group:",levels(df$Group))
    if(exists("df.2")) levels(df.2$Model) <- paste("Model:",levels(df.2$Model))
    if(exists("df.2")) levels(df.2$Group) <- paste("Group:",levels(df.2$Group))
    if(exists("df.ipred")) levels(df.ipred$Model) <- paste("Model:",levels(df.ipred$Model))
    if(exists("df.ipred")) levels(df.ipred$Group) <- paste("Group:",levels(df.ipred$Group))
    if(exists("df.ipred.samples")) levels(df.ipred.samples$Model) <- paste("Model:",levels(df.ipred.samples$Model))
    if(exists("df.ipred.samples")) levels(df.ipred.samples$Group) <- paste("Group:",levels(df.ipred.samples$Group))
  }
    
  p <- ggplot(data=df,aes(x=Time,y=PRED)) #+ labs(colour = NULL)
  #p <- ggplot(data=df,aes(x=Time,y=PRED)) #+ labs(colour = NULL)
  if(many.groups) p <- ggplot(data=df,aes(x=Time,y=PRED,group=Group,color=Group,fill=Group)) #+ labs(colour = NULL)
  if(many.models) p <- p + facet_wrap(~Model,scales=facet_scales)
  if(separate.groups && many.models) p <- ggplot(data=df,aes(x=Time,y=PRED,group=Group)) + facet_grid(Model ~ Group,scales=facet_scales)
  if(separate.groups && !many.models) p <- ggplot(data=df,aes(x=Time,y=PRED,group=Group)) + facet_wrap(~Group,scales=facet_scales)
  if(IPRED.lines) p <- p + geom_line(aes(y=IPRED,group=ID),data=df.ipred,alpha=alpha.IPRED.lines)
  if(DV.lines) p <- p + geom_line(aes(y=DV,group=ID),data=df.ipred,alpha=alpha.DV.lines)
  if(DV.points) p <- p + geom_point(aes(y=DV,group=ID),data=df.ipred,alpha=alpha.DV.points)
  if(IPRED) p <- p + stat_summary(data=df.ipred,aes(x=Time,y=IPRED,color=NULL),geom="ribbon",fun.data="median_hilow",alpha=alpha.IPRED)
  if(DV) p <- p + stat_summary(data=df.ipred,aes(x=Time,y=DV,color=NULL),geom="ribbon",fun.data="median_hilow",alpha=alpha.DV)
  if(PRED) p <- p + geom_line()
  if(sample.times) p <- p+geom_point(data=df.2,size=sample.times.size)#,aes(x=Time,y=DV,group=Group))#,color=Group))
  #if(sample.times) p <- p+geom_point(data=df.2,size=sample.times.size,position = position_jitter(w = 0, h = 10),alpha=0.8)#,aes(x=Time,y=DV,group=Group))#,color=Group))
  #if(sample.times) p <- p+geom_jitter(data=df.2,size=sample.times.size,position = position_jitter(height = 10))#,aes(x=Time,y=DV,group=Group))#,color=Group))
  
  if(sample.times.IPRED) p <- p+stat_summary(data=df.ipred.samples,aes(x=Time,y=IPRED),geom="pointrange",fun.data="median_hilow")
  if(sample.times.DV) p <- p+stat_summary(data=df.ipred.samples,aes(x=Time,y=DV),geom="pointrange",fun.data="median_hilow")
  if(sample.times.DV.lines) p <- p + geom_line(aes(y=DV,group=ID),data=df.ipred.samples,alpha=alpha.sample.times.DV.lines)
  if(sample.times.DV.points) p <- p + geom_point(aes(y=DV,group=ID),data=df.ipred.samples,alpha=alpha.sample.times.DV.points)
  p <- p+ylab(y_lab)
  
  #,scales="free",space="free"
  #ggplot(data=df,aes(x=Time,y=PRED))+  geom_line()+facet_wrap(~Model,scales="free")
  #     p1 + stat_summary(data=df.ipred,aes(x=Time,y=IPRED,fill=Group),geom="ribbon",fun.data="median_hilow",alpha=0.3)
  #     p <- ggplot(data=df.ipred,aes(x=Time,y=IPRED))
  #     p + geom_line(aes(group=ID))
  #     p + geom_line(aes(group=ID,colour=Group))
  #     p + geom_line(alpha=0.3,aes(colour=Group,group=ID))
  #     p + geom_line(aes(group=ID))+facet_grid(~Group) 
  #     
  #     p + geom_line(aes(group=ID)) + stat_summary(geom="ribbon",fun.data="median_hilow",alpha=0.3,fill="red")
  #     
  #     p + stat_summary(geom="ribbon",fun.data="median_hilow",alpha=0.3) +
  #       geom_line(aes(x=Time,y=PRED,group=Group),data=df) + geom_point(aes(x=Time,y=PRED,group=Group),data=df.2,size=4)
  #     
  #     p + stat_summary(geom="ribbon",fun.data="median_hilow",alpha=0.3)+
  #       stat_summary(fun.y = median, #fun.ymin = min, fun.ymax = max,
  #                    colour = "red")
  #     
  #     stat_summary(aes(x=wt,y=obsAUC),geom="ribbon",fun.ymin="min",fun.ymax="max")
  #     p <- ggplot(data=df,aes(x=Time,y=IPRED))
  #     p  + stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max,
  #                       colour = "red")
  #     stat_summary(aes(group=Group,fill=Group),geom="ribbon",fun.data="median_hilow",alpha=0.3) +
  #     fun.ymin="min",fun.ymax="max"
  #     
  #     p+geom_line()+facet_grid(~Group)+geom_point(data=df.2,size=4)
  
  
  return( p ) 
}
