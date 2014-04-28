#' Plot the efficience of windows 
#' 
#' Function plots the efficiency of windows around the optimal design points.
#' 
#' @param poped.db A poped database
#' @param iNumSimulations The number of efficiency calculations to make.
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param ... Extra arguments passed to \code{evaluate.fim}
#' @param xt_windows The distance on one direction from the optimal sample times.  Can be a number or a matrix of the same size as 
#' the xt matrix found in \code{poped.db$gxt}. 
#' 
#' @return A \link[ggplot2]{ggplot2} object.
#' 
#' @family evaluate_design
#' @family Simulation
#' @family Graphics
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_plot_efficiency_of_windows.R


plot_efficiency_of_windows <- function(poped.db,xt_windows,
                                       iNumSimulations=100,
                                       #a_windows=NULL,x_windows=NULL,
                                       #d_switch=poped.db$d_switch,
                                       ...){
  
  
  design = poped.db$design
  design$xt <- poped.db$gxt
  design$x <- poped.db$gx
  design$a <- poped.db$ga
  
  p = get_fim_size(poped.db) #%Size of FIM
  
  xt_val = design$xt
  a_val  = design$a
  x_val  = design$x
  
  ref_fmf <- evaluate.fim(poped.db,xt=xt_val,a=a_val,x=x_val,...)
  
  #NumRows = size(xt_windows,1)*size(xt_windows,2)/2+size(a_windows,1)*size(a_windows,2)/2+size(x_windows,1)*size(x_windows,2)/2+p+1
  
  
  eff = zeros(1,iNumSimulations)
  
  fulld = getfulld(design$d[,2],design$covd)
  fulldocc = getfulld(design$docc[,2,drop=F],design$covdocc)
  
  if(!is.null(xt_windows)){
    xt_val_min <- xt_val - xt_windows
    xt_val_max <- xt_val + xt_windows
    xt_val_min[xt_val_min<design$minxt] = design$minxt[xt_val_min<design$minxt]
    xt_val_max[xt_val_max>design$maxxt] = design$maxxt[xt_val_max>design$maxxt]
  }  
  
  for(i in 1:iNumSimulations){
    if(!is.null(xt_windows)){
      xt_new <- rand(size(xt_val,1),size(xt_val,2))*abs(xt_val_max - xt_val_min)+xt_val_min
      if((poped.db$bUseGrouped_xt)){
        for(j in 1:size(xt_val,1) ){
          for(k in min(poped.db$G[j,]):max(poped.db$G[j,])){
            tmp.lst <- xt_new[j,poped.db$G[j,]==k]
            if(length(tmp.lst)!=1){
              tmp <- sample(tmp.lst,1)
              xt_new[j,poped.db$G[j,]==k]=tmp
            }
          }
        }
      }      
      xt_val=xt_new
    }
    #     for(j in 1:size(xt_val,1) ){#Get simulated xt
    #       for(k in 1:size(xt_val,2)){
    #         xt_val[j,k] = xt_val[j,k] + rand*(xt_val_min[j,k]-xt_val_max(j,(k-1)*2+1))
    #       }
    #     }
    #     for(j in 1:size(a_windows,1) ){#Get simulated a
    #       for(k in 1:size(a_windows,2)/2){
    #         a_val(j,k) = a_windows(j,(k-1)*2+1) + rand*(a_windows(j,(k)*2)-a_windows(j,(k-1)*2+1))
    #       }
    #     }
    #     for(j in 1:size(x_windows,1) ){#Get simulated x
    #       for(k in 1:size(x_windows,2)/2){
    #         x_val(j,k) = x_windows(j,(k-1)*2+1) + rand*(x_windows(j,(k)*2)-x_windows(j,(k-1)*2+1))
    #       }
    #     }
    
    #     if((poped.db$bUseGrouped_xt)){
    #       xt_val=group_matrix(xt_val,poped.db$G)
    #     }
    #     if((poped.db$bUseGrouped_a)){
    #       a_val=group_matrix(a_val,poped.db$Ga)
    #     }
    #     if((poped.db$bUseGrouped_x)){
    #       x_val=group_matrix(x_val,poped.db$Gx)
    #     }
    #     
    #     if((!bParallel)){
    #if((poped.db$d_switch)){
    fmf <- evaluate.fim(poped.db,xt=xt_val,...)
    #       returnArgs <-  mftot(model_switch,poped.db$groupsize,ni,xt_val,x_val,a_val,bpop[,2],fulld,design$sigma,fulldocc,poped.db) 
    #       fmf <- returnArgs[[1]]
    #       poped.db <- returnArgs[[2]]
    eff[1,i] = ofv_criterion(ofv_fim(fmf,poped.db),p,poped.db)/ofv_criterion(ofv_fim(ref_fmf,poped.db),p,poped.db)
    #       } else {
    #         eff(i) = (ofv_fim(fmf,poped.db)/optdetfim)
    #       }
    #     } else {
    #       returnArgs <- ed_mftot(model_switch,poped.db$groupsize,ni,xt_val,x_val,a_val,bpop,d,poped.db$covd,design$sigma,poped.db$docc,poped.db) 
    #       fmf <- returnArgs[[1]]
    #       ED_det <- returnArgs[[2]]
    #       poped.db <- returnArgs[[3]]
    #       if((bNormalizedEfficiency)){
    #         eff(i) = ofv_efficiency(ofv_criterion(ED_det,p,poped.db),ofv_criterion(optdetfim,p,poped.db))
    #       } else {
    #         eff(i) = (ED_det/optdetfim)
    #       }
    #     }
    #       
    #       if((bStoreInFile)){
    #         tmp_xt= reshape_matlab(xt_val',size(xt_windows,1)*size(xt_windows,2)/2,1)
    #             tmp_a = reshape_matlab(a_val',size(a_windows,1)*size(a_windows,2)/2,1)
    #         tmp_x = reshape_matlab(x_val',size(x_windows,1)*size(x_windows,2)/2,1)
    #             fileData(i,1:length(tmp_xt)) = tmp_xt
    #             fileData(i,1+length(tmp_xt):length(tmp_xt)+length(tmp_a)) = tmp_a
    #             fileData(i,1+length(tmp_a)+length(tmp_xt):length(tmp_xt)+length(tmp_a)+length(tmp_x)) = tmp_x
    #             try
    #                 imf = inv(fmf)
    #                  returnArgs <- get_cv(diag_matlab(imf),bpop,d,poped.db$docc,design$sigma,poped.db) 
    # params <- returnArgs[[1]]
    #  params_cv <- returnArgs[[2]]
    #             catch
    #                 params_cv = zeros(0,p)
    #             }
    #             fileData(i,1+length(tmp_a)+length(tmp_xt)+length(tmp_x):length(tmp_a)+length(tmp_xt)+length(tmp_x)+p)=params_cv
    #             fileData(i,NumRows) = eff(i)
    #         }
    #     } else {
    #          designsin = update_designinlist(designsin,poped.db$groupsize,ni,xt_val,x_val,a_val,i,0)
    #     }
    # }
    # 
    # if((bParallel) ){#Execute everything in parallel instead
    #     designsout = execute_parallel(designsin,poped.db)
    #     for(i in 1:iNumSimulations){
    #         ED_det=designsout[[i]].ofv
    #         fmf=designsout[[i]]$FIM
    #         
    #         if((bNormalizedEfficiency)){
    #                 eff(i) = ofv_efficiency(ofv_criterion(ED_det,p,poped.db),ofv_criterion(optdetfim,p,poped.db))
    #             } else {
    #                 eff(i) = (ED_det/optdetfim)
    #          }
    #         
    #         if((bStoreInFile)){
    #             tmp_xt= reshape_matlab(designsin[[i]].xt',size(xt_windows,1)*size(xt_windows,2)/2,1)
    #         tmp_a = reshape_matlab(designsin[[i]].a',size(a_windows,1)*size(a_windows,2)/2,1)
    #             tmp_x = reshape_matlab(designsin[[i]].x',size(x_windows,1)*size(x_windows,2)/2,1)
    #         fileData(i,1:length(tmp_xt)) = tmp_xt
    #         fileData(i,1+length(tmp_xt):length(tmp_xt)+length(tmp_a)) = tmp_a
    #         fileData(i,1+length(tmp_a)+length(tmp_xt):length(tmp_xt)+length(tmp_a)+length(tmp_x)) = tmp_x
    #         try
    #         imf = inv(fmf)
    #         returnArgs <- get_cv(diag_matlab(imf),bpop,d,poped.db$docc,design$sigma,poped.db) 
    #         params <- returnArgs[[1]]
    #         params_cv <- returnArgs[[2]]
    #         catch
    #         params_cv = zeros(0,p)
    #       }
    #       fileData(i,1+length(tmp_a)+length(tmp_xt)+length(tmp_x):length(tmp_a)+length(tmp_xt)+length(tmp_x)+p)=params_cv
    #       fileData(i,NumRows) = eff(i)
    #     }
    #   }
    # }
    # 
    # 
    # min_eff = min(eff)
    # max_eff = max(eff)
    # mean_eff = mean(eff)
    # 
    # if((bStoreInFile)){
    #   csvwrite(strFileName,fileData)
    # }
    # 
    # if((bPlot)){
    #   figure
    #   ##hold on
    #   if((bNormalizedEfficiency)){
    #     title('Normalized Efficiency of samples from sampling/covariate - windows')
    #   } else {
    #     title('Efficiency of samples from sampling/covariate - windows')
    #   }
    #   ylim(matrix(c(100*min(eff)-10,max(100,100*max(eff))),nrow=1,byrow=T))
    #   plot(1:iNumSimulations,eff*100,'-b')    
    #   xlabel('Sample number')
    #   if((bNormalizedEfficiency)){
    #     ylabel('Normalized Efficiency (%)')
    #   } else {
    #     ylabel('Efficiency (%)')
    #   }
    #   ##hold off
    # }
    # 
    # return( eff min_eff mean_eff max_eff ) 
  }
  efficiency <- eff[1,]*100
  df <- data.frame(efficiency,sample=c(1:iNumSimulations))
  p <- ggplot(data=df,aes(x=sample,y=efficiency)) #+ labs(colour = NULL)
  #p <- ggplot(data=df,aes(x=Efficiency)) #+ labs(colour = NULL)
  #p+geom_histogram()
  p <- p+geom_line()+geom_point()+ylab("Normalized Efficiency (%)")
  return(p)
}