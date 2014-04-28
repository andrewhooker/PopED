#' Create global variables in the PopED database
#' 
#' Function takes design variables from input files
#' and converts them to the global variables needed
#' in PopED.  Typically not used by the user.  Instead 
#' use the function \code{\link{create.poped.database}}.
#' 
#' @param poped.db A PopED database 
#' @return A PopED database
#' @family poped_input
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R


## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

convert_variables <- function(poped.db){
                                        # takes design variables from design functions
                                        # and converts them to the global variables needed
                                        # in PopED

    poped.db$Engine = list(Type=1,Version=version$version.string)

    design = poped.db$design

    poped.db$ga = zeros(poped.db$m,poped.db$na)
    poped.db$gx = zeros(poped.db$m,poped.db$nx)
    poped.db$gxt = zeros(poped.db$m,poped.db$maxni)
    poped.db$gminxt = zeros(poped.db$m,poped.db$maxni)
    poped.db$gmaxxt = zeros(poped.db$m,poped.db$maxni)
    poped.db$gmaxa = zeros(poped.db$m,poped.db$na)
    poped.db$gmina = zeros(poped.db$m,poped.db$na)

    poped.db$gbpop=design$bpop
    poped.db$gd=design$d

    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$ni,'ni')==1)){
        poped.db$gni = design$ni
        if(is.null(dim(design$ni))) poped.db$gni = rbind(design$ni)
    }

    
    if(size(design$model_switch,1)==1 && poped.db$m!=1){
      tmp <- rep(design$model_switch,poped.db$m)
      design$model_switch <- matrix(tmp,ncol=length(design$model_switch),nrow=poped.db$m,byrow=T)
    } 
    if(!is.matrix(design$model_switch)) design$model_switch  <- matrix(design$model_switch,nrow=1)
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$model_switch,'model_switch')==1)){
        if(is.null(poped.db$global_model_switch)){
            poped.db$global_model_switch=design$model_switch
        } else {
            poped.db$global_model_switch[1:poped.db$m,1:poped.db$maxni]=design$model_switch
        }
    }

    if(size(design$xt,1)==1 && poped.db$m!=1){
      tmp <- rep(design$xt,poped.db$m)
      design$xt <- matrix(tmp,ncol=length(design$xt),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$xt,'xt')==1)){
        poped.db$gxt[1:poped.db$m,1:poped.db$maxni]=design$xt
    }
    poped.db$gxt=test_for_zeros(poped.db$gxt,poped.db$ourzero)

    if(size(design$maxxt,1)==1 && poped.db$m!=1){
      tmp <- rep(design$maxxt,poped.db$m)
      design$maxxt <- matrix(tmp,ncol=length(design$maxxt),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$maxxt,'maxxt')==1)){
        poped.db$gmaxxt[1:poped.db$m,1:poped.db$maxni] = design$maxxt
    }
    poped.db$gxt=test_for_max(poped.db$gxt,poped.db$gmaxxt)

    if(size(design$minxt,1)==1 && poped.db$m!=1){
      tmp <- rep(design$minxt,poped.db$m)
      design$minxt <- matrix(tmp,ncol=length(design$minxt),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$minxt,'minxt')==1)){
        poped.db$gminxt[1:poped.db$m,1:poped.db$maxni] = design$minxt
    }
    poped.db$gminxt=test_for_zeros(poped.db$gminxt,poped.db$ourzero)
    poped.db$gxt = test_for_min(poped.db$gxt,poped.db$gminxt)

    if(size(design$a,1)==1 && poped.db$m!=1){
      tmp <- rep(design$a,poped.db$m)
      design$a <- matrix(tmp,ncol=length(design$a),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$a,'a')==1)){
        poped.db$ga[1:poped.db$m,1:poped.db$na]=design$a
    }

    if(size(design$maxa,1)==1 && poped.db$m!=1){
      tmp <- rep(design$maxa,poped.db$m)
      design$maxa <- matrix(tmp,ncol=length(design$maxa),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$maxa,'maxa')==1)){
        poped.db$gmaxa[1:poped.db$m,1:poped.db$na]=design$maxa
    }
    poped.db$ga = test_for_max(poped.db$ga,poped.db$gmaxa)

    if(size(design$mina,1)==1 && poped.db$m!=1){
      tmp <- rep(design$mina,poped.db$m)
      design$mina <- matrix(tmp,ncol=length(design$mina),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$mina,'mina')==1)){
        poped.db$gmina[1:poped.db$m,1:poped.db$na]=design$mina
    }
    poped.db$ga= test_for_min(poped.db$ga,poped.db$gmina)

    if(size(design$x,1)==1 && poped.db$m!=1){
      tmp <- rep(design$x,poped.db$m)
      design$x <- matrix(tmp,ncol=length(design$x),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$nx),nrow=1,byrow=T),design$x,'x')==1)){
        poped.db$gx[1:poped.db$m,1:poped.db$nx]=design$x
    }

    if(size(design$discrete_x,1)==1 && poped.db$m!=1){
      tmp <- cell(poped.db$m,size(design$discrete_x,2))
      for(i in i:poped.db$m) tmp[i,]  <-  design$discrete_x
      design$discrete_x <- tmp
    } 
    
    if(length(design$groupsize)==1 && poped.db$m!=1){
      tmp <- rep(design$groupsize,poped.db$m)
      design$groupsize <- matrix(tmp,ncol=length(design$groupsize),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$groupsize,'groupsize')==1)){
        poped.db$groupsize=design$groupsize
    }

    if(size(design$G,1)==1 && poped.db$m!=1){
      if(poped.db$bUseGrouped_xt){ 
        tmp <- rep(design$G,poped.db$m)      
      } else {
        tmp <- c()
        for(i in 1:poped.db$m){
          tmp <- c(tmp,design$G+(i-1)*length(design$G))
        }
      }
      design$G <- matrix(tmp,ncol=length(design$G),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$G,'G')==1)){
        poped.db$G=design$G
    }
    
    for(i in 1:size(design$G,1)){
        for(j in 1:size(design$G,2)){
            iGroupVal = design$G[i,j]
            tmp = matrix(1,size(poped.db$gxt,1),size(poped.db$gxt,2))*iGroupVal
            inters = (design$G==tmp)
            xt_val = poped.db$gxt[i,j]
            xt_max = poped.db$gmaxxt[i,j]
            xt_min = poped.db$gminxt[i,j]
            for(k in 1:size(poped.db$gxt,1)){
                for(l in 1:size(poped.db$gxt,2)){
                    if((inters[k,l]==1)){
                        if((xt_val!=poped.db$gxt[k,l])){
                            stop(sprintf('xt[%g,%g] are grouped with xt[%g,%g] but they do not have the same initial values.\n',i,j,k,l)) 
                        }
                        if((xt_max!=poped.db$gmaxxt[k,l])){
                            stop(sprintf('xt[%g,%g] are grouped with xt[%g,%g] but they do not have the same maximal values.\n',i,j,k,l)) 
                        }
                        if((xt_min!=poped.db$gminxt[k,l])){
                            stop(sprintf('xt[%g,%g] are grouped with xt[%g,%g] but they do not have the same minimum values.\n',i,j,k,l)) 
                        }
                    }
                }
            }
        }
    }
    
    for(i in 1:max(max(max(design$G)),0)){
        tmp = matrix(1,size(poped.db$gxt,1),size(poped.db$gxt,2))*i
        inters = (design$G==tmp) 
        if((sum(sum(inters))==0)){
            if((poped.db$Engine$Type==1)){
                warning('MATLAB:noGroupMember','Time points grouped in group nr %g has no members.',i)
            } else {
                warning(sprintf('MATLAB:noGroupMember - Time points grouped in group nr %g has no members.',i))
            }
        }
    }
    
    if(size(design$Ga,1)==1 && poped.db$m!=1){
      if(poped.db$bUseGrouped_a){ 
        tmp <- rep(design$Ga,poped.db$m)      
      } else {
        tmp <- c()
        for(i in 1:poped.db$m){
          tmp <- c(tmp,design$Ga+(i-1)*length(design$Ga))
        }
      }
      design$Ga <- matrix(tmp,ncol=length(design$Ga),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$Ga,'Ga')==1)){
        if(is.null(dim(design$Ga))){
            design$Ga=matrix(data=design$Ga,nrow=1,byrow=TRUE)
        }
        poped.db$Ga=design$Ga
    }
    
    for(i in 1:size(design$Ga,1)){
        if(any(size(design$Ga)==0)) next
        for(j in 1:size(design$Ga,2)){
            iGroupVal = design$Ga[i,j]
            tmp = matrix(1,size(poped.db$ga,1),size(poped.db$ga,2))*iGroupVal
            inters = (design$Ga==tmp)
            a_val = poped.db$ga[i,j]
            a_max = poped.db$gmaxa[i,j]
            a_min = poped.db$gmina[i,j]
            for(k in 1:size(poped.db$ga,1)){
                for(l in 1:size(poped.db$ga,2)){
                    if((inters[k,l]==1)){
                        if((a_val!=poped.db$ga[k,l])){
                            stop(sprintf('a[%g,%g] are grouped with a[%g,%g] but they do not have the same initial values.\n',i,j,k,l)) 
                        }
                        if((a_max!=poped.db$gmaxa[k,l])){
                            stop(sprintf('a[%g,%g] are grouped with a[%g,%g] but they do not have the same maximal values.\n',i,j,k,l)) 
                        }
                        if((a_min!=poped.db$gmina[k,l])){
                            stop(sprintf('a[%g,%g] are grouped with a[%g,%g] but they do not have the same minimum values.\n',i,j,k,l)) 
                        }
                    }
                }
            }
        }
    }
    
    if((!isempty(design$Ga))){
        for(i in 1:max(max(max(design$Ga)),0)){
            tmp = matrix(1,size(poped.db$ga,1),size(poped.db$ga,2))*i
            inters = (design$Ga==tmp) 
            if((sum(sum(inters))==0)){
                if((poped.db$Engine$Type == 1)){
                    warning('MATLAB:noGroupMember','Covariates grouped in group nr %g has no members.',i)
                } else {
                    warning(sprintf('MATLAB:noGroupMember - Covariates grouped in group nr %g has no members.',i))
                }
            }
        }
    }
    
    if(size(design$Gx,1)==1 && poped.db$m!=1){
      if(poped.db$bUseGrouped_x){ 
        tmp <- rep(design$Gx,poped.db$m)      
      } else {
        tmp <- c()
        for(i in 1:poped.db$m){
          tmp <- c(tmp,design$Gx+(i-1)*length(design$Gx))
        }
      }
      design$Gx <- matrix(tmp,ncol=length(design$Gx),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, poped.db$nx),nrow=1,byrow=T),design$Gx,'Gx')==1)){
        poped.db$Gx=design$Gx
    }
    
    for(i in 1:size(design$Gx,1)){
        if(any(size(design$Gx)==0)) next
        for(j in 1:size(design$Gx,2)){
            iGroupVal = design$Gx[i,j]
            tmp = matrix(1,size(poped.db$gx,1),size(poped.db$gx,2))*iGroupVal
            inters = (design$Gx==tmp)
            x_val = poped.db$gx[i,j]
            discrete_val = poped.db$discrete_x[[i,j]]
            for(k in 1:size(poped.db$gx,1)){
                for(l in 1:size(poped.db$gx,2)){
                    if((inters[k,l]==1)){
                        discrete_val_2 = poped.db$discrete_x[[k,l]]
                        if((x_val!=poped.db$gx[k,l])){
                            stop(sprintf('x[%g,%g] are grouped with x[%g,%g] but they do not have the same initial values.\n',i,j,k,l)) 
                        }
                        if((size(discrete_val_2,1)!=size(discrete_val,1) || size(discrete_val_2,2)!=size(discrete_val,2) || any(discrete_val_2!=discrete_val))){
                            stop(sprintf('x[%g,%g] are grouped with x[%g,%g] but they do not have the same "possible" discrete values.\n',i,j,k,l)) 
                        }
                    }
                }
            }
        }
    }
    
    if((!isempty(design$Gx))){
        for(i in 1:max(max(max(design$Gx)),0)){
            tmp = matrix(1,size(poped.db$gx,1),size(poped.db$gx,2))*i
            inters = (design$Gx==tmp) 
            if((sum(sum(inters))==0)){
                if((poped.db$Engine$Type==1)){
                    warning('MATLAB:noGroupMember','Other design var grouped in group nr %g has no members.',i)
                } else {
                    warning(sprintf('MATLAB:noGroupMember - Other design var grouped in group nr %g has no members.',i))
                }
            }
        }
    }
    
    if(length(design$maxgroupsize)==1 && poped.db$m!=1){
      tmp <- rep(design$maxgroupsize,poped.db$m)
      design$maxgroupsize <- matrix(tmp,ncol=length(design$maxgroupsize),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$maxgroupsize,'maxgroupsize')==1)){
        poped.db$maxgroupsize=design$maxgroupsize
    }
    poped.db$groupsize=test_for_max(poped.db$groupsize,poped.db$maxgroupsize)

    if(length(design$mingroupsize)==1 && poped.db$m!=1){
      tmp <- rep(design$mingroupsize,poped.db$m)
      design$mingroupsize <- matrix(tmp,ncol=length(design$mingroupsize),nrow=poped.db$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$mingroupsize,'mingroupsize')==1)){
        poped.db$mingroupsize=design$mingroupsize
    }
    poped.db$groupsize=test_for_min(poped.db$groupsize,poped.db$mingroupsize)
    
    ## create ds_index if not already done
    if(is.null(poped.db$ds_index)){ 
      unfixed_params <- get_unfixed_params(poped.db)
      poped.db$ds_index <- t(rep(0,length(unfixed_params$all)))
      poped.db$ds_index[(length(unfixed_params$bpop)+1):length(poped.db$ds_index)] <- 1
    } else {
      if(!is.matrix(poped.db$ds_index)) poped.db$ds_index <- matrix(poped.db$ds_index,1,length(poped.db$ds_index))
    }
    
    poped.db$design <- design
    return( poped.db) 
}
