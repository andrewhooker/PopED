#' Create global variables in the PopED database
#' 
#' Function takes design variables from input files
#' and converts them to the global variables needed
#' in PopED.  Typically not used by the user.  Instead 
#' use the function \code{\link{create.poped.database}}.
#' 
#' @param globalStructure A PopED database 
#' @return A PopED database
#' @family poped_input


## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

convert_variables <- function(globalStructure){
                                        # takes design variables from design functions
                                        # and converts them to the global variables needed
                                        # in PopED

    globalStructure$Engine = list(Type=1,Version=version$version.string)

    design = globalStructure$design

    globalStructure$ga = zeros(globalStructure$m,globalStructure$na)
    globalStructure$gx = zeros(globalStructure$m,globalStructure$nx)
    globalStructure$gxt = zeros(globalStructure$m,globalStructure$maxni)
    globalStructure$gminxt = zeros(globalStructure$m,globalStructure$maxni)
    globalStructure$gmaxxt = zeros(globalStructure$m,globalStructure$maxni)
    globalStructure$gmaxa = zeros(globalStructure$m,globalStructure$na)
    globalStructure$gmina = zeros(globalStructure$m,globalStructure$na)

    globalStructure$gbpop=design$bpop
    globalStructure$gd=design$d

    if((test_mat_size(matrix(c(globalStructure$m, 1),nrow=1,byrow=T),design$ni,'ni')==1)){
        globalStructure$gni = design$ni
        if(is.null(dim(design$ni))) globalStructure$gni = rbind(design$ni)
    }

    
    if(size(design$model_switch,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$model_switch,globalStructure$m)
      design$model_switch <- matrix(tmp,ncol=length(design$model_switch),nrow=globalStructure$m,byrow=T)
    } 
    if(!is.matrix(design$model_switch)) design$model_switch  <- matrix(design$model_switch,nrow=1)
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$maxni),nrow=1,byrow=T),design$model_switch,'model_switch')==1)){
        if(is.null(globalStructure$global_model_switch)){
            globalStructure$global_model_switch=design$model_switch
        } else {
            globalStructure$global_model_switch[1:globalStructure$m,1:globalStructure$maxni]=design$model_switch
        }
    }

    if(size(design$xt,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$xt,globalStructure$m)
      design$xt <- matrix(tmp,ncol=length(design$xt),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$maxni),nrow=1,byrow=T),design$xt,'xt')==1)){
        globalStructure$gxt[1:globalStructure$m,1:globalStructure$maxni]=design$xt
    }
    globalStructure$gxt=test_for_zeros(globalStructure$gxt,globalStructure$ourzero)

    if(size(design$maxxt,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$maxxt,globalStructure$m)
      design$maxxt <- matrix(tmp,ncol=length(design$maxxt),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$maxni),nrow=1,byrow=T),design$maxxt,'maxxt')==1)){
        globalStructure$gmaxxt[1:globalStructure$m,1:globalStructure$maxni] = design$maxxt
    }
    globalStructure$gxt=test_for_max(globalStructure$gxt,globalStructure$gmaxxt)

    if(size(design$minxt,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$minxt,globalStructure$m)
      design$minxt <- matrix(tmp,ncol=length(design$minxt),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$maxni),nrow=1,byrow=T),design$minxt,'minxt')==1)){
        globalStructure$gminxt[1:globalStructure$m,1:globalStructure$maxni] = design$minxt
    }
    globalStructure$gminxt=test_for_zeros(globalStructure$gminxt,globalStructure$ourzero)
    globalStructure$gxt = test_for_min(globalStructure$gxt,globalStructure$gminxt)

    if(size(design$a,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$a,globalStructure$m)
      design$a <- matrix(tmp,ncol=length(design$a),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$na),nrow=1,byrow=T),design$a,'a')==1)){
        globalStructure$ga[1:globalStructure$m,1:globalStructure$na]=design$a
    }

    if(size(design$maxa,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$maxa,globalStructure$m)
      design$maxa <- matrix(tmp,ncol=length(design$maxa),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$na),nrow=1,byrow=T),design$maxa,'maxa')==1)){
        globalStructure$gmaxa[1:globalStructure$m,1:globalStructure$na]=design$maxa
    }
    globalStructure$ga = test_for_max(globalStructure$ga,globalStructure$gmaxa)

    if(size(design$mina,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$mina,globalStructure$m)
      design$mina <- matrix(tmp,ncol=length(design$mina),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$na),nrow=1,byrow=T),design$mina,'mina')==1)){
        globalStructure$gmina[1:globalStructure$m,1:globalStructure$na]=design$mina
    }
    globalStructure$ga= test_for_min(globalStructure$ga,globalStructure$gmina)

    if(size(design$x,1)==1 && globalStructure$m!=1){
      tmp <- rep(design$x,globalStructure$m)
      design$x <- matrix(tmp,ncol=length(design$x),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$nx),nrow=1,byrow=T),design$x,'x')==1)){
        globalStructure$gx[1:globalStructure$m,1:globalStructure$nx]=design$x
    }

    if(size(design$discrete_x,1)==1 && globalStructure$m!=1){
      tmp <- cell(globalStructure$m,size(design$discrete_x,2))
      for(i in i:globalStructure$m) tmp[i,]  <-  design$discrete_x
      design$discrete_x <- tmp
    } 
    
    if(length(design$groupsize)==1 && globalStructure$m!=1){
      tmp <- rep(design$groupsize,globalStructure$m)
      design$groupsize <- matrix(tmp,ncol=length(design$groupsize),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, 1),nrow=1,byrow=T),design$groupsize,'groupsize')==1)){
        globalStructure$groupsize=design$groupsize
    }

    if(size(design$G,1)==1 && globalStructure$m!=1){
      if(globalStructure$bUseGrouped_xt){ 
        tmp <- rep(design$G,globalStructure$m)      
      } else {
        tmp <- c()
        for(i in 1:globalStructure$m){
          tmp <- c(tmp,design$G+(i-1)*length(design$G))
        }
      }
      design$G <- matrix(tmp,ncol=length(design$G),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$maxni),nrow=1,byrow=T),design$G,'G')==1)){
        globalStructure$G=design$G
    }
    
    for(i in 1:size(design$G,1)){
        for(j in 1:size(design$G,2)){
            iGroupVal = design$G[i,j]
            tmp = matrix(1,size(globalStructure$gxt,1),size(globalStructure$gxt,2))*iGroupVal
            inters = (design$G==tmp)
            xt_val = globalStructure$gxt[i,j]
            xt_max = globalStructure$gmaxxt[i,j]
            xt_min = globalStructure$gminxt[i,j]
            for(k in 1:size(globalStructure$gxt,1)){
                for(l in 1:size(globalStructure$gxt,2)){
                    if((inters[k,l]==1)){
                        if((xt_val!=globalStructure$gxt[k,l])){
                            stop(sprintf('xt[%g,%g] are grouped with xt[%g,%g] but they do not have the same initial values.\n',i,j,k,l)) 
                        }
                        if((xt_max!=globalStructure$gmaxxt[k,l])){
                            stop(sprintf('xt[%g,%g] are grouped with xt[%g,%g] but they do not have the same maximal values.\n',i,j,k,l)) 
                        }
                        if((xt_min!=globalStructure$gminxt[k,l])){
                            stop(sprintf('xt[%g,%g] are grouped with xt[%g,%g] but they do not have the same minimum values.\n',i,j,k,l)) 
                        }
                    }
                }
            }
        }
    }
    
    for(i in 1:max(max(max(design$G)),0)){
        tmp = matrix(1,size(globalStructure$gxt,1),size(globalStructure$gxt,2))*i
        inters = (design$G==tmp) 
        if((sum(sum(inters))==0)){
            if((globalStructure$Engine$Type==1)){
                warning('MATLAB:noGroupMember','Time points grouped in group nr %g has no members.',i)
            } else {
                warning(sprintf('MATLAB:noGroupMember - Time points grouped in group nr %g has no members.',i))
            }
        }
    }
    
    if(size(design$Ga,1)==1 && globalStructure$m!=1){
      if(globalStructure$bUseGrouped_a){ 
        tmp <- rep(design$Ga,globalStructure$m)      
      } else {
        tmp <- c()
        for(i in 1:globalStructure$m){
          tmp <- c(tmp,design$Ga+(i-1)*length(design$Ga))
        }
      }
      design$Ga <- matrix(tmp,ncol=length(design$Ga),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$na),nrow=1,byrow=T),design$Ga,'Ga')==1)){
        if(is.null(dim(design$Ga))){
            design$Ga=matrix(data=design$Ga,nrow=1,byrow=TRUE)
        }
        globalStructure$Ga=design$Ga
    }
    
    for(i in 1:size(design$Ga,1)){
        if(any(size(design$Ga)==0)) next
        for(j in 1:size(design$Ga,2)){
            iGroupVal = design$Ga[i,j]
            tmp = matrix(1,size(globalStructure$ga,1),size(globalStructure$ga,2))*iGroupVal
            inters = (design$Ga==tmp)
            a_val = globalStructure$ga[i,j]
            a_max = globalStructure$gmaxa[i,j]
            a_min = globalStructure$gmina[i,j]
            for(k in 1:size(globalStructure$ga,1)){
                for(l in 1:size(globalStructure$ga,2)){
                    if((inters[k,l]==1)){
                        if((a_val!=globalStructure$ga[k,l])){
                            stop(sprintf('a[%g,%g] are grouped with a[%g,%g] but they do not have the same initial values.\n',i,j,k,l)) 
                        }
                        if((a_max!=globalStructure$gmaxa[k,l])){
                            stop(sprintf('a[%g,%g] are grouped with a[%g,%g] but they do not have the same maximal values.\n',i,j,k,l)) 
                        }
                        if((a_min!=globalStructure$gmina[k,l])){
                            stop(sprintf('a[%g,%g] are grouped with a[%g,%g] but they do not have the same minimum values.\n',i,j,k,l)) 
                        }
                    }
                }
            }
        }
    }
    
    if((!isempty(design$Ga))){
        for(i in 1:max(max(max(design$Ga)),0)){
            tmp = matrix(1,size(globalStructure$ga,1),size(globalStructure$ga,2))*i
            inters = (design$Ga==tmp) 
            if((sum(sum(inters))==0)){
                if((globalStructure$Engine$Type == 1)){
                    warning('MATLAB:noGroupMember','Covariates grouped in group nr %g has no members.',i)
                } else {
                    warning(sprintf('MATLAB:noGroupMember - Covariates grouped in group nr %g has no members.',i))
                }
            }
        }
    }
    
    if(size(design$Gx,1)==1 && globalStructure$m!=1){
      if(globalStructure$bUseGrouped_x){ 
        tmp <- rep(design$Gx,globalStructure$m)      
      } else {
        tmp <- c()
        for(i in 1:globalStructure$m){
          tmp <- c(tmp,design$Gx+(i-1)*length(design$Gx))
        }
      }
      design$Gx <- matrix(tmp,ncol=length(design$Gx),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, globalStructure$nx),nrow=1,byrow=T),design$Gx,'Gx')==1)){
        globalStructure$Gx=design$Gx
    }
    
    for(i in 1:size(design$Gx,1)){
        if(any(size(design$Gx)==0)) next
        for(j in 1:size(design$Gx,2)){
            iGroupVal = design$Gx[i,j]
            tmp = matrix(1,size(globalStructure$gx,1),size(globalStructure$gx,2))*iGroupVal
            inters = (design$Gx==tmp)
            x_val = globalStructure$gx[i,j]
            discrete_val = globalStructure$discrete_x[[i,j]]
            for(k in 1:size(globalStructure$gx,1)){
                for(l in 1:size(globalStructure$gx,2)){
                    if((inters[k,l]==1)){
                        discrete_val_2 = globalStructure$discrete_x[[k,l]]
                        if((x_val!=globalStructure$gx[k,l])){
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
            tmp = matrix(1,size(globalStructure$gx,1),size(globalStructure$gx,2))*i
            inters = (design$Gx==tmp) 
            if((sum(sum(inters))==0)){
                if((globalStructure$Engine$Type==1)){
                    warning('MATLAB:noGroupMember','Other design var grouped in group nr %g has no members.',i)
                } else {
                    warning(sprintf('MATLAB:noGroupMember - Other design var grouped in group nr %g has no members.',i))
                }
            }
        }
    }
    
    if(length(design$maxgroupsize)==1 && globalStructure$m!=1){
      tmp <- rep(design$maxgroupsize,globalStructure$m)
      design$maxgroupsize <- matrix(tmp,ncol=length(design$maxgroupsize),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, 1),nrow=1,byrow=T),design$maxgroupsize,'maxgroupsize')==1)){
        globalStructure$maxgroupsize=design$maxgroupsize
    }
    globalStructure$groupsize=test_for_max(globalStructure$groupsize,globalStructure$maxgroupsize)

    if(length(design$mingroupsize)==1 && globalStructure$m!=1){
      tmp <- rep(design$mingroupsize,globalStructure$m)
      design$mingroupsize <- matrix(tmp,ncol=length(design$mingroupsize),nrow=globalStructure$m,byrow=T)
    } 
    if((test_mat_size(matrix(c(globalStructure$m, 1),nrow=1,byrow=T),design$mingroupsize,'mingroupsize')==1)){
        globalStructure$mingroupsize=design$mingroupsize
    }
    globalStructure$groupsize=test_for_min(globalStructure$groupsize,globalStructure$mingroupsize)
    
    ## create ds_index if not already done
    if(is.null(globalStructure$ds_index)){ 
      unfixed_params <- get_unfixed_params(globalStructure)
      globalStructure$ds_index <- t(rep(0,length(unfixed_params$all)))
      globalStructure$ds_index[(length(unfixed_params$bpop)+1):length(globalStructure$ds_index)] <- 1
    }
    
    globalStructure$design <- design
    return( globalStructure) 
}
