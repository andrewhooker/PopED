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
    }    
    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$model_switch,'model_switch')==1)){
        if(is.null(poped.db$global_model_switch)){
            poped.db$global_model_switch=design$model_switch
        } else {
            poped.db$global_model_switch[1:poped.db$m,1:poped.db$maxni]=design$model_switch
        }
    }

    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$xt,'xt')==1)){
        poped.db$gxt[1:poped.db$m,1:poped.db$maxni]=design$xt
    }

    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$maxxt,'maxxt')==1)){
        poped.db$gmaxxt[1:poped.db$m,1:poped.db$maxni] = design$maxxt
    }

    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$minxt,'minxt')==1)){
        poped.db$gminxt[1:poped.db$m,1:poped.db$maxni] = design$minxt
    }
    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$a,'a')==1)){
        poped.db$ga[1:poped.db$m,1:poped.db$na]=design$a
    }

    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$maxa,'maxa')==1)){
        poped.db$gmaxa[1:poped.db$m,1:poped.db$na]=design$maxa
    }

    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$mina,'mina')==1)){
        poped.db$gmina[1:poped.db$m,1:poped.db$na]=design$mina
    }

    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$nx),nrow=1,byrow=T),design$x,'x')==1)){
        poped.db$gx[1:poped.db$m,1:poped.db$nx]=design$x
    }

    
    
    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$groupsize,'groupsize')==1)){
        poped.db$groupsize=design$groupsize
    }

    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$maxni),nrow=1,byrow=T),design$G,'G')==1)){
        poped.db$G=design$G
    }
    
    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$na),nrow=1,byrow=T),design$Ga,'Ga')==1)){
        if(is.null(dim(design$Ga))){
            design$Ga=matrix(data=design$Ga,nrow=1,byrow=TRUE)
        }
        poped.db$Ga=design$Ga
    }
    
    
    if((test_mat_size(matrix(c(poped.db$m, poped.db$nx),nrow=1,byrow=T),design$Gx,'Gx')==1)){
        poped.db$Gx=design$Gx
    }
    
    
    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$maxgroupsize,'maxgroupsize')==1)){
        poped.db$maxgroupsize=design$maxgroupsize
    }

    if((test_mat_size(matrix(c(poped.db$m, 1),nrow=1,byrow=T),design$mingroupsize,'mingroupsize')==1)){
        poped.db$mingroupsize=design$mingroupsize
    }
    
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
