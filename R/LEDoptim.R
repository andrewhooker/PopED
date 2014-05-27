#' Optimization function for D-family, E-family and Laplace approximated ED designs
#' 
#' Optimize the objective fucntion for D-family, E-family and Laplace approximated ED designs.  
#' Right now there is only one optimization algorithm used in this 
#' function 
#' \enumerate{
#' \item Adaptive random search. See \code{\link{RS_opt_gen}}.
#' }
#' This function takes information from the PopED database supplied as an argument.
#' The PopED database supplies information about the the model, parameters, design and methods to use.
#' Some of the arguments coming from the PopED database can be overwritten;  
#' if they are supplied then they are used instead of the arguments from the PopED database.
#' 
#' @inheritParams RS_opt
#' @inheritParams RS_opt_gen
#' @inheritParams create.poped.database
#' @inheritParams Doptim
#' @param fim_init The initial value of the FIM. If set to zero then it is computed.
#' @param ofv_init The inital OFV. If set to zero then it is computed.
#' @param trflag Should the optimization be output to the screen and to a file?
#' 
#' @family Optimize
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_ed.R
#' @example tests/testthat/examples_fcn_doc/examples_LEDoptim.R

LEDoptim <- function(poped.db,
                     model_switch=NULL,
                     ni=NULL,
                     xt=NULL,
                     x=NULL,
                     a=NULL,
                     bpopdescr=NULL,
                     ddescr=NULL,
                     maxxt=NULL,
                     minxt=NULL,
                     maxa=NULL,
                     mina=NULL,
                     ofv_init=0,
                     fim_init=0,
                     trflag=TRUE,
                     header_flag=TRUE,
                     footer_flag=TRUE,
                     opt_xt=poped.db$optsw[2],
                     opt_a=poped.db$optsw[4],
                     opt_x=poped.db$optsw[3],
                     out_file=NULL,
                     d_switch=FALSE,
                     use_laplace=T,
                     laplace.fim=FALSE, 
                     use_RS=poped.db$bUseRandomSearch,
                     ...){
  #+++++++++++++++++++++  D CONTINUOUS VARIABLE OPTIMIZATION FUNCTION
  
  # ------------- downsizing of general design
  downsize.list <- downsizing_general_design(poped.db)
  if(is.null(ni)) ni <- downsize.list$ni
  if(is.null(xt)) xt <- downsize.list$xt
  if(is.null(model_switch)) model_switch <- downsize.list$model_switch
  if(is.null(x)) x <- downsize.list$x
  if(is.null(a)) a <- downsize.list$a    
  if(is.null(bpopdescr)) bpopdescr <- downsize.list$bpop
  if(is.null(ddescr)) ddescr <- downsize.list$d
  if(is.null(maxxt)) maxxt <- downsize.list$maxxt
  if(is.null(minxt)) minxt <- downsize.list$minxt
  if(is.null(maxa)) maxa <- downsize.list$maxa
  if(is.null(mina)) mina <- downsize.list$mina
  
  if(sum(opt_x,opt_a,opt_xt)==0){
    cat("No optimization variables specified in input\n")
    return(invisible(list( xt= xt,x=x,a=a,fmf=fim_init,dmf=ofv_init,poped.db =poped.db ))) 
  }  
  
  # ----------------- initialization of optimization variables
  output <- calc_ofv_and_fim(poped.db,
                             ofv=ofv_init,
                             fim=fim_init,
                             d_switch=d_switch,  
                             bpopdescr=bpopdescr, 
                             ddescr=ddescr,
                             model_switch=model_switch, 
                             ni=ni, 
                             xt=xt, 
                             x=x, 
                             a=a, 
                             use_laplace=use_laplace,
                             laplace.fim=laplace.fim, 
                             ...)
  fmf <- output$fim
  dmf <- output$ofv
  fmf_init <- fmf
  dmf_init <- dmf
  xtopt=xt
  xopt=x
  aopt=a
  
  #--------------------- write out info to a file
  #if((trflag && header_flag)){
  fn=blockheader(poped.db,name="ED_Laplace_opt",e_flag=!d_switch,
                 opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,
                 fmf=fmf,dmf=dmf,
                 bpop=bpopdescr,d=ddescr,docc=poped.db$docc,sigma=poped.db$sigma,
                 out_file=out_file,
                 trflag=trflag,
                 header_flag=header_flag,
                 ...)
  #} 
  
  
  if((use_RS) ){#If we want to perform random search
    # ----------------- RANDOM SEARCH BEGINS HERE                           
    rs_output <- RS_opt_gen(poped.db,
                            d_switch=d_switch,
                            use_laplace=use_laplace,
                            fmf=fmf,
                            dmf=dmf,
                            ni=ni, 
                            xt=xtopt, 
                            model_switch=model_switch, 
                            x=xopt, 
                            a=aopt, 
                            opt_xt=opt_xt,
                            opt_a=opt_a,
                            opt_x=opt_x,
                            bpopdescr=bpopdescr, 
                            ddescr=ddescr, 
                            maxxt=maxxt, 
                            minxt=minxt,
                            maxa=maxa,
                            mina=mina,
                            trflag = trflag,
                            header_flag=FALSE,
                            footer_flag=FALSE,
                            out_file=fn,
                            laplace.fim=laplace.fim,                             
                            ...)
    
    
    xtopt=rs_output$xtopt
    xopt=rs_output$xopt
    aopt=rs_output$aopt
    fmf=rs_output$fmf
    dmf=rs_output$dmf
    poped.db=rs_output$poped.db
    
    bestxt=rs_output$xtopt
    bestx=rs_output$xopt
    besta=rs_output$aopt
    best_dmf=dmf
    
    
  } 
  if(((FALSE))){
    #if(((poped.db$bUseBFGSMinimizer))){
    # ----------------- type of optimization determination
    axt=poped.db$optsw[2]*poped.db$cfaxt*matrix(1,poped.db$m,poped.db$maxni)
    aa=poped.db$optsw[4]*poped.db$cfaa*matrix(1,poped.db$m,poped.db$na)
    optxt=opt_xt
    optx=opt_x
    opta=opt_a
    
    bfgs_init=matrix(0,0,0)
    
    
    x_k=matrix(0,0,0)
    lb=matrix(0,0,0)
    ub=matrix(0,0,0)
    if((optxt==TRUE)){
      index=t(1:numel(xtopt))
      if(poped.db$bUseGrouped_xt){
        returnArgs <- unique(poped.db$design$G) 
        temp <- returnArgs[[1]]
        index <- returnArgs[[2]]
        temp2 <- returnArgs[[3]]
      }
      index=index[minxt!=maxxt]
      x_k=t(t(xtopt[index]))
      lb=t(t(minxt[index]))
      ub=t(t(maxxt[index]))
    }
    if((opta==TRUE)){
      index=t(1:numel(aopt))
      if(poped.db$bUseGrouped_a){
        returnArgs <- unique(poped.db$design$Ga) 
        temp1 <- returnArgs[[1]]
        index <- returnArgs[[2]]
        temp2 <- returnArgs[[3]]
      }
      index=index[mina!=maxa]
      x_k=t(t(c(x_k,aopt[index])))
      lb=t(t(c(lb,mina[index])))
      ub=t(t(c(ub,maxa[index])))
      
      #x_k(end+index)=aopt[index]
      #lb(end+index)=mina[index]
      #ub(end+index)=maxa[index]
    }
    if((any(x_k<lb))){
      x_k[x_k<lb]=lb[x_k<lb]
    }
    if((isempty(bfgs_init) || any(x_k!=bfgs_init))){
      bfgs_init=x_k
      fprintf('Starting BGFS minimization with OFV of %g \n', best_dmf)
      returnArgs <- bfgsb_min('ed_laplace_ofv',
                              list(model_switch,aa,axt,poped.db$groupsize,ni,
                                   xtopt,xopt,aopt,bpop,d,poped.db$sigma,docc_full,poped.db,
                                   return_gradient=T,
                                   optxt=optxt, opta=optx, x=x),
                              x_k,lb,ub,options) 
      x_opt  <- returnArgs[[1]]
      f_k <- returnArgs[[2]]
      B_k <- returnArgs[[3]]
      if(optxt){
        notfixed=minxt!=maxxt
        if(poped.db$bUseGrouped_xt){
          xtopt[notfixed]=x_opt[poped.db$design$G[notfixed]]
          x_opt <- x_opt[-(1:numel(unique(poped.db$design$G[notfixed])))]
        } else {
          xtopt[notfixed]=x_opt[1:numel(xtopt[notfixed])]
          x_opt <- x_opt[-(1:numel(xtopt[notfixed]))]
        }
      }
      
      if(opta){
        notfixed=mina!=maxa
        
        if(poped.db$bUseGrouped_a){
          aopt[notfixed]=x_opt(poped.db$design$Ga[notfixed])
        } else {
          aopt[notfixed]=x_opt
        }
      }
      dmf=ed_laplace_ofv(model_switch,poped.db$groupsize,ni,xtopt,xopt,aopt,
                         bpopdescr,ddescr,poped.db$covd,poped.db$sigma,poped.db$docc,poped.db)
      
      fprintf('BFGS minimization finished. New OFV: %d \n', dmf)
      if((dmf>best_dmf)){
        best_dmf=dmf
        bestxt=xtopt
        bestx=xopt
        besta=aopt
      }
    }
  }
  xt=xtopt
  x=xopt
  a=aopt
  dmf=best_dmf
  
  
  #--------- Write results
  #if((trflag)){
  #  if(footer_flag){
  blockfinal(fn,fmf,dmf,poped.db$groupsize,ni,xtopt,xopt,aopt,model_switch,
             bpopdescr,ddescr,poped.db$docc,poped.db$sigma,poped.db,
             opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,
             fmf_init=fmf_init,dmf_init=dmf_init,out_file=out_file,
             trflag=trflag,
             footer_flag=footer_flag,
             ...)
  #  }
  #}
  
  return(invisible(list( xt= xt,x=x,a=a,fmf=fmf,dmf=dmf,poped.db =poped.db ))) 
}
