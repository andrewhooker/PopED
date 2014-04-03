#' Trace optimization routines
#' 
#' A helper function for writing output to the screen and files when optimizing.
#' 
#' @param fn A file to output information to. Can also be the screen if \code{''}.
#' @param it the interation number.
#' @param xtopt The matrix defining optimal sampling schedule.
#' @param xopt The cell structure defining the optimal discrete design variables.
#' @param aopt The matrix defining the optimal continuous design variables.
#' @param gxt The matrix defining the current sampling schedule.
#' @param ga The matrix defining the current continuous design variables.
#' @param dmf The current OFV.
#' @param diff The difference from the previous iteration.
#' @param ixt If xt Gradient Inversion Occured or not.
#' @param ia If a Gradient Inversion Occured or not.
#' @param itvector The iteration vector
#' @param dmfvector The dmf vector.
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim 
#' @inheritParams create.poped.database
#' @param opt_samps Are the nuber of sample times per group being optimized?
#' @param opt_inds Are the nuber of individuals per group being optimized?

## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

Dtrace <- function(fn,it,ni,xtopt,xopt,aopt,gxt,ga,dmf,diff,ixt,
                   ia, 
                   itvector,dmfvector,poped.db,
                   opt_xt=poped.db$optsw[2],
                   opt_a=poped.db$optsw[4],opt_x=poped.db$optsw[3],
                   opt_samps=poped.db$optsw[1],opt_inds=poped.db$optsw[5],
                   rsit=poped.db$rsit){


if((it==0)){
    fprintf(fn,'*******************************\nInitial Value\n ')
    fprintf(fn,'OFV(mf) = %g\n',dmf)
    fprintf(fn,'*******************************\n\n')
}
if((it!=0)){
#     if((poped.db$bShowGraphs==TRUE)){
#         if((poped.db$Engine$Type==1)){
#             #set(0,'CurrentFigure',1)
#         } else {
#             #figure(1)
#         }
#         clf
#         optSum = poped.db$optsw[2]+poped.db$optsw[3]+poped.db$optsw[4]
#         numRows = 1
#         if((optSum>1)){
#             numRows = 2
#         }
# 
#         if((poped.db$optsw[2])){
#             subplot(numRows,2,1)
#             title('The current sampling times for each group')
#             xlabel('Sampling time')
#             ylabel('Group nr')
#             ##hold on
#             for(i in 1:poped.db$m){
#                 plot(xtopt[i,1:poped.db$gni[i]],matrix(1,1,poped.db$gni[i])*i,'b*','linestyle','none')
#             }
#             ##hold off
#         }
#         if((poped.db$optsw[3]==TRUE)){
#             subplot(numRows,2,1+poped.db$optsw[2])
#             title('The current discrete var. for each group')
#             xlabel('Discrete var.-value')
#             ylabel('Group nr')
#             ##hold on
#             for(i in 1:poped.db$m){
#                 plot(xopt(i,),matrix(1,1,size(xopt,2))*i,'b*')
#             }
#             ##hold off
#         }
#         if((poped.db$optsw[4]==TRUE)){
#             subplot(numRows,2,1+poped.db$optsw[2]+poped.db$optsw[3])
#             title('The current covariates for each group')
#             xlabel('Covariate-value')
#             ylabel('Group nr')
#             ##hold on
#             for(i in 1:poped.db$m){
#                 plot(aopt[i,],matrix(1,1,size(aopt,2))*i,'b*')
#             }
#             ##hold off
#         }
# 
#         subplot(numRows,2,1+optSum)
#         ##hold on
#         title('OFV(FIM)')
#         xlabel('Iteration')
#         ylabel('Value')
#         plot(itvector,dmfvector,'-b')
#         ##hold off
#         drawnow
#     }

    if((it<=rsit)){
        fprintf(fn,'RS - It. : %g   ',it)
        fprintf(fn,'OFV : %g\n',dmf)
        fprintf('RS - It. : %g   ',it)
        fprintf('OFV : %g\n',dmf)
    } else {
        fprintf('SG - It. : %g',it-rsit)
        if((poped.db$convergence_eps!=0)){
            fprintf('  OFV : %5.4g   Diff. : %5.4g\n',dmf,diff)
        } else {
            fprintf('  OFV : %5.4g\n',dmf)
        }
    }
    if((it==rsit)){
      fprintf(fn,'\n*******************************\nR.S. Results\n ')
      fprintf(fn,'OFV(mf) = %g\n\n',dmf)
      if((opt_xt==TRUE)){
        print_xt(xtopt,poped.db$gni,poped.db$global_model_switch,fn,
                 head_txt="Optimized Sampling Schedule\n")
        
      }
      if((opt_x==TRUE)){
        tmp_txt <- "\nOptimized Discrete Variables"
        tmp_txt <- paste(tmp_txt,':\n',sep="")
        fprintf(fn,tmp_txt)
        for(ct1 in 1:poped.db$m){
          fprintf(fn,'Group %g: ', ct1)
          for(ct2 in 1:poped.db$nx){
            tmp_txt <- '%g'
            if(ct2<poped.db$nx) tmp_txt <- paste(tmp_txt,' : ',sep="")        
            fprintf(fn,tmp_txt,xopt[ct1,ct2])
            #fprintf(tmp_txt,xopt[ct1,ct2])
          }
          fprintf(fn,'\n')
        }
        fprintf(fn,'\n')
      }
      if((opt_a==TRUE)){
        tmp_txt <- "\nOptimized Covariates"
        tmp_txt <- paste(tmp_txt,':\n',sep="")
        fprintf(fn,tmp_txt)
        for(ct1 in 1:poped.db$m){
          fprintf(fn,'Group %g: ', ct1)
          for(ct2 in 1:poped.db$na){
            tmp_txt <- '%g'
            if(ct2<poped.db$na) tmp_txt <- paste(tmp_txt,' : ',sep="")
            fprintf(fn,tmp_txt,aopt[ct1,ct2])
            #fprintf(tmp_txt,aopt[ct1,ct2])
          }
          fprintf(fn,'\n')
        }
        fprintf(fn,'\n')
      }
      fprintf(fn,'*********************************\n\n')
      
    }
    if((it>rsit)){
        if((poped.db$use_logfile==TRUE || abs(diff)<poped.db$convergence_eps || it==rsit+poped.db$sgit)){
            fprintf(fn,'S$G. - Iteration %g ----------------------------------\n',it-rsit)
            if((it==rsit+poped.db$sgit || abs(diff)<=poped.db$convergence_eps)){
                fprintf(fn,'FINAL:********************************************\n')
            }
            if((poped.db$optsw[2]==TRUE)){
                if((ixt==TRUE)){
                    fprintf(fn,'xt Gradient Inversion Occured\n')
                }
                fprintf(fn,'Grad(OFV(mf)) xt:\n')
                writet(fn,ni,gxt)
                fprintf(fn,'xt opt:\n')
                writet(fn,ni,xtopt)
            }
            if((poped.db$optsw[4]==TRUE)){
                if((ia==TRUE)){
                    fprintf(fn,'a Gradient Inversion Occured\n')
                }
                fprintf(fn,'Grad(OFV(mf)) a:\n')
                write_matrix(fn,ga)
                fprintf(fn,'aopt:\n')
                write_matrix(fn,aopt)
            }
            fprintf(fn,'OFV(mf)    : %g\n',dmf)
            fprintf(fn,'diff       : %g\n',diff)
            if(((it==rsit+poped.db$sgit) || abs(diff)<=poped.db$convergence_eps)){
                fprintf(fn,'*************************************************************\n')
            }
        }
    }
}
return( ) 
}

