## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

Dtrace <- function(fn,it,ni,xtopt,xopt,aopt,gxt,ga,dmf,diff,ixt,ia,itvector,dmfvector,globalStructure,
                   opt_xt=globalStructure$optsw[2],
                   opt_a=globalStructure$optsw[4],opt_x=globalStructure$optsw[3],
                   opt_samps=globalStructure$optsw[1],opt_inds=globalStructure$optsw[5],
                   rsit=globalStructure$rsit){


if((it==0)){
    fprintf(fn,'*******************************\nInitial Value\n ')
    fprintf(fn,'OFV(mf) = %g\n',dmf)
    fprintf(fn,'*******************************\n\n')
}
if((it!=0)){
#     if((globalStructure$bShowGraphs==TRUE)){
#         if((globalStructure$Engine$Type==1)){
#             #set(0,'CurrentFigure',1)
#         } else {
#             #figure(1)
#         }
#         clf
#         optSum = globalStructure$optsw[2]+globalStructure$optsw[3]+globalStructure$optsw[4]
#         numRows = 1
#         if((optSum>1)){
#             numRows = 2
#         }
# 
#         if((globalStructure$optsw[2])){
#             subplot(numRows,2,1)
#             title('The current sampling times for each group')
#             xlabel('Sampling time')
#             ylabel('Group nr')
#             ##hold on
#             for(i in 1:globalStructure$m){
#                 plot(xtopt[i,1:globalStructure$gni[i]],matrix(1,1,globalStructure$gni[i])*i,'b*','linestyle','none')
#             }
#             ##hold off
#         }
#         if((globalStructure$optsw[3]==TRUE)){
#             subplot(numRows,2,1+globalStructure$optsw[2])
#             title('The current discrete var. for each group')
#             xlabel('Discrete var.-value')
#             ylabel('Group nr')
#             ##hold on
#             for(i in 1:globalStructure$m){
#                 plot(xopt(i,),matrix(1,1,size(xopt,2))*i,'b*')
#             }
#             ##hold off
#         }
#         if((globalStructure$optsw[4]==TRUE)){
#             subplot(numRows,2,1+globalStructure$optsw[2]+globalStructure$optsw[3])
#             title('The current covariates for each group')
#             xlabel('Covariate-value')
#             ylabel('Group nr')
#             ##hold on
#             for(i in 1:globalStructure$m){
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
        if((globalStructure$convergence_eps!=0)){
            fprintf('  OFV : %5.4g   Diff. : %5.4g\n',dmf,diff)
        } else {
            fprintf('  OFV : %5.4g\n',dmf)
        }
    }
    if((it==rsit)){
      fprintf(fn,'\n*******************************\nR.S. Results\n ')
      fprintf(fn,'OFV(mf) = %g\n\n',dmf)
      if((opt_xt==TRUE)){
        print_xt(xtopt,globalStructure$gni,poped.db$global_model_switch,fn,
                 head_txt="Optimized Sampling Schedule\n")
        
      }
      if((opt_x==TRUE)){
        tmp_txt <- "\nOptimized Discrete Variables"
        tmp_txt <- paste(tmp_txt,':\n',sep="")
        fprintf(fn,tmp_txt)
        for(ct1 in 1:globalStructure$m){
          fprintf(fn,'Group %g: ', ct1)
          for(ct2 in 1:globalStructure$nx){
            tmp_txt <- '%g'
            if(ct2<globalStructure$nx) tmp_txt <- paste(tmp_txt,' : ',sep="")        
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
        for(ct1 in 1:globalStructure$m){
          fprintf(fn,'Group %g: ', ct1)
          for(ct2 in 1:globalStructure$na){
            tmp_txt <- '%g'
            if(ct2<globalStructure$na) tmp_txt <- paste(tmp_txt,' : ',sep="")
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
        if((globalStructure$use_logfile==TRUE || abs(diff)<globalStructure$convergence_eps || it==rsit+globalStructure$sgit)){
            fprintf(fn,'S$G. - Iteration %g ----------------------------------\n',it-rsit)
            if((it==rsit+globalStructure$sgit || abs(diff)<=globalStructure$convergence_eps)){
                fprintf(fn,'FINAL:********************************************\n')
            }
            if((globalStructure$optsw[2]==TRUE)){
                if((ixt==TRUE)){
                    fprintf(fn,'xt Gradient Inversion Occured\n')
                }
                fprintf(fn,'Grad(OFV(mf)) xt:\n')
                writet(fn,ni,gxt)
                fprintf(fn,'xt opt:\n')
                writet(fn,ni,xtopt)
            }
            if((globalStructure$optsw[4]==TRUE)){
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
            if(((it==rsit+globalStructure$sgit) || abs(diff)<=globalStructure$convergence_eps)){
                fprintf(fn,'*************************************************************\n')
            }
        }
    }
}
return( ) 
}

