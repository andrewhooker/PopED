grouped_rand <- function(G,xtopt,dxt,ff1,axt){
# Calculate a new random search xt in a grouped xt manner way.
size1 = size(xtopt,1)
size2 = size(xtopt,2)
xt = zeros(size1,size2)
for(k in 1:max(max(G))){
    random_num = randn 
    tmp = matrix(1,size1,size2)*k
    inters = (G==tmp)
    xt = (xtopt+dxt/ff1*random_num*(axt>0))*(inters==1)
    
#     for(i in 1:size1){
#        for(j in 1:size2){
#           if (inters[i,j]==1)
#              xt[i,j] = xtopt[i,j]+dxt(i,j)/ff1*random_num*(axt[i,j]>0)
#           }
#        }
#     }
 }
return( xt) 
}
