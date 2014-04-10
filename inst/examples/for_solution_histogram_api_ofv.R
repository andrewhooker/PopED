tmp <- c()
for(i in 1:100){
  output <- evaluate.e.ofv.fim(poped.db,ofv_calc_type=4)
  tmp <- c(tmp,output$E_ofv)
}
API_OFV <- tmp
hist(API_OFV)