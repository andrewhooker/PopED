## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

get_vector_str <- function(vec){
   str_vector = '' 
   if((length(vec)==1)){
       str_vector = sprintf('%d',vec[1])
       return(str_vector)
   } else {
       for(i in 1:length(vec)){
           str_vector = sprintf('%s %d ',str_vector,vec[i])
       }
       str_vector = sprintf('%s',str_vector)
   }
return( str_vector) 
}
