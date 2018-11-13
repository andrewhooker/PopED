## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

grad_all <- function(func,select_par,nRow,...,subset=NULL,currentOcc=NULL,noPopED=FALSE,offdiag=FALSE){
  arg_list = list(...)
  def0 = arg_list[[select_par]]
  if (is.null(currentOcc)) {
    idx = idx0 = seq_along(def0)
    if (!is.null(subset)) {
      idx0 = as.vector(cumsum(subset)*subset)
      idx = seq(max(idx0))
    }
    if (is.matrix(def0) & select_par>6 & offdiag==FALSE) {
      idx0 = diag(idx0)
    }
    if (is.matrix(def0) & select_par>6 & offdiag==TRUE) {
      tmp = 0*def0
      tmp[lower.tri(tmp)] = idx0
      idx0 = tmp + t(tmp)
    }
  } else {
    idx  = seq(size(def0,1))
    idx0 = zeros(size(def0))
    idx0[,currentOcc] = idx
  }

  poped.db = arg_list[[length(arg_list)]]
  hlf  = poped.db$settings$hlf
  grad_all_switch = poped.db$settings$grad_all_switch[1]

  if (noPopED == TRUE) arg_list = arg_list[-length(arg_list)]
  
  gradX = zeros(nRow, length(idx))

  #Central approximation
  if (grad_all_switch == 1) {
    for (i in idx) {
      arg_list[[select_par]] = def0 + (idx0 == i)*hlf
      def_plus <- do.call(func, arg_list)
      arg_list[[select_par]] = def0 - (idx0 == i)*hlf
      def_minus <- do.call(func, arg_list)
      if (noPopED == FALSE) {
        def_plus = def_plus[[1]]        
        def_minus = def_minus[[1]]        
      }
      gradX[,i] = (def_plus - def_minus)/(2.0*hlf)
    }
  } else {
    #Complex approximation
    if (grad_all_switch == 0) {
      for (i in idx) {
        arg_list[[select_par]] = def0 + (idx0 == i)*complex(real = 0, imaginary = hlf)
        def_plus <- do.call(func, arg_list)[[1]]
        gradX[,i] = Im(def_plus)/hlf
      }
    } else {
        stop(sprintf('Unknown derivative option for grad_all'))
    }
  }
  return(gradX) 
}
