concentrat_msrd <- function (t, y, parms) {
  with(as.list(y), {
    dU <- (-Sm*U/(Km+U))  
    dC <- (-k_msrd*A*C*f_msrd)/(Henry_k_msrd*V) + 2*Sm*U/(Km+U)
    
    list(c(dU, dC))})}