concentrat <- function (t, y, parms) {
  with(as.list(y), {
    dU <- (-Sm*U/(Km+U))  
    dC <- (-k*A*C*f)/(Henry_k*V) + 2*Sm*U/(Km+U)
    
    list(c(dU, dC))})}


 