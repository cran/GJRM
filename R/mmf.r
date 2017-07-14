mmf <- function(ob){
  
  sob   <- sign(ob)
  max.p <- 0.9999999
  ob    <- abs(ob)

  res <- ifelse(ob > max.p, max.p, ob) 
  res <- res*sob

  res
      
}


