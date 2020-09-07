mmf <- function(ob, max.pr){

# 0.999999 could have one more 9 (or add 1e-06) but then it would probably be too close
 
  sob   <- sign(ob)
  ob    <- abs(ob)

  res <- ifelse(ob > max.pr, max.pr, ob) 
  res <- res*sob

  res
      
}


