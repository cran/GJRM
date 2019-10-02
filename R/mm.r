mm <- function(ob){

  epsilon <- sqrt(.Machine$double.eps)
  max.p   <- 0.9999999

  res <- ifelse(ob > max.p, max.p, ob) 
  res <- ifelse(res < epsilon, epsilon, res) 

  res
      
}


