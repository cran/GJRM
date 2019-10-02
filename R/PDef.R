PDef <- function(omega){

  eS    <- eigen(omega, symmetric = TRUE)                
  e.val <- eS$values
  e.vec <- eS$vectors
  
  check.eigen <- any(e.val <= 0)

  if(check.eigen == TRUE){
  
    n.e.val <- e.val[e.val <= 0]
    s <- sum(e.val[n.e.val])*2  
    t <- s^2*100 + 1
    
    p <- min(e.val[(e.val <= 0) == FALSE])
    
    e.val[e.val <= 0] <- p*(s - n.e.val)^2/t
    
    D     <- diag(e.val)
    D.inv <- diag(1/e.val)

    res     <- e.vec %*% D     %*% t(e.vec) 
    res.inv <- e.vec %*% D.inv %*% t(e.vec)
      
  } else {res <- omega; res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec)} 
  
  
res.inv <- (res.inv + t(res.inv) ) / 2 
 
list(check.eigen = check.eigen, res = res, res.inv = res.inv)

} 


