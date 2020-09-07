working.comp <- function(fit){
   
  G <- -(fit$gradient - fit$S.h2)
  H <-   fit$hessian  - fit$S.h
  
  # H <- PDef(H)$res
  
  tolH <- sqrt(.Machine$double.eps)
  
  W.eig <- eigen(H, symmetric = TRUE)
  
  if(min(W.eig$values) < tolH && sign( min( sign(W.eig$values) ) ) == -1) W.eig$values <- abs(W.eig$values)  
  if(min(W.eig$values) < tolH ) { pep <- which(W.eig$values < tolH); W.eig$values[pep] <- tolH } 
  
  # seems a reasonable fix, should it be required, and it worked well so far in a variety of situations
  # we could also do H <- PDef(H)$res to have something that feels a bit more theoretically founded
  
  srev    <- sqrt(W.eig$val)
  c.W     <- W.eig$vec%*%tcrossprod(diag(  srev, nrow = length(srev), ncol = length(srev)), W.eig$vec) 
  W.invsr <- W.eig$vec%*%tcrossprod(diag(1/srev, nrow = length(srev), ncol = length(srev)), W.eig$vec)
  
  X <- c.W 
  Z <- W.invsr%*%G + X%*%fit$argument 
   
  rm(G, H, W.eig, srev, c.W, W.invsr)        
          
 list( X = X , Z = Z )

}


























