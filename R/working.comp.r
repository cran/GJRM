working.comp <- function(fit){
   
  G <- -(fit$gradient - fit$S.h2)
  H <- fit$hessian - fit$S.h
  
  W.eig <- eigen(H, symmetric=TRUE)
  
  if(min(W.eig$values) < sqrt(.Machine$double.eps) && sign( min( sign(W.eig$values) ) ) == -1) W.eig$values <- abs(W.eig$values)  
  if(min(W.eig$values) < sqrt(.Machine$double.eps) ) { pep <- which(W.eig$values < sqrt(.Machine$double.eps)); W.eig$values[pep] <- 0.0000001 }  

  srev    <- sqrt(W.eig$val)
  c.W     <- W.eig$vec%*%tcrossprod(diag(srev, nrow = length(srev), ncol = length(srev) )  ,W.eig$vec) 
  W.invsr <- W.eig$vec%*%tcrossprod(diag(1/srev, nrow = length(srev), ncol = length(srev)),W.eig$vec)
  
  X <- c.W 
  Z <- W.invsr%*%G + X%*%fit$argument 
   
  rm(G, H, W.eig, srev, c.W, W.invsr)        
          
 list( X = X , Z = Z )

}


























