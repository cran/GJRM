PosDefCor <- function(Sigma, Chol = FALSE, theta12.st, theta13.st, theta23.st){



if(Chol == FALSE){
  
  eS <- eigen(Sigma)                 
  check.eigen <- any(eS$values <= 0)
  
  if(check.eigen == TRUE){
    
    D.dash <- diag(abs(eS$values), 3, 3)
    P      <- eS$vectors
    R.dash <- P %*% D.dash %*% t(P)
    D1 <- diag(1/sqrt( diag(R.dash) ), 3, 3)
    Res <- D1 %*% R.dash %*% D1
  } else Res <- Sigma 
  

} 


if(Chol == TRUE){

 C <- matrix(c(1, 0, 0, theta12.st, 1, 0, theta13.st, theta23.st, 1), nrow=3, byrow=TRUE)
 Sigma.star <- C%*%t(C)
 T <- diag(1/sqrt(diag(Sigma.star)))
 Res <- T %*% Sigma.star %*% T

}

Res
 
}
