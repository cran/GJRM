bcorrec2 <- function(VC){

###############################
# fixed quantities
###############################

n      <- VC$n
rc     <- VC$rc
margin <- VC$margins[2]

int1 <- NA; G <- matrix(NA, n, length(VC$params))

if(is.null(VC$X2)){VC$X2 <- VC$X3 <- matrix(1, n, 1); VC$X2.d2 <- VC$X3.d2 <- 1} 

if(is.null(VC$lB) && is.null(VC$uB)){

if( margin %in% c("N","N2","GU","rGU","LO","LN") )             { lB <- -Inf;      uB <- Inf}
if( margin %in% c("WEI","iG","GA","DAGUM","SM","FISK")  )      { lB <- 0.0000001; uB <- Inf}
if( margin %in% c("BE")  )                                     { lB <- 0.0000001; uB <- 0.9999999}

}else{lB <- VC$lB; uB <- VC$uB} 

###################################################################################################
###################################################################################################
int1f <- function(y){ 
   pdf <- distrHs(y, eta, sigma2, 1, nu, 1, margin)$pdf2
   log( 1 + exp( log( pdf ) + rc ) )
}
###################################################################################################
###################################################################################################
d.bpsi <- function(y, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j){ 

       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)

       pdf2                 <- dHs$pdf2
       derpdf2.dereta2      <- dHs$derpdf2.dereta2 
       derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st  
       derpdf2.dernu.st     <- dHs$derpdf2.dernu.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dbe       <- derpdf2.dereta2/pdf2
       dl.dsigma.st <- derpdf2.dersigma2.st/pdf2
       dl.dnu.st    <- derpdf2.dernu.st/pdf2


       if( margin %in% c("DAGUM","SM") ) res <- cbind( comp2*as.numeric(dl.dbe)%*%t(X1), comp2*as.numeric(dl.dsigma.st)%*%t(X2), comp2*as.numeric(dl.dnu.st)%*%t(X3) ) else
                                         res <- cbind( comp2*as.numeric(dl.dbe)%*%t(X1), comp2*as.numeric(dl.dsigma.st)%*%t(X2) )
      
       res[, j]
}
###################################################################################################
###################################################################################################
gradF <- function(params){

  G <- matrix(NA, n, length(params))

for(i in 1:n){

  X1 <- VC$X1[i,]
  X2 <- VC$X2[i,]
  X3 <- VC$X3[i,]
  nu <- nu.st <- 1
  
  eta       <- X1%*%params[1:VC$X1.d2]
  sigma2.st <- X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)]
  
     if( margin %in% c("DAGUM","SM") ){ 
  
       nu.st <- X3%*%params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)] 
 
       ss <- enu.tr(nu.st, margin)  
       nu.st <- ss$vrb.st 
       nu    <- ss$vrb   
                                      }
 
  ss        <- esp.tr(sigma2.st, margin)
  sigma2.st <- ss$vrb.st
  sigma2    <- ss$vrb
                                    
  for(j in 1:length(params)) G[i, j] <- integrate(d.bpsi, lB, uB, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j)$value

             }
             
  colSums(G)
  
}
########################################################################################################################
########################################################################################################################

for(i in 1:n){

  X1 <- VC$X1[i,]
  X2 <- VC$X2[i,]
  X3 <- VC$X3[i,]
  
  eta       <- X1%*%VC$params[1:VC$X1.d2]
  sigma2    <- esp.tr(X2%*%VC$params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)], margin)$vrb
  sigma2.st <- X2%*%VC$params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)]
  
  if( margin %in% c("DAGUM","SM") ){
  nu.st <- X3%*%VC$params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  nu    <- enu.tr(X3%*%VC$params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)], margin)$vrb
  }

  int1[i]   <- integrate(int1f, lB, uB)$value

    for(j in 1:length(VC$params)) G[i, j] <- integrate(d.bpsi, lB, uB, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j)$value

}

b  <- n - exp(-rc)*sum(int1)
bp <- -colSums(G) 
#bs <- -jacobian(gradF, VC$params, method = "simple")
 
list(b = b, bp = bp)#, bs = bs)


}    