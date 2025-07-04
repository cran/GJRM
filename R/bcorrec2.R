bcorrec2 <- function(VC, params){

###############################
# fixed quantities
###############################

n      <- VC$n
rc     <- VC$rc
margin <- VC$margins[2]
min.dn <- 1e-160 # smallest possible allowed VC$min.dn 
min.pr <- VC$min.pr
max.pr <- VC$max.pr 
nu <- nu.st <- 1

int1 <- NA; G <- matrix(NA, n, length(params))

if(is.null(VC$X2)){VC$X2 <- VC$X3 <- matrix(1, n, 1); VC$X2.d2 <- VC$X3.d2 <- 1} 

###################################################################################################

if(is.null(VC$lB) && is.null(VC$uB)){

if( margin %in% c("N","tN","GU","rGU","LO","LN") )                                     { lB <- -Inf;      uB <- Inf}
if( margin %in% c("tN"))                                                                 lB <- VC$left.trunc
if( margin %in% c("WEI","IG","GA","DAGUM","SM","FISK","GP","GPII","GPo","TW")  )  { lB <- sqrt(.Machine$double.eps); uB <- Inf} # tw should be zero here?
if( margin %in% c("BE")  )                                                        { lB <- sqrt(.Machine$double.eps); uB <- 0.999999}

}else{lB <- VC$lB; uB <- VC$uB}

###################################################################################################


for(i in 1:n){

  X1 <- VC$X1[i,]
  X2 <- VC$X2[i,]
  X3 <- VC$X3[i,]
  
  eta       <- X1%*%params[1:VC$X1.d2]
  sigma2    <- esp.tr(X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)], margin)$vrb
  sigma2.st <- X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)]
  
  if( margin %in% c("DAGUM","SM","TW") ){
  nu.st <- X3%*%params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  nu    <- enu.tr(X3%*%params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)], margin)$vrb
  }

  int1[i]  <- integrate(int1f, lB, uB, eta, sigma2, nu, margin, rc, min.dn, min.pr, max.pr, left.trunc = VC$left.trunc)$value

    for(j in 1:length(params)) G[i, j] <- integrate(d.bpsi, lB, uB, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j, min.dn, min.pr, max.pr, left.trunc = VC$left.trunc)$value

}

b  <- n - exp(-rc)*sum(int1)
bp <- -colSums(G) 
#bs <- -jacobian(func = gradF, x = params, method = "simple", n = n, VC = VC, margin = margin, lB = lB, uB = uB, rc = rc)
bs <- NULL 

list(b = b, bp = bp, bs = bs)


}    