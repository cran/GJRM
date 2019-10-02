bcorrec <- function(VC, params){ #, lB, uB){

m2sel <- c("WEI","iG","GA","BE","FISK","DAGUM","TW","SM","GP","GPII","GPo")

###############################
# fixed quantities
###############################

n      <- VC$n
rc     <- VC$rc
margin <- VC$margins[2]

weights <- VC$weights

if(is.null(VC$X2)){VC$X2 <- VC$X3 <- matrix(1, n, 1); VC$X2.d2 <- VC$X3.d2 <- 1} 

###################################################################################################
###################################################################################################

if(is.null(VC$lB) && is.null(VC$uB)){

if( margin %in% c("N","GU","rGU","LO") )                                  { lB <- -Inf;      uB <- Inf}
if( margin %in% c("LN","WEI","iG","GA","DAGUM","SM","FISK","GP","GPII","GPo","TW")  )       { lB <- sqrt(.Machine$double.eps); uB <- Inf} # tw should be zero here?
if( margin %in% c("BE")  )                                                { lB <- sqrt(.Machine$double.eps); uB <- 0.9999999}

}else{lB <- VC$lB; uB <- VC$uB}


###################################################################################################
###################################################################################################

eta <- eta.tr(VC$X1%*%params[1:VC$X1.d2], margin)
ss  <- esp.tr(VC$X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)], margin)

sigma2    <- ss$vrb
sigma2.st <- ss$vrb.st

if( margin %in% c("DAGUM","SM","TW") ){

            nus   <- enu.tr(VC$X3%*%params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)], margin)
            nu    <- nus$vrb
            nu.st <- nus$vrb.st
            
} else nu <- nu.st <- 1


###################################################################################################
###################################################################################################

#intB1   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) integrate(intB, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)$value              


intB1   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(intB, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)            
v.intB1 <- Vectorize(intB1) 
int1    <- v.intB1(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
b       <- sum(weights) - exp(-rc)*sum(weights*int1)

###################################################################################################
###################################################################################################

if( margin %in% m2sel ){

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(gradBbit1, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
gradbit1  <- -colSums( c(weights*intgrad)*VC$X1 ) # minus comes from theoretical derivation

}

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(gradBbit2, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
gradbit2  <- -colSums( c(weights*intgrad)*VC$X2 )

if( margin %in% c("DAGUM","SM","TW") ){

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(gradBbit3, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
gradbit3  <- -colSums( c(weights*intgrad)*VC$X3 )

}

###################################################################################################
###################################################################################################

if( margin %in% m2sel ){

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(hessBbit1, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
hessbit1  <- - crossprod(VC$X1*c(weights*intgrad), VC$X1)

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(hessBbit3, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
hessbit3  <- - crossprod(VC$X1*c(weights*intgrad), VC$X2)

}

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(hessBbit2, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)             
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
hessbit2  <- - crossprod(VC$X2*c(weights*intgrad), VC$X2)

if( margin %in% c("DAGUM","SM","TW") ){

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(hessBbit4, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
hessbit4  <- - crossprod(VC$X3*c(weights*intgrad), VC$X3)

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(hessBbit5, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
hessbit5  <- - crossprod(VC$X1*c(weights*intgrad), VC$X3)

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(hessBbit6, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)
hessbit6  <- - crossprod(VC$X2*c(weights*intgrad), VC$X3)

}

###################################################################################################
###################################################################################################

if(margin %in% VC$m2){
   if(margin %in% m2sel) bp <- c(gradbit1, gradbit2) else bp <- c(rep(0,length(params[1:VC$X1.d2])), gradbit2) 
}

if(margin %in% VC$m3)    bp <- c(gradbit1, gradbit2, gradbit3) 


###################################################################################################

if(margin %in% VC$m2){

   if(margin %in% m2sel) bs <- rbind( cbind( hessbit1   , hessbit3  ), 
                                      cbind( t(hessbit3), hessbit2  ) ) else{
                                      
                                            lb <- length( c( 1:VC$X1.d2 ) )
                                            ls <- length( c( (1+VC$X1.d2):(VC$X1.d2+VC$X2.d2) ) )

                                            bs <- rbind( cbind( matrix(0,lb,lb), matrix(0,lb,ls)  ), 
                                                         cbind( t(matrix(0,lb,ls)), hessbit2  ) ) 
                                                                          
                                                                             }
}

if(margin %in% VC$m3) bs <- rbind( cbind( hessbit1   , hessbit3   , hessbit5 ), 
                                   cbind( t(hessbit3), hessbit2   , hessbit6 ),
                                   cbind( t(hessbit5), t(hessbit6), hessbit4 ) ) 

###################################################################################################
###################################################################################################

list(b = b, bp = bp, bs = bs)

}    

