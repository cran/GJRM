bcorrecDiscr <- function(VC, params){

m1d  <- c("PO", "ZTP")
m2d  <- c("NBI", "NBII","PIG","DGP","DGPII")

# below there are operations that can be 
# simplified based on whether sigma2 has regressors or not
# and the set up of the matrices

###############################
# fixed quantities
###############################

n      <- VC$n
rc     <- VC$rc
margin <- VC$margins[2]
ym     <- VC$y1m
y      <- VC$y
ygrid  <- VC$ygrid
chs    <- VC$chunk.size

weights <- VC$weights


no.splits <- ceiling(length(ygrid)*n/chs)
if(no.splits > 1) num.ind <- cut(1:n, no.splits, labels = FALSE)

###################################################################################################

if(is.null(VC$X2)){VC$X2 <- VC$X3 <- matrix(1, n, 1); VC$X2.d2 <- VC$X3.d2 <- 1} 

eta <- eta.tr(VC$X1%*%params[1:VC$X1.d2], margin)

if( !(margin %in% c(m1d)) ){

  ss <- esp.tr(VC$X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)], margin)
  sigma2    <- ss$vrb
  sigma2.st <- ss$vrb.st

} else sigma2 <- sigma2.st <- 1 

if( length(sigma2)==1){ sigma2 <- rep(sigma2, n); sigma2.st <- rep(sigma2.st, n)} 

###################################################################################################


if(no.splits > 1){

sumb <- 0
gradbit1b <- rep(0, dim(VC$X1)[2])
gradbit2b <- rep(0, dim(VC$X2)[2])
hessbit1b <- matrix(0, nrow = dim(VC$X1)[2], ncol = dim(VC$X1)[2])
hessbit2b <- matrix(0, nrow = dim(VC$X2)[2], ncol = dim(VC$X2)[2])
hessbit3b <- matrix(0, nrow = dim(VC$X1)[2], ncol = dim(VC$X2)[2])

for(i in 1:no.splits){

n.red      <- length(eta[num.ind==i])
etaL       <- rep(eta[num.ind==i], each = length(ygrid))
weightsL   <- rep(weights[num.ind==i], each = length(ygrid))
sigma2L    <- rep(sigma2[num.ind==i], each = length(ygrid))
sigma2.stL <- rep(sigma2.st[num.ind==i], each = length(ygrid))
ygridL     <- rep(ygrid, n.red)

                      X1L <- apply(as.matrix(VC$X1[num.ind==i,]), 2, rep, each = length(ygrid))
if( margin %in% m2d ) X2L <- apply(as.matrix(VC$X2[num.ind==i,]), 2, rep, each = length(ygrid))

####

sumb <- sumb + sum(weightsL*intB(ygridL, etaL, sigma2L, 1, 1, 1, margin, rc, discr = TRUE, ym))

gradbit1b <- gradbit1b - colSums( c( weightsL*gradBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X1L )
hessbit1b <- hessbit1b - crossprod(X1L*c(weightsL*hessBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym)), X1L)

  if( margin %in% m2d ){

     gradbit2b <- gradbit2b - colSums(  c( weightsL*gradBbit2(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X2L )
     hessbit2b <- hessbit2b - crossprod(X2L*c(weightsL*hessBbit2(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym)), X2L)
     hessbit3b <- hessbit3b - crossprod(X1L*c(weightsL*hessBbit3(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym)), X2L)
                
                       }
                     
                     } 
                     
                     
}else{

###########################

#if(no.splits == 0){

etaL       <- rep(eta, each = length(ygrid))
weightsL   <- rep(weights, each = length(ygrid))
sigma2L    <- rep(sigma2, each = length(ygrid))
sigma2.stL <- rep(sigma2.st, each = length(ygrid))
ygridL     <- rep(ygrid, n)

                      X1L <- apply(VC$X1, 2, rep, each = length(ygrid))
if( margin %in% m2d ) X2L <- apply(VC$X2, 2, rep, each = length(ygrid))

####

sumb <- sum(weightsL*intB(ygridL, etaL, sigma2L, 1, 1, 1, margin, rc, discr = TRUE, ym))

gradbit1b <- -colSums( c( weightsL*gradBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X1L )
hessbit1b <- -crossprod(X1L*c(weightsL*hessBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym)), X1L)

if( margin %in% m2d ){

   gradbit2b <- -colSums(  c( weightsL*gradBbit2(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X2L )
   hessbit2b <- -crossprod(X2L*c(weightsL*hessBbit2(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym)), X2L)
   hessbit3b <- -crossprod(X1L*c(weightsL*hessBbit3(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym)), X2L)
                
                     }

}


###################################################################################################
###################################################################################################

b <- sum(weights) - exp(-rc)*sumb

if(margin %in% VC$m2d) bp <- c(gradbit1b, gradbit2b) else bp <- c(gradbit1b)  

if(margin %in% VC$m2d) bs <- rbind( cbind( hessbit1b   , hessbit3b  ), 
                                    cbind( t(hessbit3b), hessbit2b  ) ) else bs <- hessbit1b
                                                                
list(b = b, bp = bp, bs = bs)

}    

