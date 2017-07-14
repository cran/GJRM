TRIapprox <- function(eta1, eta2, eta3, Sigma){
  
be <- data.frame(eta1, eta2, eta3)

C <- chol(Sigma) 
C <- matrix( c( C[1, 1],  C[1, 2], C[1, 3],
                0,  C[2, 2], C[2, 3],
                0,        0, C[3, 3] ), 3 , 3)


min.p1 <- min( pnorm(eta1/sqrt(Sigma[1, 1]))  )
min.p2 <- min( pnorm(eta2/sqrt(Sigma[2, 2]))  )
min.p3 <- min( pnorm(eta3/sqrt(Sigma[3, 3]))  )

if(min(min.p1, min.p2, min.p3) == min(min.p1))  choose.i <- 1
if(min(min.p1, min.p2, min.p3) == min(min.p2))  choose.i <- 2
if(min(min.p1, min.p2, min.p3) == min(min.p3))  choose.i <- 3

beNew1 <- be
beNew1[, 1] <- be[, choose.i]
beNew2 <- beNew1
beNew2[, choose.i] <-  be[, 1]
be <- beNew2

SigmaNew1 <- Sigma
SigmaNew1[1, ] <- Sigma[choose.i, ]
SigmaNew1[choose.i, ] <- Sigma[1, ]
SigmaNew2 <- SigmaNew1
SigmaNew2[, 1] <- SigmaNew1[, choose.i]; SigmaNew2[, choose.i] <- SigmaNew1[, 1]
Sigma <- SigmaNew2

C[1, 1] <- sqrt(Sigma[1, 1])
C[2, 1] <- Sigma[2, 1]/C[1, 1]
C[3, 1] <- Sigma[3, 1]/C[1, 1]

b1.hat <- be[, 1]/C[1, 1]

mean.b1 <- -dnorm( b1.hat )/pnorm( b1.hat )
min2.p2 <- min( pnorm( (be[, 2] - C[2, 1] * mean.b1)/sqrt(Sigma[2, 2] - C[2, 1]^2))  )
min2.p3 <- min( pnorm( (be[, 3] - C[3, 1] * mean.b1)/sqrt(Sigma[3, 3] - C[3, 1]^2))  )

choose2.i <- ifelse(min(min2.p2, min2.p3) == min(min2.p2), 2, 3) 

beNew1 <- be
beNew1[, 2] <- be[, choose2.i]
beNew2 <- beNew1
beNew2[, choose2.i] <-  be[, 2]
be <- beNew2

SigmaNew1 <- Sigma
SigmaNew1[2, ] <- Sigma[choose2.i, ]
SigmaNew1[choose2.i, ] <- Sigma[2, ]
SigmaNew2 <- SigmaNew1
SigmaNew2[, 2] <- SigmaNew1[, choose2.i]; SigmaNew2[, choose2.i] <- SigmaNew1[, 2]
Sigma <- SigmaNew2

CNew1 <- C
CNew1[1, 2] <- C[choose2.i, 2]
CNew2 <- CNew1
CNew2[choose2.i, 2] <-  C[1, 2]
C <- CNew2

C[2, 2] <- sqrt(Sigma[2, 2] - C[2, 1]^2)
C[3, 2] <- (Sigma[3, 2] - C[3, 1] * C[2, 1])/C[2, 2]

b2.hat <- ( be[, 2] - C[2, 1] * mean.b1 )/C[2, 2]

mean.b2 <- -dnorm( b2.hat )/pnorm( b2.hat )

C[3, 3] <- sqrt( Sigma[3, 3] - C[3, 1]^2 - C[3, 2]^2 )

b3.hat <- ( be[, 3] - C[3, 1] * mean.b1 - C[3, 2] * mean.b2 )/C[3, 3]

C.D <- matrix( c( C[1, 1],  C[2, 1], 0,
                  0,  C[2, 2], 0,
                  0,        0, C[3, 3] ), 3 , 3)

DSigma <- C.D %*% t(C.D)
LSigma <- C %*% solve(C.D)


b1.hat2 <- be[, 1]/DSigma[1, 1]
b2.hat2 <- be[, 2]/DSigma[2, 2]

rho1 <- DSigma[1, 2]/sqrt( DSigma[1, 1] * DSigma[2, 2] )
P1   <- pbinorm(b1.hat2, b2.hat2, cov12 = rho1)

q1 <- sqrt(1 - rho1^2)

mu1 <- 1/P1 * ( - rho1 * dnorm(b2.hat2) * pnorm( (b1.hat2 - rho1 * b2.hat2)/q1 ) - dnorm(b1.hat2) * pnorm( (b2.hat2 - rho1 * b1.hat2)/q1 ) )
mu2 <- 1/P1 * ( - rho1 * dnorm(b1.hat2) * pnorm( (b2.hat2 - rho1 * b1.hat2)/q1 ) - dnorm(b2.hat2) * pnorm( (b1.hat2 - rho1 * b2.hat2)/q1 ) )

e1 <- mu1 * sqrt(DSigma[1, 1])
e2 <- mu2 * sqrt(DSigma[2, 2]) 

g3 <- LSigma[3, 1] * e1 + LSigma[3, 2] * e2

b3.hat2 <- ( be[, 3] - g3 )/sqrt(DSigma[3, 3])

p111 <- P1 * pnorm(b3.hat2)

p111

}
