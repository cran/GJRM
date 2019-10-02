cv.inform <- function(x, K = 5, data, informative = "yes"){          
        
if(x$VC$informative == "no") stop("This function can be used only for informative survival models.")        
        
Slpar <- NA              
              
mDat  <- data[sample(nrow(data)), ]

folds <- cut(seq(1,nrow(mDat)), breaks = K, labels = FALSE)

for(i in 1:K){

    test.ind    <- which(folds != i) 
    test.mDat   <- mDat[test.ind, ]
    train.mDat  <- mDat[-test.ind, ]
    cens.test   <- x$cens[test.ind] 
    cens.train  <- x$cens[-test.ind] 
    test.mDat$cens.test <- cens.test 
    


if(informative == "yes"){

out <- gamlss(x$formula, data = test.mDat, 
              surv = TRUE, margin = x$margins[1], margin2 = x$margins[2], cens = cens.test, 
              informative = "yes", inform.cov = c(x$VC$inform.cov))
              
params1 <- out$fit$params1
params2 <- out$fit$params2

Xd <- try(Xdpred(out$gam1, train.mDat, out$VC$v1pred), silent = TRUE)
if(class(Xd) == "try-error") stop("Your factor variable(s) has(have) different levels in the test and training data sets.")

Xp <- predict(out$gam1, type = "lpmatrix", newdata = train.mDat)

}

if(informative == "no"){

out1 <- gamlss(list(x$formula[[1]]), data = test.mDat, surv = TRUE, margin = x$margins[1], cens = cens.test)
out2 <- gamlss(list(x$formula[[2]]), data = test.mDat, surv = TRUE, margin = x$margins[2], cens = 1 - cens.test)   

params1 <- out1$fit$params1
params2 <- out2$fit$params1

Xd <- try(Xdpred(out1$gam1, train.mDat, out1$VC$v1pred), silent = TRUE)
if(class(Xd) == "try-error") stop("Your factor variable(s) has(have) different levels in the test and training data sets.")

Xp <- predict(out1$gam1, type = "lpmatrix", newdata = train.mDat)

}



eta1 <- Xp%*%params1
eta2 <- Xp%*%params2

Xd1P <- Xd%*%params1 
Xd1P <- ifelse(Xd1P < sqrt(.Machine$double.eps), sqrt(.Machine$double.eps), Xd1P ) # safety check
        
Xd2P <- Xd%*%params2
Xd2P <- ifelse(Xd2P < sqrt(.Machine$double.eps), sqrt(.Machine$double.eps), Xd2P )
 
pd1 <- probmS(eta1, x$VC$margins[1])
pd2 <- probmS(eta2, x$VC$margins[2])
  
p1 <- pd1$pr
p2 <- pd2$pr
  
dS1eta1 <- pd1$dS
dS2eta2 <- pd2$dS
  
Slpar[i] <- sum(  log(p1) + log(p2) + cens.train*log( Xd1P*-dS1eta1/p1) + (1-cens.train)*log( Xd2P*-dS2eta2/p2)  )

print(i)

}


list(sl = sum(Slpar), Slpar = Slpar)


}