SemiParROY.fit.post <- function(SemiParFit, Model, VC, GAM){

Ve <- R <- X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- eta1S <- eta2S <- theta12 <- theta13 <- edf <- edf1 <- theta12.a <- theta13.a <- NULL
sigma2 <- sigma2.a <- sigma3 <- sigma3.a <- nu2 <- nu2.a <- nu3 <- nu3.a <- NULL
dof2 <- dof2.a <- dof3 <- dof3.a <- nCa1 <- nCa2 <- NULL

cont1par  <- VC$m1d
cont2par  <- c(VC$m2,VC$m2d) 
cont3par  <- VC$m3 
bin.link  <- VC$bl  

if(VC$margins[2] == "LN" || VC$margins[3] == "LN") logLik <- -SemiParFit$fit$l.ln else logLik <- -SemiParFit$fit$l   

pVbres <- postVb(SemiParFit, VC)

He         <- pVbres$He    
Vb.t       <- pVbres$Vb.t

Vb         <- pVbres$Vb        
HeSh       <- pVbres$HeSh      
F          <- pVbres$F         
F1         <- pVbres$F1        
R          <- pVbres$R         
Ve         <- pVbres$Ve        
t.edf      <- pVbres$t.edf     
SemiParFit <- pVbres$SemiParFit

VC1 <- list(BivD = VC$BivD1, BivD2 = NULL)  
VC2 <- list(BivD = VC$BivD2, BivD2 = NULL)  
  
############################################################################################  


SemiParFit$fit$eta2 <- VC$X2s %*% SemiParFit$fit$argument[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
SemiParFit$fit$eta3 <- VC$X3s %*% SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]





if( VC$margins[2] %in% c(cont2par, cont3par) && VC$margins[3] %in% c(cont2par, cont3par) ){

sigma2   <- esp.tr(VC$X4s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)], VC$margins[2])$vrb  
sigma2.a <- mean(sigma2)

sigma3   <- esp.tr(VC$X5s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)], VC$margins[3])$vrb  
sigma3.a <- mean(sigma3)

}



if(VC$margins[2] %in% cont3par && VC$margins[3] %in% cont3par){

nu2   <- enu.tr(VC$X6s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2)], VC$margins[2])$vrb    
nu2.a <- mean(nu2)

nu3   <- enu.tr(VC$X7s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2)], VC$margins[3])$vrb    
nu3.a <- mean(nu3)  
} 
 

############################################################################################


if( VC$margins[2] %in% c(bin.link, cont1par) && VC$margins[3] %in% c(bin.link, cont1par) ){

teta12 <- teta.tr(VC1, VC$X4s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)])$teta
teta13 <- teta.tr(VC2, VC$X5s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)])$teta

ass.msR1 <- ass.ms(VC$BivD1, VC$nCa1, teta12)
ass.msR2 <- ass.ms(VC$BivD2, VC$nCa2, teta13)

}


if( VC$margins[2] %in% c(cont2par) && VC$margins[3] %in% c(cont2par) ){

teta12 <- teta.tr(VC1, VC$X6s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2)])$teta
teta13 <- teta.tr(VC2, VC$X7s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2)])$teta

ass.msR1 <- ass.ms(VC$BivD1, VC$nCa1, teta12)
ass.msR2 <- ass.ms(VC$BivD2, VC$nCa2, teta13)

}


if( VC$margins[2] %in% c(cont3par) && VC$margins[3] %in% c(cont3par) ){

teta12 <- teta.tr(VC1, VC$X8s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2 + VC$X8.d2)])$teta
teta13 <- teta.tr(VC2, VC$X9s%*%SemiParFit$fit$argument[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2 + VC$X8.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2 + VC$X8.d2 + VC$X9.d2)])$teta

ass.msR1 <- ass.ms(VC$BivD1, VC$nCa1, teta12)
ass.msR2 <- ass.ms(VC$BivD2, VC$nCa2, teta13)

}
 
theta12   <- ass.msR1$theta
theta12.a <- ass.msR1$theta.a
tau12     <- ass.msR1$tau  
tau12.a   <- ass.msR1$tau.a

theta13   <- ass.msR2$theta
theta13.a <- ass.msR2$theta.a
tau13     <- ass.msR2$tau  
tau13.a   <- ass.msR2$tau.a

######################

if(VC$gc.l == TRUE) gc()  

edf.loopR <- edf.loop(VC, F, F1, GAM)
 
edf  <- edf.loopR$edf
edf1 <- edf.loopR$edf1 
  
sp <- SemiParFit$sp 
  

                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, F1 = F1, t.edf = t.edf, edf = edf, Vb.t = Vb.t,
                      edf11 = edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], edf5 = edf[[5]], edf6 = edf[[6]],
                      edf7 = edf[[7]], edf8 = edf[[8]], edf9 = edf[[9]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], edf1.5 = edf1[[5]], 
                      edf1.6 = edf1[[6]], edf1.7 = edf1[[7]], edf1.8 = edf1[[8]], edf1.9 = edf1[[9]],
                      theta12 = theta12, theta12.a = theta12.a, theta13 = theta13, theta13.a = theta13.a,
                      sigma2 = sigma2, sigma2.a = sigma2.a, sigma3 = sigma3, sigma3.a = sigma3.a,
                      nu2 = nu2, nu2.a = nu2.a, nu3 = nu3, nu3.a = nu3.a, 
                      tau12 = tau12, tau12.a= tau12.a, tau13 = tau13, tau13.a = tau13.a,
                      sp = sp, R = R, Ve = Ve, 
                      dof12.a = VC$dof1, dof12 = VC$dof1, dof13.a = VC$dof2, dof13 = VC$dof2, nCa1 = nCa1, nCa2 = nCa2) 

}



