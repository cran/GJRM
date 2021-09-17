SemiParROY.fit.post <- function(SemiParFit, Model, VC, GAM){

Ve <- R <- X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- eta1S <- eta2S <- theta1 <- theta2 <- edf <- edf1 <- theta1.a <- theta2.a <- NULL
sigma1 <- sigma1.a <- sigma2 <- sigma2.a <- nu1 <- nu1.a <- nu2 <- nu2.a <- NULL
dof1 <- dof1.a <- dof2 <- dof2.a <- nCa1 <- nCa2 <- NULL

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
  
  
############################################################################################  

if(VC$margins[2] %in% c(cont2par, cont3par) ){

sigma1   <- esp.tr(SemiParFit$fit$sigma1.st, VC$margins[2])$vrb  
sigma1.a <- mean(sigma1)
  
} 

if(VC$margins[3] %in% c(cont2par, cont3par) ){

sigma2   <- esp.tr(SemiParFit$fit$sigma2.st, VC$margins[3])$vrb  
sigma2.a <- mean(sigma2)
  
}  
  

#################################################################

if(VC$margins[2] %in% cont3par ){  

nu1   <- enu.tr(SemiParFit$fit$nu1.st, VC$margins[2])$vrb    
nu1.a <- mean(nu1)

}

if(VC$margins[3] %in% cont3par ){  

nu2   <- enu.tr(SemiParFit$fit$nu2.st, VC$margins[3])$vrb    
nu2.a <- mean(nu2)

}

############################################################################################

ass.msR1 <- ass.ms(VC$BivD1, VC$nCa1, SemiParFit$fit$teta1)
ass.msR2 <- ass.ms(VC$BivD2, VC$nCa2, SemiParFit$fit$teta2)
  
theta1   <- ass.msR1$theta
theta1.a <- ass.msR1$theta.a
tau1     <- ass.msR1$tau  
tau1.a   <- ass.msR1$tau.a

theta2   <- ass.msR2$theta
theta2.a <- ass.msR2$theta.a
tau2     <- ass.msR2$tau  
tau2.a   <- ass.msR2$tau.a

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
                      theta1 = theta1, theta1.a = theta1.a, theta2 = theta2, theta2.a = theta2.a,
                      sigma2 = sigma2, sigma2.a = sigma2.a, sigma1 = sigma1, sigma1.a = sigma1.a,
                      nu1 = nu1, nu1.a = nu1.a, nu2 = nu2, nu2.a = nu2.a, 
                      tau1 = tau1, tau1.a= tau1.a, tau2 = tau2, tau2.a= tau2.a,
                      sp = sp, R = R, Ve = Ve, 
                      dof1.a = VC$dof1, dof1 = VC$dof1, dof2.a = VC$dof2, dof2 = VC$dof2, nCa1 = nCa1, nCa2 = nCa2) 

}



