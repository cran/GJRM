copulaReg.fit.post <- function(SemiParFit, VC, GAM){

Ve <- R <- theta <- edf <- edf1 <- coef.t <- NULL
theta.a <- sigma21 <- sigma22 <- sigma2.1.a <- sigma2.2.a <- nu1 <- nu2 <- nu1.a <- nu2.a <- dof <- dof.a <- nCa1 <- nCa2 <- NULL
													     
cont1par  <- c(VC$m1d,VC$bl) 
cont2par  <- c(VC$m2,VC$m2d) 
cont3par  <- c(VC$m3,VC$m3d)

if(VC$margins[1] != "LN" && VC$margins[2] != "LN") logLik <- -SemiParFit$fit$l
if(VC$margins[1] == "LN" || VC$margins[2] == "LN") logLik <- -SemiParFit$fit$l.ln


pVbres <- postVb(SemiParFit, VC)

He         <- pVbres$He        
Vb         <- pVbres$Vb   
Vb.t       <- pVbres$Vb.t
HeSh       <- pVbres$HeSh      
F          <- pVbres$F         
F1         <- pVbres$F1        
R          <- pVbres$R         
Ve         <- pVbres$Ve        
t.edf      <- pVbres$t.edf     
SemiParFit <- pVbres$SemiParFit
coef.t     <- pVbres$coef.t  


############################################################################################
# SIGMAs
############################################################################################  
# this should not bother once we set NULL
  
sigma21 <- esp.tr(SemiParFit$fit$etas1, VC$margins[1])$vrb  
sigma22 <- esp.tr(SemiParFit$fit$etas2, VC$margins[2])$vrb    
  
if( is.null(VC$X3) ) {names(sigma21) <- "sigma1"
                      names(sigma22) <- "sigma2" } 
  
  sigma2.1.a <- mean(sigma21); sigma2.2.a <- mean(sigma22)  


############################################################################################
# NUs
############################################################################################  

if(VC$margins[1] %in% cont3par ){  

nu1 <- esp.tr(SemiParFit$fit$etan1, VC$margins[1])$vrb    
if( is.null(VC$X3) ) names(nu1) <- "nu1"

nu1.a <- mean(nu1)


}


if(VC$margins[2] %in% cont3par ){  

nu2 <- esp.tr(SemiParFit$fit$etan2, VC$margins[2])$vrb    
if( is.null(VC$X3) ) names(nu2) <- "nu2"
                       
nu2.a <- mean(nu2)
 

}


############################################################################################
# THETA
############################################################################################

dep <- SemiParFit$fit$etad
if( is.null(VC$X3) ) names(dep) <- "theta"
  
theta <- teta.tr(VC, dep)$teta 


if(VC$BivD %in% VC$BivD2){

Reg2CopostR <- Reg2Copost(SemiParFit, VC, theta)

theta   <- Reg2CopostR$theta
theta.a <- Reg2CopostR$theta.a
tau     <- Reg2CopostR$tau  
tau.a   <- Reg2CopostR$tau.a

} 

if(!(VC$BivD %in% VC$BivD2)){

ass.msR <- ass.ms(VC$BivD, VC$nCa, theta)
theta <- ass.msR$theta
theta.a <- ass.msR$theta.a
tau   <- ass.msR$tau  
tau.a <- ass.msR$tau.a

}




######################

BivD <- VC$BivD

if(VC$surv == TRUE) BivD <- "N"

if(BivD == "T" && VC$margins[1] %in% c(VC$m2,VC$m3) && VC$margins[2] %in% c(VC$m2,VC$m3)){ dof <- dof.tr(SemiParFit$fit$etan)$vao  

dof.a <- mean(dof) 

if( is.null(VC$X3) ) names(dof) <- "dof" 

}else{  dof.a <- dof <- VC$dof }

######################

if(VC$gc.l == TRUE) gc()  


edf.loopR <- edf.loop(VC, F, F1, GAM)
 
edf  <- edf.loopR$edf
edf1 <- edf.loopR$edf1 

sp <- SemiParFit$sp 
   
    
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, Vb.t = Vb.t, HeSh = HeSh, Vb.t = Vb.t,
                      F = F, F1 = F1, t.edf = t.edf, edf = edf, 
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], 
                      edf5 = edf[[5]], edf6 = edf[[6]], edf7 = edf[[7]], edf8 = edf[[8]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], 
                      edf1.5 = edf1[[5]], edf1.6 = edf1[[6]], edf1.7 = edf1[[7]], edf1.8 = edf1[[8]],
                      theta = theta, theta.a = theta.a, tau = tau, tau.a= tau.a,
                      sigma21 = sigma21, sigma22 = sigma22, sigma21.a = sigma2.1.a, sigma22.a = sigma2.2.a,
                      sigma1 = sigma21, sigma2 = sigma22, sigma1.a = sigma2.1.a, sigma2.a = sigma2.2.a,
                      nu1 = nu1, nu1.a = nu1.a, nu2= nu2, nu2.a = nu2.a,
                      sp = sp, R = R, Ve = Ve, dof.a = dof.a, dof = dof, nCa1 = nCa1, nCa2 = nCa2,  coef.t = coef.t) 

}



