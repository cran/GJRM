gamlss.fit.post <- function(SemiParFit, VC, GAM){

Ve <- R <- edf <- edf1 <- coef.t <- NULL
sigma2 <- sigma2.a <- nu <- nu.a <- NULL

cont1par  <- c(VC$m1d) 
cont2par  <- c(VC$m2,VC$m2d) 
cont3par  <- c(VC$m3,VC$m3d)

if(VC$margins[1] != "LN" ) logLik <- -SemiParFit$fit$l      
if(VC$margins[1] == "LN" ) logLik <- -SemiParFit$fit$l.lnun


pVbres <- postVb(SemiParFit, VC)

He         <- pVbres$He        
Vb         <- pVbres$Vb
Vb1        <- pVbres$Vb1
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
  
if(VC$margins[1] %in% c(cont2par,cont3par) ){
  
sigma2 <- esp.tr(SemiParFit$fit$etas1, VC$margins[1])$vrb  
  
if( is.null(VC$X2) ) names(sigma2) <- "sigma"

sigma2.a <- mean(sigma2) 

}

############################################################################################
# NUs
############################################################################################  

if(VC$margins[1] %in% cont3par ){  

nu <- enu.tr(SemiParFit$fit$etan1, VC$margins[1])$vrb    
if( is.null(VC$X3) ) names(nu) <- "nu"

nu.a <- mean(nu)

}

######################

if(VC$gc.l == TRUE) gc()  

edf.loopR <- edf.loop(VC, F, F1, GAM)
 
edf  <- edf.loopR$edf
edf1 <- edf.loopR$edf1 
 
sp <- SemiParFit$sp 
   
    
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, Vb1 = Vb1, Vb.t = Vb.t, HeSh = HeSh, 
                      F = F, F1 = F1, t.edf = t.edf, edf = edf, 
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], 
                      edf5 = edf[[5]], edf6 = edf[[6]], edf7 = edf[[7]], edf8 = edf[[8]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], 
                      edf1.5 = edf1[[5]], edf1.6 = edf1[[6]], edf1.7 = edf1[[7]], edf1.8 = edf1[[8]],
                      sigma2 = sigma2, sigma2.a = sigma2.a, 
                      sigma = sigma2, sigma.a = sigma2.a,
                      nu = nu, nu.a = nu.a, sp = sp, 
                      R = R, Ve = Ve, coef.t = coef.t) 

}



