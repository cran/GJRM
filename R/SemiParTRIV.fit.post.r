SemiParTRIV.fit.post <- function(SemiParFit, VC, Model, GAM){

Ve <- R <- X2s <- X3s <- eta1S <- eta2S <- theta <- edf <- edf1 <- theta.a <- p1n <- p2n <- p3n <- NULL

logLik <- -SemiParFit$fit$l

pVbres <- postVb(SemiParFit, VC)

He         <- pVbres$He        
Vb         <- pVbres$Vb        
HeSh       <- pVbres$HeSh   
Vb.t       <- pVbres$Vb.t

F          <- pVbres$F         
F1         <- pVbres$F1        
R          <- pVbres$R         
Ve         <- pVbres$Ve        
t.edf      <- pVbres$t.edf     
SemiParFit <- pVbres$SemiParFit

if(VC$hess == FALSE) SemiParFit$fit$Fisher <- SemiParFit$fit$hessian

############################################################################################
# THETAs
############################################################################################

theta12 <- SemiParFit$fit$theta12     
theta13 <- SemiParFit$fit$theta13     
theta23 <- SemiParFit$fit$theta23  

if(is.null(VC$X4)){

names(theta12) <- "theta12"
names(theta13) <- "theta13" 
names(theta23) <- "theta23" 

}

theta12.a  <- mean(theta12) 
theta13.a  <- mean(theta13)
theta23.a  <- mean(theta23)   

############################################################################################


  if(Model %in% c("TSS","TESS")){

  SemiParFit$fit$eta2 <- VC$X2s%*%SemiParFit$fit$argument[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  SemiParFit$fit$eta3 <- VC$X3s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  
  p1n <- predict.gam(GAM$gam1, type = "response")
  p2n <- probm(VC$X2s%*%GAM$gam2$coefficients, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr 
  p3n <- probm(VC$X3s%*%GAM$gam3$coefficients, VC$margins[3], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr 
 
}



if(VC$gc.l == TRUE) gc()  


if( !(VC$penCor %in% c("unpen") && VC$l.flist == 6) ) VC$l.sp4 <- 0  
# in the previous version we had VC$l.sp4 <- 0 but with 6 eqs and unpen corrs we need l.sp4 and can't set it to 0 

edf.loopR <- edf.loop(VC, F, F1, GAM)
 
edf  <- edf.loopR$edf
edf1 <- edf.loopR$edf1 
  
sp <- SemiParFit$sp 
  
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, Vb.t = Vb.t,
                      HeSh = HeSh, F = F, F1 = F1, t.edf = t.edf, edf = edf, 
                      edf11 = edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], 
                      edf5 = edf[[5]], edf6 = edf[[6]], edf7 = edf[[7]], edf8 = edf[[8]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], 
                      edf1.5 = edf1[[5]], edf1.6 = edf1[[6]], edf1.7 = edf1[[7]], edf1.8 = edf1[[8]],
                      theta12 = theta12, theta12.a = theta12.a, 
                      theta13 = theta13, theta13.a = theta13.a,
                      theta23 = theta23, theta23.a = theta23.a,
                      sp = sp, R = R, Ve = Ve,
                      p1n = p1n, p2n = p2n, p3n = p3n) 

}



