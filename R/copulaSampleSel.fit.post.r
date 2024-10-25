copulaSampleSel.fit.post <- function(SemiParFit, VC, GAM){

Ve <- R <- X2s <- eta1S <- eta2S <- theta <- edf <- edf1 <- theta.a <-  sigma <- sigma.a <- sigma2 <- sigma2.a <- p1n <- p2n <- nu <- nu.a <- NULL

cont1par  <- VC$m1d
cont2par  <- c(VC$m2,VC$m2d) 
cont3par  <- VC$m3 
bin.link  <- VC$bl  

if(VC$margins[2] != "LN") logLik <- -SemiParFit$fit$l else logLik <- -SemiParFit$fit$l.ln 

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

############################################
# complete Matrices
############################################

SemiParFit$fit$eta2 <- VC$X2s%*%SemiParFit$fit$argument[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]  


if(is.null(VC$X3)){ # START

if(!(VC$margins[2] %in% cont1par)){ ##

sigma <- sigma2 <- sigma2.a <- sigma.a <- esp.tr(SemiParFit$fit$etas, VC$margins[2])$vrb 
names(sigma) <- names(sigma.a) <- names(sigma2) <- names(sigma2.a) <- "sigma"   


if(VC$margins[2] %in% cont3par ){

	if(VC$margins[2] %in% c("DAGUM","SM","TW")){

		nu <- nu.a <- enu.tr(SemiParFit$fit$etan, VC$margins[2])$vrb
		names(nu) <- names(nu.a) <- "nu"
                                              }  
 
                                }

} ##


dep        <- SemiParFit$fit$etad
names(dep) <- "theta" 

theta <- teta.tr(VC, dep)$teta   
   
} # FINISH  
  
############################################
############################################  
  
  
if(!is.null(VC$X3)){ # START  
  
  
if(!(VC$margins[2] %in% cont1par)){##  
  
  SemiParFit$fit$etas <- VC$X3s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]  
  
  sigma <- sigma2 <- esp.tr(SemiParFit$fit$etas, VC$margins[2])$vrb 
  sigma.a <- sigma2.a <- mean(sigma2)  
 
if(VC$margins[2] %in% cont2par){

SemiParFit$fit$etad <- VC$X4s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] 
theta <- teta.tr(VC, SemiParFit$fit$etad)$teta  

}


if(VC$margins[2] %in% cont3par){

SemiParFit$fit$etan <- VC$X4s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
SemiParFit$fit$etad <- VC$X5s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]
 
  nu    <- enu.tr(SemiParFit$fit$etan, VC$margins[2])$vrb   
  theta <- teta.tr(VC, SemiParFit$fit$etad)$teta  
  nu.a <- mean(nu) 
 
}

}##

if(VC$margins[2] %in% cont1par){

SemiParFit$fit$etad <- VC$X3s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
theta <- teta.tr(VC, SemiParFit$fit$etad)$teta  

}


}


######################
# Association measures
######################

ass.msR <- theta2tau(VC$BivD, VC$nCa, theta)
theta <- ass.msR$theta
theta.a <- ass.msR$theta.a
#tau   <- ass.msR$tau  
#tau.a <- ass.msR$tau.a


#############################################################

if(VC$gc.l == TRUE) gc()  

edf.loopR <- edf.loop(VC, F, F1, GAM)
 
edf  <- edf.loopR$edf
edf1 <- edf.loopR$edf1 

sp <- SemiParFit$sp 
  

###############################################################################

if(!is.null(theta.a))  names(theta.a)  <- "theta.a"
if(!is.null(sigma2.a)) names(sigma2.a) <- "sigma2.a"
if(!is.null(sigma.a))  names(sigma.a)  <- "sigma.a"
if(!is.null(nu.a))     names(nu.a)     <- "nu.a"

if(!is.null(theta)  && length(theta)  == 1) names(theta)  <- "theta"
if(!is.null(sigma2) && length(sigma2) == 1) names(sigma2) <- "sigma2"
if(!is.null(sigma)  && length(sigma)  == 1) names(sigma)  <- "sigma"
if(!is.null(nu)     && length(nu)     == 1) names(nu)     <- "nu"


if(!is.null(theta)  && length(theta)  > 1){ theta  <- as.matrix(theta  ); dimnames(theta )[[1]] <- dimnames(VC$X1)[[1]]; dimnames(theta )[[2]] <- "theta" }
if(!is.null(sigma2) && length(sigma2) > 1){ sigma2 <- as.matrix(sigma2) ; dimnames(sigma2)[[1]] <- dimnames(VC$X1)[[1]]; dimnames(sigma2)[[2]] <- "sigma2"}
if(!is.null(sigma)  && length(sigma)  > 1){ sigma  <- as.matrix(sigma)  ; dimnames(sigma )[[1]] <- dimnames(VC$X1)[[1]]; dimnames(sigma )[[2]] <- "sigma" }
if(!is.null(nu)     && length(nu)     > 1){ nu     <- as.matrix(nu    ) ; dimnames(nu    )[[1]] <- dimnames(VC$X1)[[1]]; dimnames(nu    )[[2]] <- "nu"    }



  
  
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, F1 = F1, t.edf = t.edf, edf = edf, Vb.t = Vb.t,
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], edf5 = edf[[5]], edf6 = edf[[6]],
                      edf7 = edf[[7]], edf8 = edf[[8]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], edf1.5 = edf1[[5]], 
                      edf1.6 = edf1[[6]], edf1.7 = edf1[[7]], edf1.8 = edf1[[8]],
                      theta = theta, theta.a = theta.a, sigma2 = sigma2, sigma2.a = sigma2.a, sigma = sigma, sigma.a = sigma.a,
                      nu = nu, nu.a = nu.a, #tau = tau, tau.a = tau.a,
                      sp = sp,  
                      p1n=p1n, p2n=p2n, R = R, Ve = Ve, dof.a = VC$dof, dof = VC$dof) 

}



