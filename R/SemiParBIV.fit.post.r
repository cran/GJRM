SemiParBIV.fit.post <- function(SemiParFit, Model, VC, GAM){

Ve <- R <- X2s <- X3s <- eta1S <- eta2S <- theta <- edf <- edf1 <- theta.a <- sigma <- sigma.a <- sigma2 <- sigma2.a <- OR <- GM <- p1n <- p2n <- nu <- nu.a <- NULL
dof <- dof.a <- nCa1 <- nCa2 <- NULL


cont1par  <- VC$m1d
cont2par  <- c(VC$m2,VC$m2d) 
cont3par  <- VC$m3 
bin.link  <- VC$bl  

if(VC$margins[2] != "LN") logLik <- -SemiParFit$fit$l else logLik <- -SemiParFit$fit$l.ln 

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
  
  
if(VC$hess == FALSE) SemiParFit$fit$Fisher <- SemiParFit$fit$hessian


##################################################################

if(VC$margins[2] %in% cont2par ){

sigma <- sigma2 <- esp.tr(SemiParFit$fit$etas, VC$margins[2])$vrb 
#if( is.null(VC$X3) ) names(sigma) <- names(sigma2) <- "sigma"

sigma.a <- sigma2.a <- mean(sigma2) 

}

#################################################################

if(VC$margins[2] %in% cont3par ){

#if(VC$margins[2] %in% c("DAGUM","SM","TW")){

sigma <- sigma2 <- esp.tr(SemiParFit$fit$etas, VC$margins[2])$vrb 
nu     <- enu.tr(SemiParFit$fit$etan, VC$margins[2])$vrb 

#if( is.null(VC$X4) && is.null(VC$X5) ) {names(sigma) <- names(sigma2) <- "sigma"; names(nu) <- "nu"}
           
#                                      }                            
  sigma.a <- sigma2.a <- mean(sigma2)
  nu.a     <- mean(nu)

}

############################################################################################
# THETA
############################################################################################

if(is.null(VC$theta.fx)){


if(VC$Model == "BPO0" ) dep <- theta <- theta.a <- 0

if(!(VC$Model %in% c("BPO0","BSS")) ){

  theta <- teta.tr(VC, SemiParFit$fit$etad)$teta
  #if( is.null(VC$X3) ) names(theta) <- "theta" 
                                     }

}


if(!is.null(VC$theta.fx)) theta <- VC$theta.fx 

############################################################################################



   
if(Model=="BSS"){

  SemiParFit$fit$eta2 <- VC$X2s%*%SemiParFit$fit$argument[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  
  p1n <- predict.gam(GAM$gam1, type="response")
  p2n <- probm(VC$X2s%*%GAM$gam2$coefficients, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr
 

if(!is.null(VC$X3)) {

SemiParFit$fit$etad <- VC$X3s%*%SemiParFit$fit$argument[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
theta <- teta.tr(VC, SemiParFit$fit$etad)$teta 


}




if(is.null(VC$X3)){

  theta <- teta.tr(VC, SemiParFit$fit$etad)$teta
  #names(theta) <- "theta" 

                   }  
                  
                  
                  
                  
}




############################################################################################




if(Model=="BSS" || Model=="BPO" || Model=="BPO0"){

  p1 <- probm(SemiParFit$fit$eta1, VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr 
  p2 <- probm(SemiParFit$fit$eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr # eta2 modified above
  
   p11 <- mm(BiCDF(p1, p2, VC$nC, theta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  
   SemiParFit$fit$p10 <- p1 - p11
   SemiParFit$fit$p11 <- p11
   SemiParFit$fit$p00 <- (1 - p2) - ( p1 - p11 )
   SemiParFit$fit$p01 <- p2 - p11
   
   SemiParFit$fit$p1 <- p1
   SemiParFit$fit$p2 <- p2

}

######################
# Association measures
######################


if(VC$BivD %in% VC$BivD2){

Reg2CopostR <- Reg2Copost(SemiParFit, VC, theta)

theta   <- Reg2CopostR$theta
theta.a <- Reg2CopostR$theta.a
#tau     <- Reg2CopostR$tau  
#tau.a   <- Reg2CopostR$tau.a

} 

if(!(VC$BivD %in% VC$BivD2)){

ass.msR <- theta2tau(VC$BivD, VC$nCa, theta)
theta <- ass.msR$theta
theta.a <- ass.msR$theta.a
#tau   <- ass.msR$tau  
#tau.a <- ass.msR$tau.a

}


######################

if(VC$margins[2] %in% bin.link && !(VC$Model %in% c("BPO0", "BSS"))){

p00 <- SemiParFit$fit$p00 
p01 <- SemiParFit$fit$p01 
p11 <- SemiParFit$fit$p11 
p10 <- SemiParFit$fit$p10 

p1 <- SemiParFit$fit$p1
p2 <- SemiParFit$fit$p2

OR <- (p00*p11)/(p01*p10)

OR  <- ifelse(OR  ==  Inf,  8.218407e+307, OR ) 
OR  <- ifelse(OR  == -Inf, -8.218407e+307, OR ) 

GM <- mean((OR - 1)/(OR + 1))
OR <- mean(OR)


rm(p00,p01,p10,p11,p1,p2)

}

if(VC$gc.l == TRUE) gc()  


edf.loopR <- edf.loop(VC, F, F1, GAM)
 
edf  <- edf.loopR$edf
edf1 <- edf.loopR$edf1 
  
sp <- SemiParFit$sp 
  
  
if(!is.null(theta.a))  names(theta.a)  <- "theta.a"
if(!is.null(nu.a))     names(nu.a)     <- "nu.a"
if(!is.null(sigma2.a)) names(sigma2.a) <- "sigma2.a"
if(!is.null(sigma.a))  names(sigma.a)  <- "sigma.a"

if(!is.null(theta) && length(theta)  == 1) names(theta)  <- "theta"
if(!is.null(nu)    && length(nu)     == 1) names(nu)     <- "nu"
if(!is.null(sigma2)&& length(sigma2) == 1) names(sigma2) <- "sigma2"
if(!is.null(sigma) && length(sigma)  == 1) names(sigma)  <- "sigma"

if(!is.null(OR) ) names(OR)  <- "OR"
if(!is.null(GM) ) names(GM)  <- "GM"


if(!is.null(theta) && length(theta)  > 1){ theta  <- as.matrix(theta);  dimnames(theta)[[1]]  <- dimnames(VC$X1)[[1]]; dimnames(theta)[[2]]  <- "theta" }
if(!is.null(nu)    && length(nu)     > 1){ nu     <- as.matrix(nu);     dimnames(nu)[[1]]     <- dimnames(VC$X1)[[1]]; dimnames(nu)[[2]]     <- "nu"    }
if(!is.null(sigma2)&& length(sigma2) > 1){ sigma2 <- as.matrix(sigma2); dimnames(sigma2)[[1]] <- dimnames(VC$X1)[[1]]; dimnames(sigma2)[[2]] <- "sigma2"}
if(!is.null(sigma) && length(sigma)  > 1){ sigma  <- as.matrix(sigma);  dimnames(sigma)[[1]]  <- dimnames(VC$X1)[[1]]; dimnames(sigma)[[2]]  <- "sigma" }


                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, F1 = F1, t.edf = t.edf, edf = edf, Vb.t = Vb.t,
                      edf11=edf1,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], edf4 = edf[[4]], edf5 = edf[[5]], edf6 = edf[[6]],
                      edf7 = edf[[7]], edf8 = edf[[8]],
                      edf1.1 = edf1[[1]], edf1.2 = edf1[[2]], edf1.3 = edf1[[3]], edf1.4 = edf1[[4]], edf1.5 = edf1[[5]], 
                      edf1.6 = edf1[[6]], edf1.7 = edf1[[7]], edf1.8 = edf1[[8]],
                      theta = theta, theta.a = theta.a, sigma2 = sigma2, sigma2.a = sigma2.a, sigma = sigma, sigma.a = sigma.a,
                      nu = nu, nu.a = nu.a, #tau = tau, tau.a= tau.a,
                      sp = sp, OR = OR, GM = GM, p1n = p1n, p2n = p2n, R = R, Ve = Ve, 
                      dof.a = VC$dof, dof = VC$dof, nCa1 = nCa1, nCa2 = nCa2) 

}



