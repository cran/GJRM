jc.probs5 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link){

######################################################################################################

CIp12 <- p12s <- nu1 <- nu1s <- nus <- nu <- NULL
CIkt <- tau <- NULL

######################################################################################################


if(type == "joint"){ 


######

if(!missing(newdata)){ 

X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")


if( !is.null(x$X3) ){

X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")
X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")

if(x$margins[1] %in% c(cont3par)) X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")
 
                    }


eta1 <- X1%*%x$coefficients[1:x$X1.d2]
eta2 <- X2%*%x$coef.t[(x$X1.d2+1):(x$X1.d2+x$X2.d2)]


if( !is.null(x$X3) ){

   sigma21 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata), x$margins[1])$vrb

   if(x$margins[1] %in% c(cont2par)) theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = 4, newdata = newdata))$teta
                                
   if(x$margins[1] %in% c(cont3par)){
     nu1   <- esp.tr(predict.SemiParBIV(x, eq = 4, newdata = newdata), x$margins[1])$vrb
     theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = 5, newdata = newdata))$teta
                                    }
                    }



if( is.null(x$X3) ){

sigma21 <- x$sigma21 
if(x$margins[1] %in% c(cont3par)) nu1 <- x$nu1 
theta <- x$theta 

}


                     }


if(missing(newdata)){

X1 <- x$X1 
X2 <- x$X2 
if( !is.null(x$X3) ){ X3 <- x$X3; X4 <- x$X4; if(x$margins[1] %in% c(cont3par)) X5 <- x$X5}

eta1  <- x$eta1
eta2  <- x$eta2
sigma21 <- x$sigma21 
if(x$margins[1] %in% c(cont3par)) nu1 <- x$nu1 
theta <- x$theta 
                     }


######

p1 <- as.numeric(distrHsAT(y1, eta1, sigma21, nu1, x$margins[1])$p2) 
p2 <- as.numeric(probmS(eta2, x$VC$margins[2])$pr)  
theta <- as.numeric(theta)
###### 





if(cond == 0){



if(!(x$BivD %in% x$BivD2)) p12 <- BiCDF(p1, p2, x$nC, theta, 3)

if(x$BivD %in% x$BivD2){

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], 3)
if(length(theta) == 1) p12[x$teta.ind1] <- BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta, 3)
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], 3)
if(length(theta) == 1) p12[x$teta.ind2] <- BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta, 3)
                         }                            
                                                                           
}


}



if(cond == 1){

if(!(x$BivD %in% x$BivD2)) p12 <- copgHsCond(p1, p2, theta, dof = 3, x$BivD)$c.copula.be1


if(x$BivD %in% x$BivD2){

p12 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta[x$teta.ind1], dof = 3, x$Cop1)$c.copula.be1
if(length(theta) == 1) p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta, dof = 3, x$Cop1)$c.copula.be1
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta[x$teta.ind2], dof = 3, x$Cop2)$c.copula.be1
if(length(theta) == 1) p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta, dof = 3, x$Cop2)$c.copula.be1
                         }                            
                                                                           
}



}



if(cond == 2){

if(!(x$BivD %in% x$BivD2)) p12 <- copgHsCond(p1, p2, theta, dof = 3, x$BivD)$c.copula.be2


if(x$BivD %in% x$BivD2){

p12 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta[x$teta.ind1], dof = 3, x$Cop1)$c.copula.be2
if(length(theta) == 1) p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta, dof = 3, x$Cop1)$c.copula.be2
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta[x$teta.ind2], dof = 3, x$Cop2)$c.copula.be2
if(length(theta) == 1) p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta, dof = 3, x$Cop2)$c.copula.be2
                         }                            
                                                                           
}


}



# kendalls' tau

if(x$BivD %in% x$BivD2)    tau <- Reg2Copost(x$SemiParFit, x$VC, theta)$tau 
if(!(x$BivD %in% x$BivD2)) tau <- ass.ms(x$VC$BivD, x$VC$nCa, theta)$tau



############################## 


if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
if(!is.null(x$VC$mono.sm.pos)) mono.sm.pos <- x$VC$mono.sm.pos else mono.sm.pos <- c(NULL, x$VC$mono.sm.pos2 + x$VC$X1.d2)  
bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos]) 

lf <- length(x$coefficients)

eta1s <- eta.tr( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1]) 


if( !is.null(x$X3) ){
sigma21s <- esp.tr( X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]), x$VC$margins[1])$vrb 
if(x$margins[1] %in% cont3par) nu1s <- esp.tr(X4%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)]), x$VC$margins[1])$vrb 
                    } 
                    
if(  is.null(x$X3) ){
   if(x$margins[1] %in% cont3par){
            sigma21s <- esp.tr(bs[, lf - 2], x$VC$margins[1])$vrb; nu1s <- esp.tr(bs[, lf - 1], x$VC$margins[1])$vrb 
            sigma21s <- matrix(rep(sigma21s, each = dim(eta1s)[1]), ncol = n.sim, byrow = FALSE)
            nu1s     <- matrix(rep(nu1s, each = dim(eta1s)[1]), ncol = n.sim, byrow = FALSE)
                                 }
   
   if(x$margins[1] %in% cont2par){
                     sigma21s <- esp.tr(bs[, lf - 1], x$VC$margins[1])$vrb 
                     sigma21s <- matrix(rep(sigma21s, each = dim(eta1s)[1]), ncol = n.sim, byrow = FALSE)
                                 }
                    }                     
   
                       
p1s <- distrHsAT(y1, eta1s, sigma21s, nu1s, x$margins[1])$p2
p2s <- probmS( X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])$pr

  
if( !is.null(x$X3) ){
  if(x$margins[1] %in% cont2par) epds <- X4%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
  if(x$margins[1] %in% cont3par) epds <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                    }
if( is.null(x$X3)  ) epds <- bs[, lf]
       
                         
est.RHOb <- teta.tr(x$VC, epds)$teta
if( is.null(x$X3) ) est.RHOb <- matrix(rep(est.RHOb, each = dim(p1s)[1]), ncol = n.sim, byrow = FALSE)

###########################   



if(cond == 0){


if(x$VC$BivD %in% c("N","T")) p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, 3, test = FALSE), dim(p1s)[1], n.sim) else{

if(!(x$BivD %in% x$BivD2)) p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, par2 = 3, test = FALSE), dim(p1s)[1], n.sim)

if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
                      
                        }

                                                                                                                        }
}

                                                                               
      
if(cond == 1){


if(!(x$BivD %in% x$BivD2)) p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = 3, x$BivD)$c.copula.be1

if(x$BivD == "T") p12s <- matrix(p12s, dim(p1s)[1], n.sim)


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])
 
if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- copgHsCond(p1s[x$teta.ind1,], p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = 3, x$Cop1)$c.copula.be1                                               
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- copgHsCond(p1s[x$teta.ind2,], p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = 3, x$Cop2)$c.copula.be1
                                                                                                     
                       }

              }




if(cond == 2){


if(!(x$BivD %in% x$BivD2)) p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = 3, x$BivD)$c.copula.be2

if(x$BivD == "T") p12s <- matrix(p12s, dim(p1s)[1], n.sim)


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])
 
if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- copgHsCond(p1s[x$teta.ind1,], p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = 3, x$Cop1)$c.copula.be2                                               
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- copgHsCond(p1s[x$teta.ind2,], p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = 3, x$Cop2)$c.copula.be2
                                                                                                     
                       }
                       
             }







nCa   <- x$VC$nCa
BivDt <- x$VC$BivD

  if(x$BivD %in% x$BivD2){
  
  if(x$BivD %in% x$BivD2[1:4]) { BivDt <- "C0"; nCa <- 3} 
  if(x$BivD %in% x$BivD2[5:8]) { BivDt <- "J0"; nCa <- 6}
  if(x$BivD %in% x$BivD2[9:12]){ BivDt <- "G0"; nCa <- 4}
  
                         }
  
ass.msR <- ass.ms(BivDt, nCa, est.RHOb)
taus    <- ass.msR$tau
if(!is.null(x$X3) && BivDt %in% c("AMH", "FGM")) taus <- matrix(taus, nrow(x$X3), nrow(bs))
CIkt <- rowQuantiles(taus, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
#if( is.null(x$X3) ) CIkt <- t(CIkt) 



 if(x$BivD %in% x$BivD2){ 
 
   if(length(x$theta) > 1){
 
     if( length(x$teta2) != 0) CIkt[x$teta.ind2, ] <- -CIkt[x$teta.ind2, ]; CIkt[x$teta.ind2, c(1,2)] <- CIkt[x$teta.ind2, c(2,1)] 
                                 
                          }else{
 
     if( length(x$teta2) != 0) CIkt <- -CIkt; CIkt[, c(1,2)] <- CIkt[, c(2,1)]
                                 
                                }
 }






} # interv


}## biv
         
         
         
         
         
######################################################################################################
######################################################################################################

if(type == "independence"){

if(!missing(newdata)){

X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")

if( !is.null(x$X3) ){

X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")
if(x$margins[1] %in% c(cont3par)) X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")
 
                    }


eta1 <- X1%*%x$gamlss1$coefficients[1:x$X1.d2]
eta2 <- X2%*%x$gamlss2$coef.t[1:x$X2.d2]

if( !is.null(x$X3) ){
   sigma2 <- esp.tr(predict.SemiParBIV(x$gamlss1, eq = 2, newdata = newdata), x$margins[1])$vrb
   if(x$margins[1] %in% c(cont3par)) nu <- esp.tr(predict.SemiParBIV(x$gamlss1, eq = 3, newdata = newdata), x$margins[1])$vrb
                    }

if( is.null(x$X3) ){
   sigma2 <- x$gamlss1$sigma2 
   if(x$margins[1] %in% cont3par) nu <- x$gamlss1$nu 

                    }
}


if(missing(newdata)){

X1 <- x$X1
X2 <- x$X2
if( !is.null(x$X3) ){ X3 <- x$X3; if(x$margins[1] %in% c(cont3par)) X4 <- x$X4}

eta1   <- x$gamlss1$eta1
eta2   <- x$gamlss2$eta1
sigma2 <- x$gamlss1$sigma2
if(x$margins[1] %in% c(cont3par)) nu <- x$gamlss1$nu 

}



p1 <- as.numeric(distrHsAT(y1, eta1, sigma2, nu, x$margins[1])$p2) 
p2 <- as.numeric(probmS(eta2, x$VC$margins[2])$pr)  
  

p12 <- p1*p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1


if(intervals == TRUE){

bs1 <- rMVN(n.sim, mean = x$gamlss1$coefficients, sigma=x$gamlss1$Vb)

bs2 <- rMVN(n.sim, mean = x$gamlss2$coef.t, sigma=x$gamlss2$Vb.t)
mono.sm.pos <- x$gamlss2$VC$mono.sm.pos 
bs2[, mono.sm.pos] <- ifelse(bs2[, mono.sm.pos] < 0, 0, bs2[, mono.sm.pos]) 




lf1 <- length(x$gamlss1$coefficients)

eta1s <- eta.tr( X1%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1]) 


if( !is.null(x$X3) ){
sigma2s <- esp.tr( X3%*%t(bs1[,(x$X1.d2+1):(x$X1.d2+x$X3.d2)]), x$VC$margins[1])$vrb 
if(x$margins[1] %in% cont3par) nus <- esp.tr(X4%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X4.d2)]), x$VC$margins[1])$vrb 
                    } 
                    
if(  is.null(x$X3) ){
   if(x$margins[1] %in% cont2par){
            sigma2s <- esp.tr(bs1[, lf1], x$VC$margins[1])$vrb 
            sigma2s <- matrix(rep(sigma2s, each = dim(eta1s)[1]), ncol = n.sim, byrow = FALSE)
                                 }
   
   if(x$margins[1] %in% cont3par){
            sigma2s <- esp.tr(bs1[, lf1-1], x$VC$margins[1])$vrb
            nus     <- esp.tr(bs1[, lf1], x$VC$margins[1])$vrb
            sigma2s <- matrix(rep(sigma2s, each = dim(eta1s)[1]), ncol = n.sim, byrow = FALSE)
            nus     <- matrix(rep(nus, each = dim(eta1s)[1]), ncol = n.sim, byrow = FALSE)

                                 }
                    }                     
   
                       
p1s <- distrHsAT(y1, eta1s, sigma2s, nus, x$margins[1])$p2
p2s <- probmS( X2%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2])$pr 

p12s <- p1s*p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p1s

                      } # intervals


} # independence


list(p12 = p12, p12s = p12s, p1 = p1, p2 = p2, p3 = NULL, CIkt = CIkt, tau = tau)


}



