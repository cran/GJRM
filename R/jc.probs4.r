jc.probs4 <- function(x, y1, y2, newdata, type, cond, sf, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = FALSE){

######################################################################################################

CIp12 <- p12s <- NULL
CIkt <- tau <- theta <- CItheta <- NULL

p12s  <- matrix(0, 1, 2) 

param1 <- x$coef.t[1:x$X1.d2]
param2 <- x$coef.t[(x$X1.d2+1):(x$X1.d2+x$X2.d2)]

newdata[, x$v1[1]] <- ifelse(newdata[, x$v1[1]] < 0.0001, 0.0001, newdata[, x$v1[1]]) 
newdata[, x$v2[1]] <- ifelse(newdata[, x$v2[1]] < 0.0001, 0.0001, newdata[, x$v2[1]])

######################################################################################################




nC   <- x$nC
BivD <- x$BivD





if(type == "joint"){ 

######

if(!missing(newdata)){  

X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")
if( !is.null(x$X3) ) X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")


eta1 <- X1%*%param1
eta2 <- X2%*%param2

if(sf == "hazard"){
                    Xd1 <- Xdpred(x$gam1, newdata, as.character(x$formula[[1]][2]))
                    Xd2 <- Xdpred(x$gam2, newdata, as.character(x$formula[[2]][2]))
}

if( !is.null(x$X3) ) theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = 3, newdata = newdata))$teta
if( is.null(x$X3) )  theta <- x$theta 

                     }


if(missing(newdata)){

X1 <- x$X1 
X2 <- x$X2 
if( !is.null(x$X3) ) X3 <- x$X3

eta1  <- x$eta1
eta2  <- x$eta2
theta <- x$theta 

Xd1 <- x$VC$Xd1
Xd2 <- x$VC$Xd2
                     }


######

pS1 <- probmS(eta1, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)
pS2 <- probmS(eta2, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)

p1 <- as.numeric(pS1$pr) 
p2 <- as.numeric(pS2$pr)  
theta <- as.numeric(theta)
###### 





if(x$BivD %in% x$BivD2){ 

if(theta >= 0){ nC <- x$fit$nC1; BivD <- x$fit$Cop1 } else{ nC <- x$fit$nC2; BivD <- x$fit$Cop2} 

}





if(cond == 0){


p12 <- mm(BiCDF(p1, p2, nC, theta, x$dof), min.pr = min.pr, max.pr = max.pr)

#if(sf == "cumhaz") p12 <- -log(p12)
#
#if(sf == "hazard"){
#
#                    cp <- copgHs2(p1, p2, eta1 = NULL, eta2 = NULL, theta, theta, BivD, x$dof, min.pr = min.pr, max.pr = max.pr, min.dn = x$VC$min.dn)$c.copula2.be1be2
#
#                    Gp1 <- pS1$dS
#                    Gp2 <- pS2$dS
#                    
#                    Xt1 <- Xd1 %*% param1
#                    Xt2 <- Xd2 %*% param2
#                                       
#                    p12 <- (cp*Gp1*Gp2*Xt1*Xt2)/p12
#
#                                      
#                  }
#


}







if(cond == 1){




p12 <- copgHsCond(p1, p2, theta, dof = x$dof, BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be1


#if(sf == "cumhaz") p12 <- -log(p12)
#
#if(sf == "hazard"){
#
#                    cp <- copgHs2(p1, p2, eta1 = NULL, eta2 = NULL, theta, theta, BivD, x$dof, min.pr = min.pr, max.pr = max.pr, min.dn = x$VC$min.dn)$c.copula2.be1
#
#                    Gp1 <- pS1$dS                  
#                    Xt1 <- Xd1 %*% param1
#                    
#                    
#                    p12 <- (cp*Gp1*Xt1) / p12
#
#
#}




}



if(cond == 2){



p12 <- copgHsCond(p1, p2, theta, dof = x$dof, BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2


#if(sf == "cumhaz") p12 <- -log(p12)
#
#if(sf == "hazard"){
#
#                    cp <- copgHs3(p1, p2, eta1 = NULL, eta2 = NULL, theta, theta, BivD, x$dof)$c.copula2.be2
#
#                    Gp2 <- pS2$dS                  
#                    Xt2 <- Xd2 %*% param2
#                    
#                    p12 <- (cp*Gp2*Xt2)/p12
#
#
#}


}




# kendalls' tau

 tau   <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = tau.res, t.theta = theta)$tau
 theta <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = FALSE,   t.theta = theta)$theta
     

############################## 


if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
if(!is.null(x$VC$mono.sm.pos)) mono.sm.pos <- x$VC$mono.sm.pos else mono.sm.pos <- c(x$VC$mono.sm.pos1, x$VC$mono.sm.pos2 + x$VC$X1.d2)  
bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos]) 

lf <- length(x$coefficients)

param1s <- t(bs[,1:x$X1.d2])
param2s <- t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)])  
  
P1s <- probmS( X1%*%param1s, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr) 
P2s <- probmS( X2%*%param2s, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)
     
p1s <- P1s$pr 
p2s <- P2s$pr

  
if( !is.null(x$X3) ) epds <- X3%*%t(bs[,(x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
if( is.null(x$X3)  ) epds <- bs[, lf]
       
                         
est.RHOb <- teta.tr(x$VC, epds)$teta




#if(sf == "hazard"){
#
#Gp1s <- P1s$dS
#Gp2s <- P2s$dS
#
#Xt1s <- Xd1 %*% param1s
#Xt2s <- Xd2 %*% param2s
#
#                                     
#}



###########################   



if(cond == 0){


p12s <- mm(BiCDF(p1s, p2s, nC, est.RHOb, x$dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)

#if(sf == "cumhaz") p12s <- -log(p12s) 
#
#if(sf == "hazard"){
#
#            cps <- copgHs2(p1s, p2s, eta1 = NULL, eta2 = NULL, est.RHOb, est.RHOb, BivD, x$dof, min.pr = min.pr, max.pr = max.pr, min.dn = x$VC$min.dn)$c.copula2.be1be2
#
#            p12s <- (cps*Gp1s*Gp2s*Xt1s*Xt2s)/p12s
#
#                  }
#



                   
}         
         
         
         
         
         
if(cond == 1){


p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = x$dof, BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be1


#if(sf == "cumhaz") p12s <- -log(p12s) 
#
#
#if(sf == "hazard"){
#
#                    cps <- copgHs2(p1s, p2s, eta1 = NULL, eta2 = NULL, est.RHOb, est.RHOb, BivD, x$dof, min.pr = min.pr, max.pr = max.pr, min.dn = x$VC$min.dn)$c.copula2.be1
#
#                    p12s <- (cps*Gp1s*Xt1s) / p12s
#
#
#}
#



}




if(cond == 2){


p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = x$dof, BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2


#if(sf == "cumhaz") p12s <- -log(p12s)
#
#if(sf == "hazard"){
#
#                    cps <- copgHs3(p1s, p2s, eta1 = NULL, eta2 = NULL, est.RHOb, est.RHOb, BivD, x$dof)$c.copula2.be2
#                    
#                    p12s <- (cps*Gp2s*Xt2s)/p12s
#
#
#}




}






  
ass.msR <- theta2tau(x$VC$BivD, x$VC$nCa, est.RHOb, tau.res = tau.res, t.theta = theta)
if(tau.res == TRUE) taus <- ass.msR$tau
thetas <- ass.msR$theta




                                                  
if(tau.res == TRUE){taus <- matrix(taus, 1, n.sim); CIkt <- rowQuantiles(taus, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)}

thetas <- matrix(thetas, 1, n.sim)
CItheta <- rowQuantiles(thetas, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)





} # interv


}## biv
                        
######################################################################################################
######################################################################################################

if(type == "independence"){

param1 <- x$gamlss1$coef.t[1:x$X1.d2]
param2 <- x$gamlss2$coef.t[1:x$X2.d2]


if(!missing(newdata)){

X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")

eta1 <- X1%*%param1
eta2 <- X2%*%param2


#if(sf == "hazard"){
#                    Xd1 <- Xdpred(x$gam1, newdata, as.character(x$formula[[1]][2]))
#                    Xd2 <- Xdpred(x$gam2, newdata, as.character(x$formula[[2]][2]))
#                   }
#

}

if(missing(newdata)){

X1 <- x$X1
X2 <- x$X2

eta1 <- x$gamlss1$eta1
eta2 <- x$gamlss2$eta1

Xd1 <- x$VC$Xd1
Xd2 <- x$VC$Xd2

}


pS1 <- probmS(eta1, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)
pS2 <- probmS(eta2, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)

p1  <- as.numeric(pS1$pr) 
p2  <- as.numeric(pS2$pr) 




if(cond == 0){

p12 <- mm(p1*p2, min.pr = min.pr, max.pr = max.pr)

#if(sf == "cumhaz") p12 <- -log(p12)
#
#if(sf == "hazard"){
#                    Gp1 <- pS1$dS
#                    Gp2 <- pS2$dS
#                    
#                    Xt1 <- Xd1 %*% param1
#                    Xt2 <- Xd2 %*% param2
#                                       
#                    p12 <- (Gp1*Gp2*Xt1*Xt2)/p12
#                   }
#

}


if(cond == 1){

               p12 <- p2

#if(sf == "cumhaz") p12 <- -log(p12)
#
#
#if(sf == "hazard"){
#                    Gp2 <- pS2$dS                  
#                    Xt2 <- Xd2 %*% param2                   
#                    p12 <- -(Gp2*Xt2)/p12
#                   }
#


}





if(cond == 2){

               p12 <- p1

#if(sf == "cumhaz") p12 <- -log(p12)
#
#
#if(sf == "hazard"){
#                    Gp1 <- pS1$dS                  
#                    Xt1 <- Xd1 %*% param1                   
#                    p12 <- -(Gp1*Xt1)/p12
#                   }
#


}






if(intervals == TRUE){

bs1 <- rMVN(n.sim, mean = x$gamlss1$coef.t, sigma=x$gamlss1$Vb.t)
mono.sm.pos <- x$gamlss1$VC$mono.sm.pos 
bs1[, mono.sm.pos] <- ifelse(bs1[, mono.sm.pos] < 0, 0, bs1[, mono.sm.pos]) 

bs2 <- rMVN(n.sim, mean = x$gamlss2$coef.t, sigma=x$gamlss2$Vb.t)
mono.sm.pos <- x$gamlss2$VC$mono.sm.pos 
bs2[, mono.sm.pos] <- ifelse(bs2[, mono.sm.pos] < 0, 0, bs2[, mono.sm.pos]) 

#############  
# etas
############# 

param1s <- t(bs1[,1:x$X1.d2])
param2s <- t(bs2[,1:x$X2.d2])

pS1s <- probmS( X1%*%param1s, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)
pS2s <- probmS( X2%*%param2s, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr) 

p1s <- pS1s$pr
p2s <- pS2s$pr



if(cond == 0){

p12s <- mm(p1s*p2s, min.pr = min.pr, max.pr = max.pr)

#if(sf == "cumhaz") p12s <- -log(p12s)
#
#
#if(sf == "hazard"){
#                    Gp1s <- pS1s$dS
#                    Gp2s <- pS2s$dS
#                    
#                    Xt1s <- Xd1 %*% param1s
#                    Xt2s <- Xd2 %*% param2s
#                                       
#                    p12s <- (Gp1s*Gp2s*Xt1s*Xt2s)/p12s
#                   }
#
#



}



if(cond == 1){ p12s <- p2s

#if(sf == "cumhaz") p12s <- -log(p12s)
#
#
#if(sf == "hazard"){
#                    Gp2s <- pS2s$dS
#                    
#                    Xt2s <- Xd2 %*% param2s
#                                       
#                    p12s <- -(Gp2s*Xt2s)/p12s
#                   }
#


}




if(cond == 2){ p12s <- p1s

#if(sf == "cumhaz") p12s <- -log(p12s)
#
#if(sf == "hazard"){
#                    Gp1s <- pS1s$dS
#                    
#                    Xt1s <- Xd1 %*% param1s
#                                       
#                    p12s <- -(Gp1s*Xt1s)/p12s
#                    
#                   }

}





} # intervals




} # independence


list(p12 = p12, p12s = matrix(p12s, 1, length(p12s)), p1 = p1, p2 = p2, p3 = NULL, CItau = CIkt, tau = tau, CItheta = CItheta, theta = theta)


}



