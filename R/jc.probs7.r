jc.probs7 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = FALSE){

if(!(y1 >= range(x$y1)[1] && y1 <= range(x$y1)[2])) stop("y1 is not within its observed range.") 
if(!(y2 >= range(x$y2)[1] && y2 <= range(x$y2)[2])) stop("y2 is not within its observed range.") 

if(x$BivD %in% x$BivD2) stop("The chosen copula is not supported by this function.")


# Ordinal-Continuous case

nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- CIkt <- tau <- theta <- CItheta <- NULL
dof <- x$dof

pk1 <- p12.s <- p12s.s <- pks1 <- 0

infty <- 1e+25

p12s  <- matrix(0, 1, 2) 


if(type == "joint"){


pkt <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum 
pk  <- pkt[, y1]
if(y1 > 1) pk1 <- pkt[, y1 - 1]

p1m <- pk - pk1                         
                                                    

eta2 <- predict(x, eq = 2, newdata = newdata, type = "link")


if( !is.null(x$X3) && x$margins[2] %in% c(cont2par, cont3par) ){

sigma2 <- esp.tr(predict(x, eq = 3, newdata = newdata), x$margins[2])$vrb

if(x$margins[2] %in% cont2par)  eq.th <- 4
if(x$margins[2] %in% cont3par){ eq.nu <- 4; eq.th <- 5} # case of 3 par margin not implemented yet

if(x$margins[2] %in% cont3par)  nu <- enu.tr(predict(x, eq = eq.nu, newdata = newdata), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict(x, eq = eq.th, newdata = newdata))$teta

                                                                }

if( is.null(x$X3) ){

   sigma2 <- x$sigma2
   nu     <- x$nu 
   theta  <- x$theta 

                   }






if(cond == 0 || cond == 1 ){ 

p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2

if(y1 > 1){
            p12.f <- mm(BiCDF(pk,  p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr)
            p12.s <- mm(BiCDF(pk1, p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr)
          }
          
if(y1 == 1) p12.f <- mm(BiCDF(pk,  p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr)
         
p12 <- p12.f - p12.s          

if(cond == 1) p12 <- p12/p1m

                     
                            }  




if(cond == 2){ 

   
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2

if(y1 > 1){
            p12.f <- copgHsCond(pk,  p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
            p12.s <- copgHsCond(pk1, p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
            
          }

if(y1 == 1) p12.f <- copgHsCond(pk,  p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2


p12 <- p12.f - p12.s

                   
             } 





#######################################################################################
### kendalls' tau and theta ###

if(tau.res == TRUE) tau   <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = tau.res)$tau
                    theta <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = FALSE)$theta
     



#############################
##### intervals == TRUE #####
#############################

if(intervals == TRUE){

# Cut points need to be transformed first...

cut.sim <- x$coefficients[1 : (x$VC$K1 - 1)]
	cut.sim.ti <- rep(0, x$VC$K1 - 1)
	cut.sim.ti[1] <- cut.sim[1] ; for(i in 2 : (x$VC$K1 - 1)) {cut.sim.ti[i] <- sqrt(cut.sim[i] - cut.sim[i - 1])}

coefficients.s <- x$coefficients
	coefficients.s[1 : (x$VC$K1 - 1)] <- cut.sim.ti


bs <- rMVN(n.sim, mean = coefficients.s, sigma = x$Vb)  


# ... and then transformed back

cut.sim.ti <- bs[, 1 : (x$VC$K1 - 1)]
	cut.sim <- matrix(nrow = n.sim, ncol = x$VC$K1 - 1, 0)	
	cut.sim[, 1] <- cut.sim.ti[, 1] ; for (i in 2 : (x$VC$K1 - 1)) cut.sim[, i] <- cut.sim[, i - 1] + cut.sim.ti[, i]^2

bs[, 1 : (x$VC$K1 - 1)] <- cut.sim


#############  
# etas
############# 


X1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix") 

              
# Adjustments for the ordinal-continuous model

if (!is.null(x$VC$K1)) {
	CLM.shift  <- x$VC$K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
                       } else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
                               }


eta1s <- X1 %*% t(bs[, (CLM.shift2 + 1) : (x$X1.d2 + CLM.shift2)])
n.s   <- 1


lp1s <- lp1s1 <- matrix(nrow = 1, ncol = n.sim, 0)

for (i in 1 : n.sim) {

	c1s.m   <- t(matrix(nrow = x$VC$K1 - 1, ncol = n.s        , bs   [i, 1 : CLM.shift2])) 
	eta1s.m <-   matrix(nrow = n.s        , ncol = x$VC$K1 - 1, eta1s[ , i]              )

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	           lp1s[, i]  <- lp1s_i[, y1] 
	if(y1 > 1) lp1s1[, i] <- lp1s_i[, y1 - 1] 

                      }

###


           pks  <- probm(lp1s,  x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
if(y1 > 1) pks1 <- probm(lp1s1, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

eta2s <- eta.tr( X2 %*% t(bs[, (x$X1.d2 + CLM.shift2 + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift2)]), x$VC$margins[2] )


#############  
# thetas
#############  

if( is.null(x$X3) ) epds <- bs[, length(x$coefficients)]
  
if(!is.null(x$X3) ){ 

  if(x$VC$margins[2] %in% cont2par){   
  
                X4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")
                epds <- X4%*%t(bs[,((x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)) + CLM.shift2])
  
                                   }
         
         
  if(x$VC$margins[2] %in% cont3par){
  
                X5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")  
                epds <- X5%*%t(bs[,((x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)) + CLM.shift2])
                
                                   }
                   }



                    est.RHOb <- teta.tr(x$VC, epds)$teta   


                    ass.msR <- theta2tau(x$VC$BivD, x$VC$nCa, est.RHOb, tau.res = tau.res)
if(tau.res == TRUE) taus    <- ass.msR$tau
                    thetas  <- ass.msR$theta




                                                 
if(tau.res == TRUE){ taus <- matrix(taus, 1, n.sim)
                     CIkt <- rowQuantiles(taus,   probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
                    }
                    thetas  <- matrix(thetas, 1, n.sim) 
                    CItheta <- rowQuantiles(thetas, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   

#############  
# sigmas
#############  

if( x$VC$margins[2] %in% cont2par ){  

      if( is.null(x$X3) )  sigma2.star <- bs[, CLM.shift2 + x$X1.d2 + x$X2.d2 + 1] 
      if(!is.null(x$X3) ){
      
                       X3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                       sigma2.star <- X3 %*% t(bs[,((x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)) + CLM.shift2]) 
  
                          }
  
       sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
                                    }   


##### cond == 0 and cond == 1 #####

if(cond == 0 || cond == 1){  
       
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2

if(y1 > 1){ 
            p12s.f <- mm(BiCDF(pks,  p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr )
            p12s.s <- mm(BiCDF(pks1, p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr )

           }
           
if(y1 == 1) p12s.f <- mm(BiCDF(pks,  p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr )
          
          
p12s <- p12s.f - p12s.s          


if(cond == 1) p12s <- p12s/(pks - pks1)


                           }   


 
if(cond == 2){
     
     
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2


if(y1 > 1){
            p12s.f <- copgHsCond(pks,  p2s, est.RHOb, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
            p12s.s <- copgHsCond(pks1, p2s, est.RHOb, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
            
          }

if(y1 == 1) p12s.f <- copgHsCond(pks,  p2s, est.RHOb, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2


p12s <- p12s.f - p12s.s

              }



} # int


  
} # biv





###########################################
###########################################
###########################################


if(type == "independence"){


#coefficients.ind <- x$coeff.ind   # c(x$coefficients.ind$c1, x$fit_ind$fit$argument[-c(1 : (x$VC$K1 - 1))], 0)
#Vb.ind           <- x$Vb.ind      # rbind(cbind(x$Vb.ind,0),0)

xx <- x
xx$coefficients <- x$coeff.ind
xx$Vb           <- x$Vb.ind

           pkt <- predict(xx, eq = 1, newdata = newdata, type = "response")$p1.cum 
           pk  <- pkt[, y1]
if(y1 > 1) pk1 <- pkt[, y1 - 1]

p1m <- pk - pk1                         
                                                    

eta2 <- predict(xx, eq = 2, newdata = newdata, type = "link")


if( !is.null(xx$X3) && xx$margins[2] %in% c(cont2par, cont3par) ){

sigma2 <- esp.tr(predict(xx, eq = 3, newdata = newdata), xx$margins[2])$vrb

if(xx$margins[2] %in% cont2par)  eq.th <- 4
if(xx$margins[2] %in% cont3par){ eq.nu <- 4; eq.th <- 5} # case of 3 par margin not implemented yet

if(xx$margins[2] %in% cont3par)  nu <- enu.tr(predict(xx, eq = eq.nu, newdata = newdata), xx$margins[2])$vrb

}

if( is.null(xx$X3) ){

   sigma2 <- esp.tr(xx$coeff.ind[length(xx$coeff.ind)], xx$margins[2])$vrb 

}



p2 <- distrHsAT(y2, eta2, sigma2, nu, xx$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2  

if(cond == 0) p12 <- p1m*p2                             
if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1m



#############################
##### intervals == TRUE #####
#############################

if(intervals == TRUE){

# Cut points need to be transformed first...

cut.sim <- xx$coefficients[1 : (xx$VC$K1 - 1)]
	cut.sim.ti <- rep(0, xx$VC$K1 - 1)
	cut.sim.ti[1] <- cut.sim[1] ; for(i in 2 : (xx$VC$K1 - 1)) {cut.sim.ti[i] <- sqrt(cut.sim[i] - cut.sim[i - 1])}

coefficients.s <- xx$coefficients
coefficients.s[1 : (xx$VC$K1 - 1)] <- cut.sim.ti


bs <- rMVN(n.sim, mean = coefficients.s, sigma = xx$Vb)  


# ... and then transformed back

cut.sim.ti <- bs[, 1 : (xx$VC$K1 - 1)]
	cut.sim <- matrix(nrow = n.sim, ncol = xx$VC$K1 - 1, 0)	
	cut.sim[, 1] <- cut.sim.ti[, 1] ; for (i in 2 : (xx$VC$K1 - 1)) cut.sim[, i] <- cut.sim[, i - 1] + cut.sim.ti[, i]^2

bs[, 1 : (xx$VC$K1 - 1)] <- cut.sim


#############  
# etas
############# 


X1 <- predict(xx, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict(xx, eq = 2, newdata = newdata, type = "lpmatrix") 

              
# Adjustments for the ordinal-continuous model

if (!is.null(xx$VC$K1)) {
	CLM.shift  <- xx$VC$K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}


eta1s <- X1 %*% t(bs[, (CLM.shift2 + 1) : (xx$X1.d2 + CLM.shift2)])
n.s   <- 1


lp1s <- lp1s1 <- matrix(nrow = 1, ncol = n.sim, 0)

for (i in 1 : n.sim) {
	c1s.m   <- t(matrix(nrow = xx$VC$K1 - 1, ncol = n.s        , bs   [i, 1 : CLM.shift2])) 
	eta1s.m <-   matrix(nrow = n.s        , ncol = xx$VC$K1 - 1, eta1s[ , i]              )

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	           lp1s[, i]  <- lp1s_i[, y1] 
	if(y1 > 1) lp1s1[, i] <- lp1s_i[, y1 - 1] 

                      }

###


           pks  <- probm(lp1s,  xx$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
if(y1 > 1) pks1 <- probm(lp1s1, xx$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

eta2s <- eta.tr( X2 %*% t(bs[, (xx$X1.d2 + CLM.shift2 + 1) : (xx$X1.d2 + xx$X2.d2 + CLM.shift2)]), xx$VC$margins[2] )


#############  
# sigmas
#############  

if( xx$VC$margins[2] %in% cont2par ){  

      if( is.null(xx$X3) )  sigma2.star <- bs[, CLM.shift2 + xx$X1.d2 + xx$X2.d2 + 1] 
      if(!is.null(xx$X3) ){
      
                       X3 <- predict(xx, eq = 3, newdata = newdata, type = "lpmatrix")
                       sigma2.star <- X3 %*% t(bs[,((xx$X1.d2 + xx$X2.d2 + 1):(xx$X1.d2 + xx$X2.d2 + xx$X3.d2)) + CLM.shift2]) 
  
                          }
  
sigma2 <- esp.tr(sigma2.star, xx$VC$margins[2])$vrb   
    
                         }   

p2s  <- distrHsAT(y2, eta2s, sigma2, nu, xx$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2  

if(cond == 0) p12s <- (pks - pks1)*p2s                             
if(cond == 1) p12s <- p2s 
if(cond == 2) p12s <- (pks - pks1)



       

} # int

rm(xx)

} # indep




list(p12 = p12, p12s = matrix(p12s, 1, length(p12s)), p1 = pk, p2 = p2, p3 = NULL, CItau = CIkt, tau = tau, CItheta = CItheta, theta = theta)


}