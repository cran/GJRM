res.check <- function(x, intervals = FALSE, n.sim = 100, prob.lev = 0.05){

if(x$triv == TRUE ) stop("This function is not suitable for trivariate probit models.")
if(x$VC$ccss == "yes" && intervals == TRUE) stop("Intervals not available yet for residuals from sample selection models.") 
if(x$Model == "ROY" && intervals == TRUE) stop("Intervals not available yet for residuals from Roy models.")  


y1m <- y2m <- y3m <- NA
qr <- qr1 <- qr2 <- NULL


if(x$univar.gamlss == FALSE){###


if(x$surv.flex == TRUE){ ###

if((!is.null(x$type.cens1) && x$type.cens1 != "R") || (!is.null(x$type.cens2) && x$type.cens2 != "R")) stop("Only case of right censoring implemented. Get in touch to check progress. ")


if(x$end.surv != TRUE) par(mfrow = c(1, 2))

if(x$margins[1] %in% x$bl && x$end.surv == FALSE){

	qr1 <- -log(x$fit$p1)
	H1  <- -log(survfit(Surv(qr1, x$cens1) ~ 1,  type = "kaplan-meier", conf.type = "none")$surv)
	
	#hist(qr1, freq = FALSE, xlab = "Cox-Snell residuals")
	#lines(density(qr1, adjust = 2), lwd = 2)
	
	qqplot(qr1, H1, xlab = "Cox-Snell residuals", ylab = "Cumulative hazards of residuals")
	abline(0, 1)

			}

if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) ){


p1  <- distrHsAT(x$y1, x$eta1, x$sigma21, x$nu1, x$margins[1], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2
qr1 <- qnorm(p1)
#hist(qr1, freq = FALSE, xlab = xlab)
#lines(density(qr1, adjust = 2), lwd = 2)

if(intervals == FALSE){qqnorm(qr1); abline(0, 1)}
if(intervals == TRUE) int.rescheck(x, x$VC$margins[1], n.rep = n.sim, prob.lev = prob.lev, y2m = NULL, eq = 1)

                                         }


qr2 <- -log(x$fit$p2)
H2  <- -log(survfit(Surv(qr2, x$cens2) ~ 1,  type = "kaplan-meier", conf.type = "none")$surv)

#hist(qr2, freq = FALSE, xlab = "Cox-Snell residuals")
#lines(density(qr2, adjust = 2), lwd = 2)

qqplot(qr2, H2, xlab = "Cox-Snell residuals", ylab = "Cumulative hazards of residuals")
abline(0, 1)



}###





if(x$surv.flex == FALSE){##

mbin <- c("probit", "logit", "cloglog")


if(x$VC$margins[1] %in% mbin &&  x$VC$margins[2] %in% mbin ) stop("Diagnostic plot(s) not available given the chosen margin(s).")


if(x$Cont == "NO"){


y2     <- x$y2
eta2   <- x$eta2
sigma2 <- x$sigma2
nu     <- x$nu


if(x$VC$ccss == "yes"){ # for the time being, to be coherent with eta1, I use eta2 etc. calculated from X2s etc. (more relevant for splines)
                        # same reasoning for Roy models
   
   p1  <- x$p1[x$inde]                                            
   p0  <- mm( 1 - p1, min.pr = 1e-300, max.pr = 0.9999999999999999)     
   
   eta2   <- eta2[x$inde]
   theta  <- x$theta
   
   if(length(theta) > 1) theta <- theta[x$inde]

   if(!is.null(x$X3) && !(x$VC$margins[2] %in% x$VC$m1d) ){
   sigma2 <- sigma2[x$inde]
   nu     <- nu[x$inde]
                                                          }
                                                          
                      }
                      
                      




if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3) && x$VC$ccss != "yes" && x$Model != "ROY"){

                                            p <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                              if(x$VC$margins[2] %in% c("TW")){ if( any(y2 == 0) == TRUE ) p[y2 == 0] <- runif(sum(y2 == 0), min = 1e-300, max = p[y2 == 0]) } 
                                            qr <- qnorm(p)
 
                                            } 
                                            
                                            
if(x$VC$margins[2] %in% c(x$VC$m2, x$VC$m3) && x$VC$ccss == "yes" && x$Model != "ROY"){

                                            p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0, p2, x$nC, theta, x$dof), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm( (p2 - p12)/p1, min.pr = 1e-300, max.pr = 0.9999999999999999)                                                  
                                              if(x$VC$margins[2] %in% c("TW")){ if( any(y2 == 0) == TRUE ) p[y2 == 0] <- runif(sum(y2 == 0), min = 1e-300, max = p[y2 == 0]) }
                                            qr  <- qnorm(p)
                                                                                        
                                            } 
                                            
                                            
if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d) && x$VC$ccss != "yes" && x$Model != "ROY"){ 
                                            pd <- distrHsATDiscr(y2, eta2, sigma2, nu, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                 max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2) 
                                            p <- pd$p2
                                            d <- pd$pdf2  
                                              if(intervals == TRUE) set.seed(100)
                                            qr <- qnorm( runif(y2, p - d, p))                                           
                                            } 
                                            
                                            
if(x$VC$margins[2] %in% c(x$VC$m1d, x$VC$m2d, x$VC$m3d) && x$VC$ccss == "yes" && x$Model != "ROY"){ 

                                            p2  <- distrHsATDiscr(y2, eta2, sigma2, nu, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0, p2, x$nC, theta, x$dof), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm( (p2 - p12)/p1, min.pr = 1e-300, max.pr = 0.9999999999999999)  
                                            
                                            y2a <- y2; y2a <- y2 - 1; ind.y2a <- y2a < 0; y2a <- ifelse(y2a < 0, 0, y2a)
                                            p2  <- distrHsATDiscr(y2a, eta2, sigma2, nu, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2 
                                            p12 <- mm(BiCDF(p0, p2, x$nC, theta, x$dof), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            pm1 <- mm( (p2 - p12)/p1, min.pr = 1e-300, max.pr = 0.9999999999999999)   
                                            pm1[ind.y2a] <- 0  
                                            
                                            pm1 <- ifelse(pm1 < 0, 0, pm1) # these needed for very bad model fits
                                            p   <- ifelse(p   < 0, 0, p  )
                                            pm1 <- ifelse(pm1 > p, p, pm1)
                                            ru  <- runif(y2, pm1, p)
                                            ru  <- ifelse(ru == 0, 1e-16, ru)      
                                            ru  <- ifelse(ru > 0.999999, 0.999999, ru)                                     
                                            qr <- qnorm( ru)  
                                            
                                            }                                             
                                            
                                             
   
if(x$Model == "ROY"){                                           

   y2 <- x$y2
   y3 <- x$y3
   
   eta2 <- x$eta2[x$inde0]
   eta3 <- x$eta3[x$inde1]
   
   sigma1 <- x$sigma2[x$inde0]
   sigma2 <- x$sigma3[x$inde1]

   nu1 <- x$nu2[x$inde0]
   nu2 <- x$nu3[x$inde1] 
   
   theta1 <- x$theta12[x$inde0]
   theta2 <- x$theta13[x$inde1]
   
   p1eq2 <- x$p1[x$inde0]  
   p0eq2 <- 1 - p1eq2
   p1eq3 <- x$p1[x$inde1]   
   p0eq3 <- 1 - p1eq3 
   

if(x$VC$margins[2] %in% c(x$VC$m2, x$VC$m3)){

                                            p2  <- distrHsAT(y2, eta2, sigma1, nu1, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2
                                            p12 <- mm(BiCDF(p0eq2, p2, x$nC1, theta1, x$dof12), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)
                                            p   <- mm(p12/p0eq2, min.pr = 1e-300, max.pr = 0.9999999999999999)                                              
                                              if(x$VC$margins[2] %in% c("TW")){ if( any(y2 == 0) == TRUE ) p[y2 == 0] <- runif(sum(y2 == 0), min = 1e-300, max = p[y2 == 0]) }
                                            qr  <- qnorm(p)
                                                                                        
                                            }

if(x$VC$margins[3] %in% c(x$VC$m2, x$VC$m3)){

                                            p2  <- distrHsAT(y3, eta3, sigma2, nu2, x$margins[3], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0eq3, p2, x$nC2, theta2, x$dof13), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm((p2 - p12)/p1eq3, min.pr = 1e-300, max.pr = 0.9999999999999999)                                             
                                              if(x$VC$margins[3] %in% c("TW")){ if( any(y3 == 0) == TRUE ) p[y3 == 0] <- runif(sum(y3 == 0), min = 1e-300, max = p[y3 == 0]) }
                                            qr2 <- qnorm(p)
                                                                                        
                                            }                                             
                                            
if(x$VC$margins[2] %in% c(x$VC$m1d, x$VC$m2d, x$VC$m3d)){ 

                                            p2  <- distrHsATDiscr(y2, eta2, sigma1, nu1, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2
                                            p12 <- mm(BiCDF(p0eq2, p2, x$nC1, theta1, x$dof12), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm(p12/p0eq2, min.pr = 1e-300, max.pr = 0.9999999999999999) 
                                            
                                            y2a <- y2; y2a <- y2 - 1; ind.y2a <- y2a < 0; y2a <- ifelse(y2a < 0, 0, y2a)
                                            p2  <- distrHsATDiscr(y2a, eta2, sigma1, nu1, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2 
                                            p12 <- mm(BiCDF(p0eq2, p2, x$nC1, theta1, x$dof12), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            pm1 <- mm(p12/p0eq2, min.pr = 1e-300, max.pr = 0.9999999999999999)  
                                            pm1[ind.y2a] <- 0  
                                            
                                            pm1 <- ifelse(pm1 < 0, 0, pm1) 
                                            p   <- ifelse(p   < 0, 0, p  )
                                            pm1 <- ifelse(pm1 > p, p, pm1)
                                            ru  <- runif(y2, pm1, p)
                                            ru  <- ifelse(ru == 0, 1e-16, ru)  
                                            ru  <- ifelse(ru > 0.999999, 0.999999, ru)                                                                                
                                            qr <- qnorm(ru)                                              
                                            
                                            }

if(x$VC$margins[3] %in% c(x$VC$m1d, x$VC$m2d, x$VC$m3d)){ 

                                            p2  <- distrHsATDiscr(y3, eta3, sigma2, nu2, x$margins[3], y2m = y3m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0eq3, p2, x$nC2, theta2, x$dof13), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm((p2 - p12)/p1eq3, min.pr = 1e-300, max.pr = 0.9999999999999999) 
                                            
                                            y3a <- y3; y3a <- y3 - 1; ind.y3a <- y3a < 0; y3a <- ifelse(y3a < 0, 0, y3a)
                                            p2  <- distrHsATDiscr(y3a, eta3, sigma2, nu2, x$margins[3], y2m = y3m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2 
                                            p12 <- mm(BiCDF(p0eq3, p2, x$nC2, theta2, x$dof13), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            pm1 <- mm((p2 - p12)/p1eq3, min.pr = 1e-300, max.pr = 0.9999999999999999)  
                                            pm1[ind.y3a] <- 0  
                                            
                                            pm1 <- ifelse(pm1 < 0, 0, pm1) 
                                            p   <- ifelse(p   < 0, 0, p  )  
                                            pm1 <- ifelse(pm1 > p, p, pm1)
                                            ru  <- runif(y3, pm1, p)
                                            ru  <- ifelse(ru == 0, 1e-16, ru)   
                                            ru  <- ifelse(ru > 0.999999, 0.999999, ru)                                            
                                            qr2 <- qnorm(ru)  
                                            
                                            }  


                                          
}                                            
                                            
                                            
                                                    


if(x$Model != "ROY"){  

#par(mfrow = c(1, 1))

#hist(qr, freq = FALSE, xlab = xlab)
#lines(density(qr, adjust = 2), lwd = 2)

if(intervals == FALSE){qqnorm(qr); abline(0, 1)}
if(intervals == TRUE) int.rescheck(x, x$VC$margins[2], n.rep = n.sim, prob.lev = prob.lev, y2m = y2m, eq = 1)

}  


if(x$Model == "ROY"){  


   par(mfrow = c(1, 2))

   #hist(qr, freq = FALSE, xlab = xlab)
   #lines(density(qr, adjust = 2), lwd = 2)
   qqnorm(qr); abline(0, 1)

   #hist(qr2, freq = FALSE, xlab = xlab)
   #lines(density(qr2, adjust = 2), lwd = 2)
   qqnorm(qr2); abline(0, 1)
   
                     } # ROY

#if(intervals == FALSE){qqnorm(qr); abline(0, 1)}
#if(intervals == TRUE){xx <- x; xx$n <- x$n.se0; xx$eta2 <- x$eta2[x$inde0]; xx$sigma2 <- x$sigma1[x$inde0]; xx$nu <- x$nu1[x$inde0]; int.postcheck(xx, xx$VC$margins[2], n.rep = n.sim, prob.lev = prob.lev, y2m = y2m); rm(xx)}
#if(intervals == FALSE){}
#
#if(intervals == TRUE){
#
#xx   <- x
#xx$n <- x$n.se1
#
#xx$eta2   <- x$eta3[x$inde1]
#xx$y2     <- x$y3
#xx$sigma2 <- x$sigma2[x$inde1]
#xx$nu     <- x$nu2[x$inde1]
#
#int.postcheck(xx, x$VC$margins[3], n.rep = n.sim, prob.lev = prob.lev, y2m = y3m)
#
#rm(xx)
#
#}


  
  




} # x$Cont == "NO"




if(x$Cont == "YES"){

y1 <- x$y1
y2 <- x$y2

y1m <- y2m <- NA




if(x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3)) {p1 <- distrHsAT(x$y1, x$eta1, x$sigma21, x$nu1, x$margins[1], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2 

                                            if(x$VC$margins[1] %in% c("TW")){ if( any(x$y1 == 0) == TRUE ) p1[x$y1 == 0] <- runif(sum(x$y1 == 0), min = 1e-300, max = p1[x$y1 == 0]) }


}


if(x$VC$margins[1] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){

pd <- distrHsATDiscr(x$y1, x$eta1, x$sigma21, x$nu1, x$margins[1], y2m = y1m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                     max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1) 
p <- pd$p2
d <- pd$pdf2   

set.seed(100)
p1 <- runif(y1, p - d, p) 

}



if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3)){   p2 <-      distrHsAT(x$y2, x$eta2, x$sigma22, x$nu2, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2 

                                            if(x$VC$margins[2] %in% c("TW")){ if( any(x$y2 == 0) == TRUE ) p2[x$y2 == 0] <- runif(sum(x$y2 == 0), min = 1e-300, max = p2[x$y2 == 0])   }

}

if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){

pd <- distrHsATDiscr(x$y2, x$eta2, x$sigma22, x$nu2, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                     max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)
p <- pd$p2
d <- pd$pdf2   

set.seed(100)
p2 <- runif(y2, p - d, p) 

}


par(mfrow = c(1, 2))

qr1 <- qnorm(p1)
#hist(qr1, freq = FALSE, xlab = xlab)
#lines(density(qr1, adjust = 2), lwd = 2)


if(intervals == FALSE){qqnorm(qr1); abline(0, 1)}
if(intervals == TRUE) int.rescheck(x, x$VC$margins[1], n.rep = n.sim, prob.lev = prob.lev, y2m = y1m, eq = 1)


qr2 <- qnorm(p2)
#hist(qr2, freq = FALSE, xlab = xlab2)
#lines(density(qr2, adjust = 2),lwd=2)

if(intervals == FALSE){qqnorm(qr2); abline(0, 1)}
if(intervals == TRUE) int.rescheck(x, x$VC$margins[2], n.rep = n.sim, prob.lev = prob.lev, y2m = y2m, eq = 2)

}


}##

}###




if(x$univar.gamlss == TRUE){

if(x$surv.flex == FALSE){ ###

if(x$VC$margins[1] %in% c("GEVlink") ) stop("Diagnostic plot(s) not available given the chosen margin(s).")

y1 <- x$y1

                         


if(x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3)) p1 <- distrHsAT(x$y1, x$eta1, x$sigma2, x$nu, x$margins[1], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc)$p2 

        if(x$VC$margins[1] %in% c("TW")){ if( any(x$y1 == 0) == TRUE )  p1[x$y1 == 0] <- runif(sum(x$y1 == 0), min = 1e-300, max = p1[x$y1 == 0]) }

if(x$VC$margins[1] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){
      pd <- distrHsATDiscr(x$y1, x$eta1, x$sigma2, x$nu, x$margins[1], y2m = y1m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                           max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc) 
      p <- pd$p2
      d <- pd$pdf2   
      p1 <- runif(y1, p - d, p) 
                                                      }

#par(mfrow = c(1, 1))

qr <- qnorm(p1)
#hist(qr, freq = FALSE, xlab = xlab)
#lines(density(qr, adjust = 2), lwd = 2)

if(intervals == FALSE){qqnorm(qr); abline(0, 1)}
if(intervals == TRUE) int.rescheck(x, x$VC$margins[1], n.rep = n.sim, prob.lev = prob.lev, y2m = y1m, eq = 1)

}###


if(x$surv.flex == TRUE){ ###


if(x$type.cens != "R") stop("Only case of right censoring implemented. Get in touch to check progress. ")


#par(mfrow = c(1, 1))

qr <- -log(x$fit$p1)
H  <- -log(survfit(Surv(qr, x$cens) ~ 1,  type = "kaplan-meier", conf.type = "none")$surv)

#hist(qr, freq = FALSE, xlab = "Cox-Snell residuals")
#lines(density(qr, adjust = 2), lwd = 2)

qqplot(qr, H, xlab = "Cox-Snell residuals", ylab = "Cumulative hazards of residuals")
abline(0, 1)

}###


}


L <- list(qr = qr, qr1 = qr1, qr2 = as.numeric(qr2))

invisible(L)


}

