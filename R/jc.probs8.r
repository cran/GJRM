jc.probs8 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = FALSE){

if(!(y1 >= range(x$y1)[1] && y1 <= range(x$y1)[2])) stop("y1 is not within its observed range.") 
if(!(y2 >= range(x$y2)[1] && y2 <= range(x$y2)[2])) stop("y2 is not within its observed range.") 

if(x$BivD %in% x$BivD2) stop("The chosen copula is not supported by this function.")

# Ordinal-Ordinal case

CIp12 <- CIkt <- tau <- theta <- CItheta <- NULL
dof   <- x$dof
infty <- 1e+25
p12.f <- p12.s <- p12.t <- p12.z <- 0
p12sf <- p12ss <- p12st <- p12sz <- 0

n1.s <- n2.s <- 1

pk1s1 <- pk2s1 <- pk11 <- pk22 <- 0

theta <- NULL

p12s  <- matrix(0, 1, 2) 


if(type == "joint"){

           pk1t  <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum   
           pk1   <- pk1t[, y1]                                                       
if(y1 > 1) pk11  <- pk1t[, y1 - 1]

           pk2t  <- predict(x, eq = 2, newdata = newdata, type = "response")$p2.cum  
           pk2   <- pk2t[, y2] 						             
if(y2 > 1) pk22  <- pk2t[, y2 - 1] 


p1m <- pk1 - pk11                         
p2m <- pk2 - pk22     
                       

if(!is.null(x$X3)) theta <- teta.tr(x$VC, predict(x, eq = 3, newdata = newdata, type = "link"))$teta
if( is.null(x$X3)) theta <- x$theta 


if(y1 > 1 && y2 > 1){

p12.f <- mm(BiCDF(pk1,  pk2,  x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 
p12.s <- mm(BiCDF(pk11, pk2,  x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 
p12.t <- mm(BiCDF(pk1,  pk22, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 
p12.z <- mm(BiCDF(pk11, pk22, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 

                    }

if(y1 == 1 && y2 == 1){

p12.f <- mm(BiCDF(pk1,  pk2,  x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 

                       }


if(y1 > 1 && y2 == 1){

p12.f <- mm(BiCDF(pk1,  pk2,  x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 
p12.s <- mm(BiCDF(pk11, pk2,  x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 

                     }


if(y1 == 1 && y2 > 1){

p12.f <- mm(BiCDF(pk1,  pk2,  x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 
p12.t <- mm(BiCDF(pk1,  pk22, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) 

                     }

p12 <- p12.f - p12.s - p12.t + p12.z   

if(cond == 1) p12 <- p12/p1m
if(cond == 2) p12 <- p12/p2m

#######################################################################################
### kendalls' tau and theta ###

if(tau.res == TRUE) tau   <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = tau.res)$tau
                    theta <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = FALSE)$theta
     
     
     
#############################
##### intervals == TRUE ##### 
#############################

if(intervals == TRUE){



# Cut points need to be transformed first... not non increasing sequence and then sample

cut1.sim <- x$coefficients[1 : (x$VC$K1 - 1)]
	cut1.sim.ti <- rep(0, x$VC$K1 - 1)
	cut1.sim.ti[1] <- cut1.sim[1] ; for(i in 2 : (x$VC$K1 - 1)) {cut1.sim.ti[i] <- sqrt(cut1.sim[i] - cut1.sim[i - 1])}

cut2.sim <- x$coefficients[x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)]
	cut2.sim.ti <- rep(0, x$VC$K2 - 1)
	cut2.sim.ti[1] <- cut2.sim[1] ; for(i in 2 : (x$VC$K2 - 1)) {cut2.sim.ti[i] <- sqrt(cut2.sim[i] - cut2.sim[i - 1])}

coefficients.s <- x$coefficients
	coefficients.s[1 : (x$VC$K1 - 1)]                 <- cut1.sim.ti
	coefficients.s[x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)] <- cut2.sim.ti




bs <- rMVN(n.sim, mean = coefficients.s, sigma = x$Vb)  

# ... and then transformed back to increasing sequence as this is what predict is based on
# x$Vb is on non increasing scale

cut1.sim.ti <- bs[, 1 : (x$VC$K1 - 1)]
	cut1.sim <- matrix(nrow = n.sim, ncol = x$VC$K1 - 1, 0)	
	cut1.sim[, 1] <- cut1.sim.ti[, 1] ; for (i in 2 : (x$VC$K1 - 1)) cut1.sim[, i] <- cut1.sim[, i - 1] + cut1.sim.ti[, i]^2

cut2.sim.ti <- bs[, x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)]
	cut2.sim <- matrix(nrow = n.sim, ncol = x$VC$K2 - 1, 0)	
	cut2.sim[, 1] <- cut2.sim.ti[, 1] ; for (i in 2 : (x$VC$K2 - 1)) cut2.sim[, i] <- cut2.sim[, i - 1] + cut2.sim.ti[, i]^2

bs[, 1 : (x$VC$K1 + x$VC$K2 - 2)] <- cbind(cut1.sim, cut2.sim)


################


X1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix") 

              
# Adjustments for the ordinal-ordinal model

CLM.shift  <- x$VC$K1 + x$VC$K2 - 2

eta1s <- X1 %*% t(bs[, (CLM.shift + 1) : (x$X1.d2 + CLM.shift)])
eta2s <- X2 %*% t(bs[, (x$X1.d2 + CLM.shift + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift)])


lp1s <- lp2s <- lp1s1 <- lp2s1 <- matrix(nrow = 1, ncol = n.sim, 0)


for (i in 1 : n.sim) {

	c1s.m   <- t(matrix(nrow = x$VC$K1 - 1, ncol = 1, bs[i, 1 : (x$VC$K1 - 1)])) 
	c2s.m   <- t(matrix(nrow = x$VC$K2 - 1, ncol = 1, bs[i, x$VC$K1 : CLM.shift])) 

	eta1s.m <- matrix(nrow = 1, ncol = x$VC$K1 - 1, eta1s[, i])
	eta2s.m <- matrix(nrow = 1, ncol = x$VC$K2 - 1, eta2s[, i])

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	lp2s_i <- cbind(c2s.m - eta2s.m, infty)

	           lp1s[, i]  <- lp1s_i[, y1]                                                          
	if(y1 > 1) lp1s1[, i] <- lp1s_i[, y1-1]                                                        
	           lp2s[, i]  <- lp2s_i[, y2]                                                          
	if(y2 > 1) lp2s1[, i] <- lp2s_i[, y2-1]                                                        
	
                      }

### 

           pk1s  <- probm(lp1s,  x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
           pk2s  <- probm(lp2s,  x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

if(y1 > 1) pk1s1 <- probm(lp1s1, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
if(y2 > 1) pk2s1 <- probm(lp2s1, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

##

if(  is.null(x$X3) ) epds <- bs[,length(x$coefficients)]

if( !is.null(x$X3) ){ 
	             X3   <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")                           
	             epds <- X3%*%t(bs[,((x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)) + CLM.shift])
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


##



################# 

if(y1 > 1 && y2 > 1){

p12sf <- mm(BiCDF(pk1s,  pk2s,  x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
p12ss <- mm(BiCDF(pk1s1, pk2s,  x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
p12st <- mm(BiCDF(pk1s,  pk2s1, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
p12sz <- mm(BiCDF(pk1s1, pk2s1, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
	
}	
	

if(y1 == 1 && y2 == 1){

p12sf <- mm(BiCDF(pk1s,  pk2s,  x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)      
	
}	
		

if(y1 > 1 && y2 == 1){

p12sf <- mm(BiCDF(pk1s,  pk2s,  x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
p12ss <- mm(BiCDF(pk1s1, pk2s,  x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)      
	
}


if(y1 == 1 && y2 > 1){

p12sf <- mm(BiCDF(pk1s,  pk2s,  x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
p12st <- mm(BiCDF(pk1s,  pk2s1, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
	
}	
	
	
	
p12s <- p12sf - p12ss - p12st + p12sz   	
	
	
if(cond == 1) p12s <- p12s/(pk1s - pk1s1)
if(cond == 2) p12s <- p12s/(pk2s - pk2s1)



} # int




     
} # biv





###########################################
###########################################
###########################################


if(type == "independence"){ 

#coefficients.ind <- x$coeff.ind # c(x$coefficients.ind$c1, x$coefficients.ind$c2, x$coefficients.ind$beta1, x$coefficients.ind$beta2, 0)
#Vb.ind <- x$Vb.ind              # rbind(cbind(x$Vb.ind,0),0)

xx <- x
xx$coefficients <- x$coeff.ind # coefficients.ind
xx$Vb           <- x$Vb.ind    # Vb.ind

           pk1t  <- predict(xx, eq = 1, newdata = newdata, type = "response")$p1.cum   
           pk1   <- pk1t[, y1]                                                       
if(y1 > 1) pk11  <- pk1t[, y1 - 1]

           pk2t  <- predict(xx, eq = 2, newdata = newdata, type = "response")$p2.cum  
           pk2   <- pk2t[, y2] 						             
if(y2 > 1) pk22  <- pk2t[, y2 - 1] 

p1m <- pk1 - pk11                        
p2m <- pk2 - pk22     
                      
p12 <- p1m*p2m  

if(cond == 1) p12 <- p2m
if(cond == 2) p12 <- p1m

#############################
##### intervals == TRUE ##### 
#############################

if(intervals == TRUE){

# Cut points need to be transformed first...

cut1.sim <- xx$coefficients[1 : (xx$VC$K1 - 1)]
	cut1.sim.ti <- rep(0, xx$VC$K1 - 1)
	cut1.sim.ti[1] <- cut1.sim[1] ; for(i in 2 : (xx$VC$K1 - 1)) {cut1.sim.ti[i] <- sqrt(cut1.sim[i] - cut1.sim[i - 1])}

cut2.sim <- xx$coefficients[xx$VC$K1 : (xx$VC$K1 + xx$VC$K2 - 2)]
	cut2.sim.ti <- rep(0, xx$VC$K2 - 1)
	cut2.sim.ti[1] <- cut2.sim[1] ; for(i in 2 : (xx$VC$K2 - 1)) {cut2.sim.ti[i] <- sqrt(cut2.sim[i] - cut2.sim[i - 1])}

coefficients.s <- xx$coefficients
	coefficients.s[1 : (xx$VC$K1 - 1)]                   <- cut1.sim.ti
	coefficients.s[xx$VC$K1 : (xx$VC$K1 + xx$VC$K2 - 2)] <- cut2.sim.ti


# ... and then transformed back

bs <- rMVN(n.sim, mean = coefficients.s, sigma = xx$Vb)  


cut1.sim.ti <- bs[, 1 : (xx$VC$K1 - 1)]
	cut1.sim <- matrix(nrow = n.sim, ncol = xx$VC$K1 - 1, 0)	
	cut1.sim[, 1] <- cut1.sim.ti[, 1] ; for (i in 2 : (xx$VC$K1 - 1)) cut1.sim[, i] <- cut1.sim[, i - 1] + cut1.sim.ti[, i]^2

cut2.sim.ti <- bs[, xx$VC$K1 : (xx$VC$K1 + xx$VC$K2 - 2)]
	cut2.sim <- matrix(nrow = n.sim, ncol = xx$VC$K2 - 1, 0)	
	cut2.sim[, 1] <- cut2.sim.ti[, 1] ; for (i in 2 : (xx$VC$K2 - 1)) cut2.sim[, i] <- cut2.sim[, i - 1] + cut2.sim.ti[, i]^2

bs[, 1 : (xx$VC$K1 + xx$VC$K2 - 2)] <- cbind(cut1.sim, cut2.sim)


################



X1 <- predict(xx, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict(xx, eq = 2, newdata = newdata, type = "lpmatrix") 

              
# Adjustments for the ordinal-ordinal model

CLM.shift  <- xx$VC$K1 + xx$VC$K2 - 2

eta1s <- X1 %*% t(bs[, (CLM.shift + 1) : (xx$X1.d2 + CLM.shift)])
eta2s <- X2 %*% t(bs[, (xx$X1.d2 + CLM.shift + 1) : (xx$X1.d2 + xx$X2.d2 + CLM.shift)])


lp1s  <- lp2s <- lp1s1 <- lp2s1 <- matrix(nrow = 1, ncol = n.sim, 0)


for (i in 1 : n.sim) {

	c1s.m   <- t(matrix(nrow = xx$VC$K1 - 1, ncol = 1, bs[i, 1 : (xx$VC$K1 - 1)])) 
	c2s.m   <- t(matrix(nrow = xx$VC$K2 - 1, ncol = 1, bs[i, xx$VC$K1 : CLM.shift])) 

	eta1s.m <- matrix(nrow = 1, ncol = xx$VC$K1 - 1, eta1s[, i])
	eta2s.m <- matrix(nrow = 1, ncol = xx$VC$K2 - 1, eta2s[, i])

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	lp2s_i <- cbind(c2s.m - eta2s.m, infty)

	           lp1s[, i]  <- lp1s_i[, y1]                                                          
	if(y1 > 1) lp1s1[, i] <- lp1s_i[, y1-1]                                                        
	           lp2s[, i]  <- lp2s_i[, y2]                                                          
	if(y2 > 1) lp2s1[, i] <- lp2s_i[, y2-1]                                                        
	
}

### 

pk1s  <- probm(lp1s, xx$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
pk2s  <- probm(lp2s, xx$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

if(y1 > 1) pk1s1 <- probm(lp1s1, xx$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
if(y2 > 1) pk2s1 <- probm(lp2s1, xx$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

p1ms <- pk1s - pk1s1                
p2ms <- pk2s - pk2s1

p12s <- p1ms*p2ms
		
if(cond == 1) p12s <- p2ms
if(cond == 2) p12s <- p1ms



} # int

rm(xx)

} # indep


list(p12 = p12, p12s = matrix(p12s, 1, length(p12s)), p1 = pk1, p2 = pk2, p3 = NULL, CItau = CIkt, tau = tau, CItheta = CItheta, theta = theta)


}