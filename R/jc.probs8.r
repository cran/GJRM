jc.probs8 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr){

# Prediction for ordinal-ordinal case

nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- NULL
CIkt <- tau <- NULL
dof <- x$dof

infty <- 1e+25

p12s <- NULL


if(type == "joint"){


##### newdata included as inputs #####

if(!missing(newdata)){

nu <- sigma2 <- NA

#type <- "response"
pk1 <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum # cdf
	pk1 <- diag(pk1[, y1]) # pk1 <- pk1[, y1] # @ GIAMPIERO: I think this is what you want... maybe we should change it in jc.probs7
pk2 <- predict(x, eq = 2, newdata = newdata, type = "response")$p2.cum # cdf
	pk2 <- diag(pk2[, y2]) # pk2 <- pk2[, y2]

#eta2 <- predict(x, eq = 2, newdata = newdata, type = "response")


if (!is.null(x$X3)) {
	sigma2 <- NA #
	nu     <- NA # These parameters are non supported by the ordinal-ordinal model
	theta <- teta.tr(x$VC, predict(x, eq = 3, newdata = newdata, type = "response"))$teta
}



if( is.null(x$X3) ){

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}
#
#
}


##### newdata NOT included as inputs #####

if(missing(newdata)){

infty <- 1e+25

#cut <- t(matrix(nrow = x$VC$K1 - 1, ncol = x$n, x$coefficients[1 : (x$VC$K1 - 1)])) #
#eta1 <- matrix(nrow = x$n, ncol = x$VC$K1 - 1, x$eta1) ##############################

lp1 <- cbind(x$fit$lp1, infty)
lp2 <- cbind(x$fit$lp2, infty) 

lp1.sel <- diag(lp1[, y1]) # lp1.sel <- lp1[, y1] # The choice of y1 selects the relevant cut point
lp2.sel <- diag(lp2[, y2]) # lp2.sel <- lp2[, y2] # The choice of y1 selects the relevant cut point

pk1 <- probm(lp1.sel, x$VC$margins[1], only.pr = FALSE, bc = TRUE, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr # cumulative distribution function
	pk1[y1 == x$VC$K1] <- 1
pk2 <- probm(lp2.sel, x$VC$margins[2], only.pr = FALSE, bc = TRUE, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr # cumulative distribution function
	pk2[y2 == x$VC$K2] <- 1


#p1   <- x$p1
#eta2 <- x$eta2

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

} 


##### cond == 0, cond == 1 and cond == 2 ##### QUESTION TO GM: I had to re-write bits below for the ordinal-ordinal case

p0 <- pk1
p2 <- pk2

if(!(x$BivD %in% x$BivD2)) p12 <- mm(BiCDF(p0, p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr) # This is a bivariate cdf

if(x$BivD %in% x$BivD2) {
	nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1), 2] 
	nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2), 2]
	p12 <- NA

	if(length(x$teta1) != 0) {
		if(length(theta) > 1 ) p12[x$teta.ind1] <- mm(BiCDF(p0[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr)
		if(length(theta) == 1) p12[x$teta.ind1] <- mm(BiCDF(p0[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof), min.pr = min.pr, max.pr = max.pr)
	}

	if( length(x$teta2) != 0) {
		if(length(theta) > 1 ) p12[x$teta.ind2] <- mm(BiCDF(p0[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr)
		if(length(theta) == 1) p12[x$teta.ind2] <- mm(BiCDF(p0[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof), min.pr = min.pr, max.pr = max.pr)
	}                       
}

if (cond == 1) p12 <- p12 / p0
if (cond == 2) p12 <- p12 / p2









# kendalls' tau

if(x$BivD %in% x$BivD2)   {x$SemiParFit <- x; tau <- Reg2Copost(x$SemiParFit, x$VC, theta)$tau }
if(!(x$BivD %in% x$BivD2)) tau <- ass.ms(x$VC$BivD, x$VC$nCa, theta)$tau 


#############################
##### intervals == TRUE ##### 
#############################

if(intervals == TRUE){

# Cut points need to be transformed first...

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


# ... and then transformed back

cut1.sim.ti <- bs[, 1 : (x$VC$K1 - 1)]
	cut1.sim <- matrix(nrow = n.sim, ncol = x$VC$K1 - 1, 0)	
	cut1.sim[, 1] <- cut1.sim.ti[, 1] ; for (i in 2 : (x$VC$K1 - 1)) cut1.sim[, i] <- cut1.sim[, i - 1] + cut1.sim.ti[, i]^2

cut2.sim.ti <- bs[, x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)]
	cut2.sim <- matrix(nrow = n.sim, ncol = x$VC$K2 - 1, 0)	
	cut2.sim[, 1] <- cut2.sim.ti[, 1] ; for (i in 2 : (x$VC$K2 - 1)) cut2.sim[, i] <- cut2.sim[, i - 1] + cut2.sim.ti[, i]^2

bs[, 1 : (x$VC$K1 + x$VC$K2 - 2)] <- cbind(cut1.sim, cut2.sim)


#############  
# etas
############# 

##### newdata included as inputs #####

#type <- "lpmatrix"

if(!missing(newdata)){ 
	pred.1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
		X1  <- pred.1 ; colnames(X1) <- colnames(x$X1)
		X2s <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix") 
}
              

##### newdata NOT included as inputs #####
         
if( missing(newdata)){ X1 <- x$X1
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       }


### 

# Adjustments for the ordinal-ordinal model

CLM.shift  <- x$VC$K1 + x$VC$K2 - 2

eta1s <- X1  %*% t(bs[, (CLM.shift + 1) : (x$X1.d2 + CLM.shift)])
eta2s <- X2s %*% t(bs[, (x$X1.d2 + CLM.shift + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift)])

n1.s <- dim(X1)[1]
n2.s <- dim(X2s)[1]

lp1s <- matrix(nrow = n1.s, ncol = n.sim, 0)
lp2s <- matrix(nrow = n2.s, ncol = n.sim, 0)

for (i in 1 : n.sim) {
	c1s.m   <- t(matrix(nrow = x$VC$K1 - 1, ncol = n1.s, bs[i, 1 : (x$VC$K1 - 1)])) 
	c2s.m   <- t(matrix(nrow = x$VC$K2 - 1, ncol = n2.s, bs[i, x$VC$K1 : CLM.shift])) 

	eta1s.m <-   matrix(nrow = n1.s, ncol = x$VC$K1 - 1, eta1s[, i])
	eta2s.m <-   matrix(nrow = n2.s, ncol = x$VC$K2 - 1, eta2s[, i])

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	lp2s_i <- cbind(c2s.m - eta2s.m, infty)

	lp1s[, i] <- diag(lp1s_i[, y1]) # The choice of y1 selects the relevant cut point
	lp2s[, i] <- diag(lp2s_i[, y2]) # The choice of y2 selects the relevant cut point
}

### 

pk1s <- probm(lp1s, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
pk2s <- probm(lp2s, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
#eta2s <- eta.tr( X2s %*% t(bs[, (x$X1.d2 + CLM.shift2 + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift2)]), x$VC$margins[2] )

#############  
# thetas
#############  

if( is.null(x$X3) ) epds <- bs[,length(x$coefficients)]

if ( !is.null(x$X3) ) { 
	if(!missing(newdata)) { X3s <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix") }
	if( missing(newdata)) { if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3 }    
	epds <- X3s%*%t(bs[,((x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)) + CLM.shift])
}




est.RHOb <- teta.tr(x$VC, epds)$teta
   

    
    




################# 


if( is.null(x$X3) ){



if(is.null(sigma2)) sigma2 <- rep(1, length(est.RHOb)) 
if(is.null(nu))     nu     <- rep(1, length(est.RHOb))


est.RHOb <- matrix(rep(est.RHOb, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }


##### cond == 0, cond == 1 and cond == 2 ##### @ GIAMPIERO: as before, I have replaced the parts commented with the bits below

p0s <- pk1s
p2s <- pk2s

if (x$VC$BivD %in% c("N","T")) {                                                 # @ GIAMPIERO: this line is not included in the previous bits on cond...
	p12s <- matrix(mm(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  ), dim(pk1s)[1], n.sim) 
} else {
	if (x$BivD %in% x$BivD2) {
		p12s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])
		if(length(x$teta1) != 0) p12s[x$teta.ind1, ] <- mm(BiCDF(p0s[x$teta.ind1, ], p2s[x$teta.ind1, ], nC1,  est.RHOb[x$teta.ind1, ]), min.pr = min.pr, max.pr = max.pr)                  
		if(length(x$teta2) != 0) p12s[x$teta.ind2, ] <- mm(BiCDF(p0s[x$teta.ind2, ], p2s[x$teta.ind2, ], nC2, -est.RHOb[x$teta.ind2, ]), min.pr = min.pr, max.pr = max.pr)            
	}

	if(!(x$BivD %in% x$BivD2)) {
		p12s <- mm(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr)    
		p12s <- matrix(nrow = n1.s, ncol = n.sim, p12s) # This part is not in jc.probs2
	}

	if (cond == 1) p12s <- p12s / p0s
	if (cond == 2) p12s <- p12s / p2s
} 








 
   
   
   
   
  
   
   
   

nCa   <- x$VC$nCa
BivDt <- x$VC$BivD

  if(x$BivD %in% x$BivD2){
  
  if(x$BivD %in% x$BivD2[c(1:4,13:16)]) { BivDt <- "C0"; nCa <- 3} 
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

   
   
   
   
   
   

} # int




     
} # biv





###########################################
###########################################
###########################################


if(type == "independence"){ 

##### newdata included as inputs ##### 

if(!missing(newdata)){

nu  <- sigma2 <- NA
pk1 <- probm( predict(x, eq = 1, newdata = newdata, type = "lpmatrix") %*% x$coefficients.ind$beta1, x$VC$margin[1], 
	      min.dn = min.pr, min.pr = min.pr, max.pr = max.pr )$pr # cdf
pk2 <- probm( predict(x, eq = 2, newdata = newdata, type = "lpmatrix") %*% x$coefficients.ind$beta2, x$VC$margin[2], 
	      min.dn = min.pr, min.pr = min.pr, max.pr = max.pr )$pr # cdf
#pk1 <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum # cdf
#	pk1 <- diag(pk1[, y1]) # pk1 <- pk1[, y1]
#pk2 <- predict(x, eq = 2, newdata = newdata, type = "response")$p2.cum # cdf
#	pk2 <- diag(pk2[, y2]) # pk2 <- pk2[, y2]
#eta2 <- predict(x, eq = 2, newdata = newdata, type = "response")





} 


##### newdata NOT included as inputs ##### 

if(missing(newdata)){

infty <- 1e+25

#cut <- t(matrix(nrow = x$VC$K1 - 1, ncol = x$n, x$coefficients[1 : (x$VC$K1 - 1)])) #
#eta1 <- matrix(nrow = x$n, ncol = x$VC$K1 - 1, x$eta1) ##############################

eta1 <- x$VC$X1 %*% x$coefficients.ind$beta1 ; eta1.m <- matrix(nrow = x$n, ncol = x$VC$K1 - 1, eta1)
eta2 <- x$VC$X2 %*% x$coefficients.ind$beta2 ; eta2.m <- matrix(nrow = x$n, ncol = x$VC$K2 - 1, eta2)

c1.m <- t(matrix(nrow = x$VC$K1 - 1, ncol = x$n, x$coefficients.ind$c1))
c2.m <- t(matrix(nrow = x$VC$K2 - 1, ncol = x$n, x$coefficients.ind$c2))

lp1 <- cbind(c1.m - eta1.m, infty)
lp2 <- cbind(c2.m - eta2.m, infty) 

lp1.sel <- diag(lp1[, y1]) # lp1.sel <- lp1[, y1] # The choice of y1 selects the relevant cut point
lp2.sel <- diag(lp2[, y2]) # lp2.sel <- lp2[, y2] # The choice of y2 selects the relevant cut point

pk1 <- probm(lp1.sel, x$VC$margins[1], only.pr = FALSE, bc = TRUE, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr # cumulative distribution function
	pk1[y1 == x$VC$K1] <- 1
pk2 <- probm(lp2.sel, x$VC$margins[2], only.pr = FALSE, bc = TRUE, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr # cumulative distribution function
	pk2[y2 == x$VC$K2] <- 1


#p1   <- x$p1
#eta2 <- x$eta2

sigma2 <- x$sigma2
nu     <- x$nu 

} 



##### cond == 0, cond == 1 and cond == 2 ##### 

p0 <- pk1
p2 <- pk2                    
p12 <- p0 * p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p0

        





#############################
##### intervals == TRUE ##### 
############################# 

if(intervals == TRUE){

# Cut points need to be transformed first...

cut1.sim <- x$coefficients.ind$c1
	cut1.sim.ti <- rep(0, x$VC$K1 - 1)
	cut1.sim.ti[1] <- cut1.sim[1] ; for(i in 2 : (x$VC$K1 - 1)) {cut1.sim.ti[i] <- sqrt(cut1.sim[i] - cut1.sim[i - 1])}

cut2.sim <- x$coefficients.ind$c2
	cut2.sim.ti <- rep(0, x$VC$K2 - 1)
	cut2.sim.ti[1] <- cut2.sim[1] ; for(i in 2 : (x$VC$K2 - 1)) {cut2.sim.ti[i] <- sqrt(cut2.sim[i] - cut2.sim[i - 1])}

#coefficients.s <- x$coefficients
#	coefficients.s[1 : (x$VC$K1 - 1)]                 <- cut1.sim.ti
#	coefficients.s[x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)] <- cut2.sim.ti
#coefficients.s <- c(cut1.sim.ti, cut2.sim.ti, x$gam1$coefficients, x$gam2$coefficients)

bs.c1 <- rMVN(n.sim, mean = cut1.sim.ti, sigma = x$Vb.ind[1 : (x$VC$K1 - 1), 1 : (x$VC$K1 - 1)])
bs.c2 <- rMVN(n.sim, mean = cut2.sim.ti, sigma = x$Vb.ind[x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2), x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)])
bs1 <- rMVN(n.sim, mean = x$coefficients.ind$beta1, sigma = x$Vb.ind[(x$VC$K1 + x$VC$K2 - 1) : (x$VC$K1 + x$VC$K2 + x$VC$X1.d2 - 2), (x$VC$K1 + x$VC$K2 - 1) : (x$VC$K1 + x$VC$K2 + x$VC$X1.d2 - 2)]) #
bs2 <- rMVN(n.sim, mean = x$coefficients.ind$beta2, sigma = x$Vb.ind[(x$VC$K1 + x$VC$K2 + x$VC$X1.d2 - 1) : (x$VC$K1 + x$VC$K2 + x$VC$X1.d2 + x$VC$X2.d2 - 2), (x$VC$K1 + x$VC$K2 + x$VC$X1.d2 - 1) : (x$VC$K1 + x$VC$K2 + x$VC$X1.d2 + x$VC$X2.d2 - 2)]) # The intercepts are deleted from the covariance matrices
bs <- cbind(bs.c1, bs.c2, bs1, bs2)
#bs <- rMVN(n.sim, mean = coefficients.s, sigma = x$Vb)  


# ... and then transformed back 

cut1.sim.ti <- bs[, 1 : (x$VC$K1 - 1)]
	cut1.sim <- matrix(nrow = n.sim, ncol = x$VC$K1 - 1, 0)	
	cut1.sim[, 1] <- cut1.sim.ti[, 1] ; for (i in 2 : (x$VC$K1 - 1)) cut1.sim[, i] <- cut1.sim[, i - 1] + cut1.sim.ti[, i]^2

cut2.sim.ti <- bs[, x$VC$K1 : (x$VC$K1 + x$VC$K2 - 2)]
	cut2.sim <- matrix(nrow = n.sim, ncol = x$VC$K2 - 1, 0)	
	cut2.sim[, 1] <- cut2.sim.ti[, 1] ; for (i in 2 : (x$VC$K2 - 1)) cut2.sim[, i] <- cut2.sim[, i - 1] + cut2.sim.ti[, i]^2

bs[, 1 : (x$VC$K1 + x$VC$K2 - 2)] <- cbind(cut1.sim, cut2.sim)


#############  
# etas
#############   


##### newdata included as inputs ##### 

#type <- "lpmatrix"

if(!missing(newdata)){ 
	pred.1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
		X1  <- pred.1 ; colnames(X1) <- colnames(x$X1)
		X2s <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix") 
}
     

##### newdata NOT included as inputs #####
         
if( missing(newdata)){ X1 <- x$X1
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       }


### 

# Adjustments for the ordinal-ordinal model

CLM.shift  <- x$VC$K1 + x$VC$K2 - 2

eta1s <- X1  %*% t(bs[, (CLM.shift + 1) : (x$X1.d2 + CLM.shift)])
eta2s <- X2s %*% t(bs[, (x$X1.d2 + CLM.shift + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift)])

n1.s <- dim(X1)[1]
n2.s <- dim(X2s)[1]

lp1s <- matrix(nrow = n1.s, ncol = n.sim, 0)
lp2s <- matrix(nrow = n2.s, ncol = n.sim, 0)

for (i in 1 : n.sim) {
	c1s.m   <- t(matrix(nrow = x$VC$K1 - 1, ncol = n1.s, bs[i, 1 : (x$VC$K1 - 1)])) 
	c2s.m   <- t(matrix(nrow = x$VC$K2 - 1, ncol = n2.s, bs[i, x$VC$K1 : CLM.shift])) 

	eta1s.m <-   matrix(nrow = n1.s, ncol = x$VC$K1 - 1, eta1s[, i])
	eta2s.m <-   matrix(nrow = n2.s, ncol = x$VC$K2 - 1, eta2s[, i])

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	lp2s_i <- cbind(c2s.m - eta2s.m, infty)

	lp1s[, i] <- diag(lp1s_i[, y1]) # The choice of y1 selects the relevant cut point
	lp2s[, i] <- diag(lp2s_i[, y2]) # The choice of y2 selects the relevant cut point
}

### 

pk1s <- probm(lp1s, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
pk2s <- probm(lp2s, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
#eta2s <- eta.tr( X2s %*% t(bs[, (x$X1.d2 + CLM.shift2 + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift2)]), x$VC$margins[2] )


    
    




#################


if( is.null(x$X3) ){ 

if(is.null(sigma2)) sigma2 <- rep(1, dim(eta1s)[2]) 
if(is.null(nu))     nu     <- rep(1, dim(eta1s)[2])

sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }


##### cond == 0, cond == 1 and cond == 2 ##### 

p0s <- pk1s
p2s <- pk2s                    
p12s <- p0s * p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p0s

    


 
 
 
 
}# int


} # indep


list(p12 = p12, p12s = p12s, p1 = pk1, p2 = p2, p3 = NULL, CIkt = CIkt, tau = tau)


}