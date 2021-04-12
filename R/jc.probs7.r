jc.probs7 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr){


# Prediction for ordinal-continuous case

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
pk   <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum # cdf
	pk <- pk[, y1]
eta2 <- predict(x, eq = 2, newdata = newdata, type = "response")


if( !is.null(x$X3) && x$margins[2] %in% c(cont2par,cont3par) ){

sigma2 <- esp.tr(predict(x, eq = 3, newdata = newdata, type = "response"), x$margins[2])$vrb

if(x$margins[2] %in% cont2par) eq.th <- 4
if(x$margins[2] %in% cont3par){ eq.nu <- 4; eq.th <- 5}

if(x$margins[2] %in% cont3par) nu <- enu.tr(predict(x, eq = eq.nu, newdata = newdata, type = "response"), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict(x, eq = eq.th, newdata = newdata, type = "response"))$teta

}


# new part for discrete
#if( !is.null(x$X3) && x$margins[2] %in% cont1par ) theta <- teta.tr(x$VC, 
#	predict.SemiParBIV(x, eq = 3, newdata = newdata))$teta              # DO I NEED TO INCLUDE THIS PART!?!??!




if( is.null(x$X3) ){

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}


}


##### newdata NOT included as inputs #####

if(missing(newdata)){

infty <- 1e+25

#cut <- t(matrix(nrow = x$VC$K1 - 1, ncol = x$n, x$coefficients[1 : (x$VC$K1 - 1)])) #
#eta1 <- matrix(nrow = x$n, ncol = x$VC$K1 - 1, x$eta1) ##############################

lp1 <- cbind(x$fit$lp1, infty)

lp1.sel <- lp1[, y1] # The choice of y1 selects the relevant cut point

if (y1 != x$VC$K1){
	pk <- probm(lp1.sel, x$VC$margins[1], only.pr = FALSE, bc = TRUE, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr # cumulative distribution function
} else {
	pk <- 1
}

#p1   <- x$p1
eta2 <- x$eta2

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}


##### cond == 0 and cond == 1 #####

if((cond == 0 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)) || (cond == 1 && x$margins[2] %in% c(x$VC$m2, x$VC$m3))){#*# ord - cont

p0  <- pk       
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$p2

if(!(x$BivD %in% x$BivD2)) p12 <- mm(BiCDF(p0, p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )

if(x$BivD %in% x$BivD2){

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- NA

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- mm(BiCDF(p0[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  )
if(length(theta) == 1) p12[x$teta.ind1] <- mm(BiCDF(p0[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )
                          }
                       
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- mm(BiCDF(p0[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  )
if(length(theta) == 1) p12[x$teta.ind2] <- mm(BiCDF(p0[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )
                          }                       

}

if(cond == 1 && y1 %in% seq.int(1 : x$VC$K1)) p12 <- p12 / p0
                      
}#*#


##### cond == 2 #####

if(cond == 2 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*# ord - cont

p0  <- pk     
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$p2


if(!(x$BivD %in% x$BivD2)) p12 <- copgHsCond(p0, p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2


if(x$BivD %in% x$BivD2){

p12 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- copgHsCond(p0[x$teta.ind1], p2[x$teta.ind1], theta[x$teta.ind1], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12[x$teta.ind1] <- copgHsCond(p0[x$teta.ind1], p2[x$teta.ind1], theta, dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- copgHsCond(p0[x$teta.ind2], p2[x$teta.ind2], theta[x$teta.ind2], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12[x$teta.ind2] <- copgHsCond(p0[x$teta.ind2], p2[x$teta.ind2], theta, dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }                            
                                                                           
                       }
                   
}#*#


##### Some other cases #####

#if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*# ord - discr
#
#p0   <- pk
#ppdf <- distrHsATDiscr(y2, eta2, sigma2, nu = 1, x$margins[2], x$VC$y2m)
#p2   <- ppdf$p2 
#pdf2 <- ppdf$pdf2
#
#
#if(x$BivD %in% x$BivD2){
#
#nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
#nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]
#
#C1 <- C2 <- NA
#
#if( length(x$teta1) != 0){
#
#if(length(theta) > 1){
#C1[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], p2[x$teta.ind1],          nC1, theta[x$teta.ind1], dof)
#C2[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], mm(p2[x$teta.ind1]-pdf2[x$teta.ind1]), nC1, theta[x$teta.ind1], dof)
#                     }
#
#if(length(theta) == 1){
#C1[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], p2[x$teta.ind1],          nC1, theta, dof)
#C2[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], mm(p2[x$teta.ind1]-pdf2[x$teta.ind1]), nC1, theta, dof)
#                       }
#
#                          }
#                       
#                       
#                       
#if( length(x$teta2) != 0){
#
#if(length(theta) > 1){
#C1[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], p2[x$teta.ind2],          nC2, theta[x$teta.ind2], dof)
#C2[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], mm(p2[x$teta.ind2]-pdf2[x$teta.ind2]), nC2, theta[x$teta.ind2], dof)
#                     }
#
#if(length(theta) == 1){
#C1[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], p2[x$teta.ind2],          nC2, theta, dof)
#C2[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], mm(p2[x$teta.ind2]-pdf2[x$teta.ind2]), nC2, theta, dof)
#                       }
#
#                          }                       
#
#}
#
#
#
#if(!(x$BivD %in% x$BivD2)){
#
#C1 <- BiCDF(p0, p2,          x$nC, theta, dof)
#C2 <- BiCDF(p0, mm(p2-pdf2), x$nC, theta, dof)
#
#}
#
#
#
#A <- ifelse(C1 - C2 < min.pr, min.pr, C1 - C2)
#
##if(y1 == 0){
##
##p12 <- A
##if(cond == 1) p12 <- p12/p0
##if(cond == 2) p12 <- p12/pdf2
##
##          }
##
##if(y1 == 1){
##
##p12 <- ifelse( pdf2 - A < min.pr, min.pr, pdf2 - A)
##if(cond == 1) p12 <- p12/p1
##if(cond == 2) p12 <- p12/pdf2
##
##           }
#
#if (y1 %in% seq.int(1, x$VC$K1)) {
#	p12 <- A
#	if (cond == 1) p12 / p1
#	if (cond == 2) p12 / pdf2
#}
#
#}#*#


# kendalls' tau

if(x$BivD %in% x$BivD2)   {x$SemiParFit <- x; tau <- Reg2Copost(x$SemiParFit, x$VC, theta)$tau }
if(!(x$BivD %in% x$BivD2)) tau <- ass.ms(x$VC$BivD, x$VC$nCa, theta)$tau 


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

# Adjustments for the ordinal-continuous model

if (!is.null(x$VC$K1)) {
	CLM.shift  <- x$VC$K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}


eta1s <- X1 %*% t(bs[, (CLM.shift2 + 1) : (x$X1.d2 + CLM.shift2)])
n.s   <- dim(X1)[1]


lp1s <- matrix(nrow = n.s, ncol = n.sim, 0)

for (i in 1 : n.sim) {
	c1s.m   <- t(matrix(nrow = x$VC$K1 - 1, ncol = n.s        , bs   [i, 1 : CLM.shift2])) 
	eta1s.m <-   matrix(nrow = n.s        , ncol = x$VC$K1 - 1, eta1s[ , i]              )

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	lp1s[, i] <- lp1s_i[, y1] # The choice of y1 selects the relevant cut point

}

###


pks   <- probm( lp1s, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
eta2s <- eta.tr( X2s %*% t(bs[, (x$X1.d2 + CLM.shift2 + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift2)]), x$VC$margins[2] )


#############  
# thetas
#############  

if( is.null(x$X3) ) epds <- bs[,length(x$coefficients)]
  
if( !is.null(x$X3) ){ 


  if(x$VC$margins[2] %in% cont1par){   
  
  
if(!missing(newdata)){ X3s <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}    
  
                epds <- X3s%*%t(bs[,((x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)) + CLM.shift2])
  
                                   }


  if(x$VC$margins[2] %in% cont2par){   
  
  
if(!missing(newdata)){ X4s <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
                epds <- X4s%*%t(bs[,((x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)) + CLM.shift2])
  
                                   }
         
         
  if(x$VC$margins[2] %in% cont3par){
  
if(!missing(newdata)){ X5s <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X5s <- x$X5s else X5s <- x$X5}    
  
                epds <- X5s%*%t(bs[,((x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)) + CLM.shift2])
                
                                   }
  
  
  
  }

est.RHOb <- teta.tr(x$VC, epds)$teta
   

#############  
# sigmas
#############  

if( x$VC$margins[2] %in% cont2par ){  

      if( is.null(x$X3) )   sigma2.star <- bs[, CLM.shift2 + x$X1.d2 + x$X2.d2 + 1] 
      if( !is.null(x$X3) ) {
      
if(!missing(newdata)){ X3s <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
                
            sigma2.star <- X3s %*% t(bs[,((x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)) + CLM.shift2]) 
  
                           }
  
sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
}    
    

#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs[, CLM.shift2 + x$X1.d2 + x$X2.d2 + 2] # t(as.matrix(bs[,  x$X1.d2 + x$X2.d2 + 2]))
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")}                    
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
              nu.st <- X4s %*% t(bs[, ((x$X1.d2 + x$X2.d2 + x$X3.d2 + 1) : (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)) + CLM.shift2]) 
  
                      }
  
 nu <- enu.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


if( is.null(x$X3) ){

est.RHOb <- matrix(rep(est.RHOb, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }


##### cond == 0 and cond == 1 #####

if((cond == 0 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)) || (cond == 1 && x$margins[2] %in% c(x$VC$m2, x$VC$m3))){#*# ord - cont

p0s  <- pks         
p2s  <- matrix( distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$p2 , dim(p0s)[1], n.sim)

if(x$VC$BivD %in% c("N","T")) p12s <- matrix(mm(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  ), dim(pks)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- mm(BiCDF(p0s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,]), min.pr = min.pr, max.pr = max.pr  )                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- mm(BiCDF(p0s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,]), min.pr = min.pr, max.pr = max.pr  )
                      
                        }

if(!(x$BivD %in% x$BivD2)) {
	p12s <- mm(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  )    
	p12s <- matrix(nrow = n.s, ncol = n.sim, p12s) # This part is not in jc.probs2
}

}

if(cond == 1 && y1 %in% seq.int(1 : x$VC$K1)) p12s <- p12s / p0s
 
}#*# 


##### cond == 2 #####
 
if(cond == 2 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*# ord - cont

p0s  <- pks         
p2s  <- matrix( distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$p2 , dim(p0s)[1], n.sim)


if(!(x$BivD %in% x$BivD2)) p12s <- matrix( copgHsCond(p0s, p2s, est.RHOb, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2 , dim(p0s)[1], n.sim) 
         if(x$BivD == "T") p12s <- matrix(p12s, dim(pks)[1], n.sim)


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(pks)[1])
 
if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- copgHsCond(p0s[x$teta.ind1,], p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2                                               
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- copgHsCond(p0s[x$teta.ind2,], p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                                                                                                     
                       }

}


##### Some other cases #####

#if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#
#
#p0s   <- 1 - p1s
#ppdf  <- distrHsATDiscr(y2, eta2s, sigma2, nu = 1, x$margins[2], x$VC$y2m) # not sure about y2m
#p2s   <- ppdf$p2 
#pdf2s <- ppdf$pdf2
#
#if(x$VC$BivD %in% c("N","T")) C1s <- matrix(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), dim(p1s)[1], n.sim) else{
#
#
#if(x$BivD %in% x$BivD2){
#
#C1s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])
#
#if( length(x$teta1) != 0) C1s[x$teta.ind1,] <- BiCDF(p0s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
#if( length(x$teta2) != 0) C1s[x$teta.ind2,] <- BiCDF(p0s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
#                      
#                        }
#
#if(!(x$BivD %in% x$BivD2)) C1s <- BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE)
#
#
#
#}
#
#
#
#if(x$VC$BivD %in% c("N","T")) C2s <- matrix(BiCDF(p0s, mm(p2s-pdf2s), x$nC, est.RHOb, dof, test = FALSE), dim(p1s)[1], n.sim) else{
#
#
#if(x$BivD %in% x$BivD2){
#
#C2s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])
#
#if( length(x$teta1) != 0) C2s[x$teta.ind1,] <- BiCDF(p0s[x$teta.ind1,], mm(p2s[x$teta.ind1,]-pdf2s[x$teta.ind1,]), nC1,  est.RHOb[x$teta.ind1,])                  
#if( length(x$teta2) != 0) C2s[x$teta.ind2,] <- BiCDF(p0s[x$teta.ind2,], mm(p2s[x$teta.ind2,]-pdf2s[x$teta.ind2,]), nC2, -est.RHOb[x$teta.ind2,])
#                      
#                        }
#
#if(!(x$BivD %in% x$BivD2)) C2s <- BiCDF(p0s, mm(p2s-pdf2s), x$nC, est.RHOb, dof, test = FALSE)
#
#
#}
#
#
#
#
#As <- ifelse(C1s - C2s < min.pr, min.pr, C1s - C2s)
#
#
#
#
#if(y1 == 0){
#p12s <- As
#if(cond == 1) p12s <- p12s/p0s
#if(cond == 2) p12s <- p12s/pdf2s
#          }
#
#if(y1 == 1){
#p12s <- ifelse( pdf2s - As < min.pr, min.pr, pdf2s - As)
#if(cond == 1) p12s <- p12s/p1s
#if(cond == 2) p12s <- p12s/pdf2s
#           }
#       
#}#*# 
 
   
   
   
   
   
   
   
   

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

nu   <- sigma2 <- NA
pk   <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum # cdf
	pk <- pk[, y1]
eta2 <- predict(x, eq = 2, newdata = newdata, type = "response")


if( !(x$VC$margins[2] %in% cont1par) ){

if( !is.null(x$X3) ){

sigma2 <- esp.tr(predict(x, eq = 3, newdata = newdata, type = "response"), x$margins[2])$vrb
if(x$margins[2] %in% cont3par) nu <- enu.tr(predict(x, eq = 4, newdata = newdata, type = "response"), x$margins[2])$vrb

                    }

if( is.null(x$X3) ){

sigma2 <- x$sigma2
nu     <- x$nu 

                   }


                                 }


}


##### newdata NOT included as inputs #####

if(missing(newdata)){

infty <- 1e+25

#cut <- t(matrix(nrow = x$VC$K1 - 1, ncol = x$n, x$coefficients[1 : (x$VC$K1 - 1)])) #
#eta1 <- matrix(nrow = x$n, ncol = x$VC$K1 - 1, x$eta1) ##############################

lp1 <- cbind(x$fit$lp1, infty)

lp1.sel <- lp1[, y1] # The choice of y1 selects the relevant cut point

if (y1 != x$VC$K1){
	pk <- probm(lp1.sel, x$VC$margins[1], only.pr = FALSE, bc = TRUE, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr # cumulative distribution function
} else {
	pk <- 1
}

#p1   <- x$p1
eta2 <- x$eta2

sigma2 <- x$sigma2
nu     <- x$nu 

}


##### cond == 0, cond == 1 and cond == 2 #####

if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*# ord - cont

#if(y1 == 0 || y1 == 1){  

p0  <- pk                   
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$p2
p12 <- p0 * p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p0

#                      }

#if(y1 == 1){                                          
#
#p12 <- p1*p2 
#
#if(cond == 1) p12 <- p2
#if(cond == 2) p12 <- p1
#
#           }
                     
}  #*#         


##### Some other cases #####

#if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#
#
#p0   <- 1 - p1
#pdf2 <- distrHsATDiscr(y2, eta2, sigma2, nu = 1, x$margins[2], x$VC$y2m)$pdf2  
#
#if(y1 == 0){
#p12 <- p0*pdf2
#if(cond == 1) p12 <- pdf2
#if(cond == 2) p12 <- p0
#          }
#
#if(y1 == 1){
#p12 <- p1*pdf2
#if(cond == 1) p12 <- pdf2
#if(cond == 2) p12 <- p1
#           }
#       
#}#*#  


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


#bs1 <- rMVN(n.sim, mean = x$gam1$coefficients, sigma=x$gam1$Vp)
#bs2 <- rMVN(n.sim, mean = x$gamlss$coefficients, sigma=x$gamlss$Vb)

bs <- rMVN(n.sim, mean = coefficients.s, sigma = x$Vb)  


# ... and then transformed back

cut.sim.ti <- bs[, 1 : (x$VC$K1 - 1)]
	cut.sim <- matrix(nrow = n.sim, ncol = x$VC$K1 - 1, 0)	
	cut.sim[, 1] <- cut.sim.ti[, 1] ; for (i in 2 : (x$VC$K1 - 1)) cut.sim[, i] <- cut.sim[, i - 1] + cut.sim.ti[, i]^2

bs[, 1 : (x$VC$K1 - 1)] <- cut.sim


#############  
# etas
#############  

##### newdata included as inputs #####

if(!missing(newdata)){ pred.1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
				X1 <- pred.1 ; colnames(X1) <- colnames(x$X1)
                       X2s <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix") }

                       
##### newdata NOT included as inputs #####

if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       } 
  

###

# Adjustments for the ordinal-continuous model

if (!is.null(x$VC$K1)) {
	CLM.shift  <- x$VC$K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}


eta1s <- X1 %*% t(bs[, (CLM.shift2 + 1) : (x$X1.d2 + CLM.shift2)])
n.s   <- dim(X1)[1]


lp1s <- matrix(nrow = n.s, ncol = n.sim, 0)

for (i in 1 : n.sim) {
	c1s.m   <- t(matrix(nrow = x$VC$K1 - 1, ncol = n.s        , bs   [i, 1 : CLM.shift2])) 
	eta1s.m <-   matrix(nrow = n.s        , ncol = x$VC$K1 - 1, eta1s[ , i]              )

	lp1s_i <- cbind(c1s.m - eta1s.m, infty)
	lp1s[, i] <- lp1s_i[, y1] # The choice of y1 selects the relevant cut point

}

###


pks   <- probm( lp1s, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr   
eta2s <- eta.tr( X2s %*% t(bs[, (x$X1.d2 + CLM.shift2 + 1) : (x$X1.d2 + x$X2.d2 + CLM.shift2)]) , x$VC$margins[2])


#############  
# sigmas   
#############  

if( x$VC$margins[2] %in% cont2par ){ 

      if( is.null(x$X3) )   sigma2.star <- bs[, CLM.shift2 + x$X1.d2 + x$X2.d2 + 1] 
      
      if( !is.null(x$X3) ){
      
if(!missing(newdata)){ X3s <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")}                     
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
      
                      sigma2.star <- X3s %*% t(bs[, ((x$X1.d2 + x$X2.d2 + 1) : (x$X1.d2 + x$X2.d2 + x$X3.d2)) + CLM.shift2]) 

                          }

sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
}    
    

#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs[, CLM.shift2 + x$X1.d2 + x$X2.d2 + 2] # t(as.matrix(bs2[,  x$X2.d2 + 2]))
  
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")}                      
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}     
  
             nu.st <- X4s %*% t(bs[, ((x$X1.d2 + x$X2.d2 + x$X3.d2 + 1) : (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)) + CLM.shift2]) 
                      }
                      
 nu <- enu.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


if( is.null(x$X3) ){

sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }


##### cond == 0, cond == 1 and cond == 2 #####

if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*# ord - cont

#if(y1 == 0 || y1 == 1){  

p0s  <- pks                   
p2s  <- matrix( distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$p2 , dim(p0s)[1], n.sim)
p12s <- p0s * p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p0s

#                      }

#if(y1 == 1){                                          
#
#p12s <- p1s*p2s
#
#if(cond == 1) p12s <- p2s
#if(cond == 2) p12s <- p1s
#
#           }
                     
}  #*#     


##### Some other cases #####

#if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#
#
#p0s   <- 1 - p1s
#pdf2s <- distrHsATDiscr(y2, eta2s, sigma2, nu = 1, x$margins[2], x$VC$y2m)$pdf2 # not sure about y2m
#
#if(y1 == 0){
#p12s <- p0s*pdf2s
#if(cond == 1) p12s <- pdf2s
#if(cond == 2) p12s <- p0s
#          }
#
#if(y1 == 1){
#p12s <- p1s*pdf2s
#if(cond == 1) p12s <- pdf2s
#if(cond == 2) p12s <- p1s
#           }
#       
#}#*#  
 
 
 
}# int


} # indep


list(p12 = p12, p12s = p12s, p1 = pk, p2 = p2, p3 = NULL, CIkt = CIkt, tau = tau)


}