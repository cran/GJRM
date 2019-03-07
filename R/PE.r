PE <- function(x1, idx, n.sim = 100, prob.lev = 0.05, 
                  hd.plot = FALSE, 
                  main = "Histogram and Kernel Density of Simulated Average Effects", 
                  xlab = "Simulated Average Effects", ...){


etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- C.11 <- C.10 <- C.11s <- C.10s <- v0s <- v1s <- ATs <- AUXs <- NULL
nu  <- 1 
nus <- rep(1, n.sim)



if(!(x1$margins[1] %in% x1$bl) ) stop("Error in first margin value. It should be one of:\nprobit, logit, cloglog.")
if(!(x1$margins[2] %in% x1$bl) ) stop("Error in second margin value. It should be one of:\nprobit, logit, cloglog.")



#contd <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","DAGUM","SM","GP","BE")

#if( !(x2$margins[2] %in% contd) ) stop("This function currently works only for continuous margins.")


if( missing(idx) ) stop("You must provide a value for idx.")  


#index1 <- 1:length(x1$y2)
#index1 <- index1[x1$y2==1]

#if( idx < 1 || idx > length(index1) ) stop("idx out of range.")  


#t1 <- x1$y1[x1$y2==1]
#t2 <- x2$y1

#if(all.equal(t1,t2) == FALSE) stop("The first equations of the two models must have the same (binary) response variable.") 


m2       <- x1$VC$m2 
m3       <- x1$VC$m3 
bin.link <- x1$VC$bl  

end <- 0
epsilon <- 0.0000001 # 0.9999999 sqrt(.Machine$double.eps)
max.p   <- 0.9999999
est.ATb <- NA
indD <- list()


#if(x2$margins[2] == "DAGUM") { if( min(sqrt(x2$sigma2[idx])) <= 1) stop("sigma parameter has value(s) smaller than 1, hence the mean is indeterminate.")}
#if(x2$margins[2] == "SM"   ) { if( min(sqrt(x2$sigma2[idx])*x2$nu[idx]) <= 1) stop("sigma*nu has value(s) smaller than 1, hence the mean is indeterminate.")}


########################################################
# model x1 - binary binary
########################################################

ind1  <- 1:x1$X1.d2 
ind2  <- x1$X1.d2 + (1:x1$X2.d2)

#X1 <- as.matrix((x1$X1[index1,])[idx,])
#X2 <- as.matrix((x1$X2[index1,])[idx,])
#eta1 <- (x1$eta1[index1])[idx]
#eta2 <- (x1$eta2[index1])[idx]#


X1 <- as.matrix(x1$X1[idx,])
X2 <- as.matrix(x1$X2[idx,])

eta1 <- t(X1)%*%x1$coefficients[ind1] 
eta2 <- t(X2)%*%x1$coefficients[ind2]  


p.1 <- probm(eta1, x1$margins[1])$pr 
p.2 <- probm(eta2, x1$margins[1])$pr

if( is.null(x1$X3)  )   ass.p <- x1$theta    

if( !is.null(x1$X3) ) { X3    <- as.matrix( x1$X3[idx,] )
                        etd   <- t(X3)%*%x1$coefficients[(x1$X1.d2+x1$X2.d2+1):(x1$X1.d2+x1$X2.d2+x1$X3.d2)] 
                        ass.p <- teta.tr(x1$VC, etd)$teta
                       }                                                        #  ass.p <- (x1$theta[index1])[idx]  



AUX   <- pmax( BiCDF(p.1, p.2, x1$nC, ass.p) , epsilon ) # this is p11
C.11  <- AUX / p.1
C.10  <- pmax( p.2 - AUX, epsilon )   / (1 - p.1)   # this is actually C01, this has been corrected
              
              
######
# CIs
######
              
bs <- rMVN(n.sim, mean = x1$coefficients, sigma=x1$Vb)

eta1s <- t(X1)%*%t(bs[,ind1]) 
eta2s <- t(X2)%*%t(bs[,ind2])  

p.1s <- probm(eta1s, x1$margins[1])$pr 
p.2s <- probm(eta2s, x1$margins[1])$pr

if( !is.null(x1$X3) ) etds <- t(X3)%*%t(bs[,(x1$X1.d2+x1$X2.d2+1):(x1$X1.d2+x1$X2.d2+x1$X3.d2)])
if(  is.null(x1$X3) ) etds <- bs[, length(x1$coefficients)]
   

resT    <- teta.tr(x1$VC, etds)
ass.ps  <- resT$teta 



if(x1$BivD == "N") {

	for(i in 1:n.sim){ 
	
	         AUXs[i]  <- pmax( BiCDF(p.1s[i], p.2s[i], x1$nC, ass.ps[i], test = FALSE), epsilon )
		 C.11s[i] <- AUXs[i]  / p.1s[i]
		 C.10s[i] <- ( p.2s[i] - AUXs[i] ) / (1 - p.1s[i])  
	                  }                 
}

if(x1$BivD != "N") {


 AUXs  <- pmax( BiCDF(p.1s, p.2s, x1$nC, ass.ps, test = FALSE), epsilon )
 C.11s <- AUXs / p.1s 
 C.10s <- ( p.2s - AUXs ) / (1 - p.1s)
             
}





########################################################
# model x2
########################################################

# idx is for a specific individual


#ind1  <- 1:x2$X1.d2 
#ind2  <- x2$X1.d2+(1:x2$X2.d2)
#
#X1 <- as.matrix(x2$X1[idx,])
#X2 <- as.matrix(x2$X2[idx,])
#
#
#eta1 <- x2$eta1[idx]
#eta2 <- x2$eta2[idx] 
#
#
#ass.p  <- x2$theta
#sigma2 <- x2$sigma2
#nu     <- x2$nu
#
#
#
#if( x2$margins[2] %in% m2){
#     if( !is.null(x2$X3) && !is.null(x2$X4) ) { ass.p <- ass.p[idx]; sigma2 <- sigma2[idx] }
#}
#
#
#
#if( x2$margins[2] %in% m3){
#    if( !is.null(x2$X3) && !is.null(x2$X4) && !is.null(x2$X5)) { ass.p <- ass.p[idx]; sigma2 <- sigma2[idx]; nu <- nu[idx] }
#}
#
#########################################################
#
#
#
#
#if( x2$margins[2] %in% c("N","GU","rGU","LO") )                             { lil <- -Inf;    uil <- Inf}
#if( x2$margins[2] %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","GP")  )  { lil <- epsilon; uil <- Inf}
#if( x2$margins[2] %in% c("BE")  )                                           { lil <- epsilon; uil <- max.p}
#
#
#j <- 1; inds <- posi <- NULL 
#tol <- 1e-4
#
#lo <- length(x2$y2)#*2 #1000
##rlo <- (max(x2$y2) - min(x2$y2))/lo
##if(rlo > 1) lo <- round(lo*rlo) 
#
#
#if( x2$margins[2] %in% c("N","N2","GU","rGU","LO","LN") )               seq.y <- seq(min(x2$y2) - ((max(x2$y2) - min(x2$y2))/2), max(x2$y2) + ((max(x2$y2) - min(x2$y2))/2), length.out = lo)
#if( x2$margins[2] %in% c("WEI","iG","GA","DAGUM","SM","FISK","GP")  )  seq.y <- seq(1e-12, max(x2$y2) + ((max(x2$y2) - min(x2$y2))/2), length.out = lo)
#if( x2$margins[2] %in% c("BE")  )                                       seq.y <- seq(1e-12, 0.9999999,  length.out = lo)      
#
#
#seq.y <- seq(1e-12, max(x2$y2), length.out = lo)
#iresS <- NA
#
#for(i in 1:lo){ 
#
#  ires  <- distrHsAT(seq.y[i], eta2, sigma2, nu, x2$margins[2])
#  ires0 <- ires$p2
#  ires1 <- ires$pdf2
#  p.1 <- 1 - probm(eta1, x2$margins[1])$pr 
#  p.1 <- pmax(p.1, epsilon ) 
#  p.1 <- ifelse(p.1 > max.p, max.p, p.1)  
#  ires2 <- copgHsAT(p1 = p.1, p2 = ires0, teta = ass.p, BivD = x2$BivD)$c.copula.be2
#  
#  ires <- iresS[i] <- ires2*ires1/p.1  
#  
#  
#  if(min(ires) == max(ires) && min(ires) < tol) inds[i] <- TRUE else inds[i] <- FALSE 
#  if(i > 1 && inds[i] != inds[i-1]) { if(j==1){ posi[j] <- i-1; j <- j + 1} else posi[2] <- i}  
#
#}
#
#if((is.null(posi) || length(posi) < 2) && x2$margins[2] %in% c("N","N2","GU","rGU","LO","LN")) stop("Increase the tolerance value or try a different range.")
#
#
#
#if((is.null(posi) || length(posi) < 2) && x2$margins[2] %in% c("WEI","iG","GA","DAGUM","SM","FISK","GP") ){
#
#
#posi <- NULL; j <- 1
#
#  for(i in 1:lo){ 
#
#  ires  <- distrHsAT(seq.y[i], eta2, sigma2, nu, x2$margins[2])
#  ires0 <- ires$p2
#  ires1 <- ires$pdf2
#  p.1 <- 1 - probm(eta1, x2$margins[1])$pr 
#  p.1 <- pmax(p.1, epsilon ) 
#  p.1 <- ifelse(p.1 > max.p, max.p, p.1)  
#  ires2 <- copgHsAT(p1 = p.1, p2 = ires0, teta = ass.p, BivD = x2$BivD)$c.copula.be2
#  
#  ires <- iresS[i] <- ires2*ires1/p.1 
#  
#    if(min(ires) == max(ires) && min(ires) < tol) inds[i] <- TRUE else inds[i] <- FALSE 
#    
#    if(i > 2 && inds[i] != inds[i-1]) { if(j == 1) posi[2] <- i; j <- j + 1}
#    
#    if(!is.null(posi[2])) posi[1] <- 1 
#
#    }
#
#  if(is.null(posi) || length(posi) < 2) stop("Increase the tolerance value or try a different range.")
#
#
#}
#
#
#
#lil <- seq.y[posi[1]]  # to check again
#uil <- seq.y[posi[2]]  # max(x2$y2)*100 # CAREFUL HERE!! seq.y[posi[2]]
#
#
#
#########################################################
#########################################################
#
#ConExp0 <- function(y2, eta2, sigma2, nu, margin2, p.1, ass.p, BivD){   # to check 
#
#	pp0 <- distrHsAT(y2, eta2, sigma2, nu, margin2) 
#	p2.0   <- pp0$p2
#	pdf2.0 <- pp0$pdf2 
#	dc0 <- copgHsAT(p1 = p.1, p2=p2.0, teta=ass.p, BivD=BivD)$c.copula.be2	
#	cond0 <- y2*dc0*pdf2.0/p.1
#	cond0
#                                                                   }
#
#ConExp1 <- function(y2, eta2, sigma2, nu, margin2, p.1, ass.p, BivD){   
#
#	pp1 <- distrHsAT(y2, eta2, sigma2, nu, margin2)
#	p2.1   <- pp1$p2
#	pdf2.1 <- pp1$pdf2 	
#	dc1 <- copgHsAT(p1 = p.1, p2=p2.1, teta=ass.p, BivD=BivD)$c.copula.be2	
#	cond1 <- y2 * ( (1 - dc1)*pdf2.1/(1-p.1))
#	cond1
#                                                                   }
#                                                                   
#
#integr0 <- function(eta2, sigma2, nu, margin2, p.1, ass.p, BivD, lil, uil){
#
#  distrExIntegrate(ConExp0, lower=lil, upper=uil, eta2=eta2, 
#            sigma2=sigma2, nu = nu, margin2=margin2, p.1=p.1, 
#            ass.p=ass.p, BivD=BivD)    #$value
#                         
#                       }
#
#integr1 <- function(eta2, sigma2, nu, margin2, p.1, ass.p, BivD, lil, uil){
#
#  distrExIntegrate(ConExp1, lower=lil, upper=uil, eta2=eta2,  
#            sigma2=sigma2, nu= nu, margin2=margin2, p.1 = p.1, 
#            ass.p=ass.p, BivD=BivD)    #$value
#                         
#                       }
#
#v.integr0 <- Vectorize(integr0)  
#v.integr1 <- Vectorize(integr1) 
#
#
########################################################
#
#
#p.1 <- 1 - probm(eta1, x2$margins[1])$pr 
#
#p.1 <- pmax(p.1, epsilon ) 
#p.1 <- ifelse(p.1 > max.p, max.p, p.1)
#
#
#
#v0 <- v.integr0(eta2 = eta2, sigma2=sigma2, nu = nu, margin2=x2$margins[2], p.1=p.1, ass.p=ass.p, BivD=x2$BivD, lil=lil, uil=uil)            
#v1 <- v.integr1(eta2=eta2, sigma2=sigma2, nu = nu, margin2=x2$margins[2], p.1=p.1, ass.p=ass.p, BivD=x2$BivD, lil=lil, uil=uil)
#
#####################
## Effect of interest
#####################
#
#
#AT <- C.11*v1 - C.10*v0  


AT <- C.11 - C.10  

   
   
##################

#bs <- rMVN(n.sim, mean = x2$coefficients, sigma=x2$Vb)
#
#eta1s <- t(X1)%*%t(bs[,ind1])
#eta2s <- t(X2)%*%t(bs[,ind2])  
#
#p.1s <- pmax(  1 - probm(eta1s, x2$margins[1])$pr , epsilon )
#p.1s <- ifelse(p.1s > max.p, max.p, p.1s)
#
#
#
#
#   
#   if(x2$margins[2] %in% m2 ){
#   
#   if( !is.null(x2$X3) ) sigma2.st <- t(x2$X3[idx,])%*%t(bs[,(x2$X1.d2+x2$X2.d2+1):(x2$X1.d2+x2$X2.d2+x2$X3.d2)]) 
#   if(  is.null(x2$X3) ) sigma2.st <- bs[, length(x2$coefficients) - 1]
#   
#   if( !is.null(x2$X4) ) etds <- t(x2$X4[idx,])%*%t(bs[,(x2$X1.d2+x2$X2.d2+x2$X3.d2 + 1):(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2)])
#   if(  is.null(x2$X4) ) etds <- bs[, length(x2$coefficients)]  
#   
#   sigma2s <- esp.tr(sigma2.st, x2$margins[2])$vrb
#   
#   
#  
#  }
#  
#   if(x2$margins[2] %in% m3 ){
#   
#   if( !is.null(x2$X3) ) sigma2.st <- t(x2$X3[idx,])%*%t(bs[,(x2$X1.d2+x2$X2.d2+1):(x2$X1.d2+x2$X2.d2+x2$X3.d2)]) 
#   if(  is.null(x2$X3) ) sigma2.st <- bs[,length(x2$coefficients) - 2]
#   
#   if( !is.null(x2$X4) ) nu.st <- t(x2$X4[idx,])%*%t(bs[,(x2$X1.d2+x2$X2.d2+x2$X3.d2+1):(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2)]) 
#   if(  is.null(x2$X4) ) nu.st <- bs[, length(x2$coefficients) - 1]   
#   
#   if( !is.null(x2$X5) ) etds <- t(x2$X5[idx, ])%*%t(bs[,(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2 + 1):(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2+x2$X5.d2)])
#   if(  is.null(x2$X5) ) etds <- bs[,length(x2$coefficients)]  
#   
#   sigma2s <- esp.tr(sigma2.st, x2$margins[2])$vrb
#   
#   nus <- esp.tr(nu.st, x2$margins[2])$vrb
#      
#      
#  }  
# 
# 
# 
#
#
#resT    <- teta.tr(x2$VC, etds)
#ass.ps  <- resT$teta 
#
#
#
#
#
#for(i in 1:n.sim){ 
#
#v0s[i] <- v.integr0(eta2 = eta2s[i], sigma2=sigma2s[i], nu=nus[i], margin2=x2$margins[2], p.1 = p.1s[i], ass.p=ass.ps[i], BivD=x2$BivD, lil=lil, uil=uil)
#v1s[i] <- v.integr1(eta2 = eta2s[i], sigma2=sigma2s[i], nu=nus[i], margin2=x2$margins[2], p.1 = p.1s[i], ass.p=ass.ps[i], BivD=x2$BivD, lil=lil, uil=uil)
#
#
#                 }
#
#
#
#
#
#for(i in 1:n.sim) ATs[i] <- C.11s[i]*v1s[i] - C.10s[i]*v0s[i]


ATs <- C.11s - C.10s



 if(hd.plot == TRUE){
  
  hist(ATs, freq = FALSE, main = main, xlab = xlab, 
       ylim = c(0, max(density(ATs)$y, hist(ATs, plot = FALSE)$density)), ...)
  lines(density(ATs))
 
 }




CIs <- as.numeric(quantile(ATs, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

res <- c(CIs[1], AT, CIs[2])

out <- list(res = res, prob.lev = prob.lev) # , lims = c(lil, uil) ) 

class(out) <- "PE"

out

}

