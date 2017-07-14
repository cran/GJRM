jc.probs4 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link){

######################################################################################################

CIp12 <- p12s <- NULL

######################################################################################################


if(type == "bivariate"){ 


######

if(!missing(newdata)){ 

X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")
if( !is.null(x$X3) ) X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")


eta1 <- X1%*%x$coef.t[1:x$X1.d2]
eta2 <- X2%*%x$coef.t[(x$X1.d2+1):(x$X1.d2+x$X2.d2)]


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
                     }


######

p1 <- probmS(eta1, x$VC$margins[1])$pr 
p2 <- probmS(eta2, x$VC$margins[2])$pr  
  
###### 


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



if(!(x$BivD %in% x$BivD2)) p12 <- BiCDF(p1, p2, x$nC, theta, 3)


if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2

############################## 


if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
if(!is.null(x$VC$mono.sm.pos)) mono.sm.pos <- x$VC$mono.sm.pos else mono.sm.pos <- c(x$VC$mono.sm.pos1, x$VC$mono.sm.pos2 + x$VC$X1.d2)  
bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos]) 

lf <- length(x$coefficients)

# try with 1 number
                       
p1s <- probmS( X1%*%t(bs[,1:x$X1.d2])                     , x$VC$margins[1])$pr 
p2s <- probmS( X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])$pr

  
if( !is.null(x$X3) ) epds <- X3%*%t(bs[,(x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
if( is.null(x$X3)  ) epds <- bs[, lf]
       
                         
est.RHOb <- teta.tr(x$VC, epds)$teta
if( is.null(x$X3) ) est.RHOb <- matrix(rep(est.RHOb, each = dim(p1s)[1]), ncol = n.sim, byrow = FALSE)

###########################   


if(x$VC$BivD %in% c("N","T")) p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, 3, test = FALSE), dim(p1s)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
                      
                        }

if(!(x$BivD %in% x$BivD2)) p12s <- BiCDF(p1s, p2s, x$nC, est.RHOb, par2 = 3, test = FALSE)
                                                                                                                        }
                                                                               
if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s


} # interv


}## biv
                        
######################################################################################################
######################################################################################################

if(type == "independence"){

if(!missing(newdata)){

X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")

eta1 <- X1%*%x$gamlss1$coef.t[1:x$X1.d2]
eta2 <- X2%*%x$gamlss2$coef.t[1:x$X2.d2]

}

if(missing(newdata)){

X1 <- x$X1
X2 <- x$X2

eta1 <- x$gamlss1$eta1
eta2 <- x$gamlss2$eta1

}


p1 <- probmS(eta1, x$VC$margins[1])$pr 
p2 <- probmS(eta2, x$VC$margins[2])$pr  
  

p12 <- p1*p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1


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

p1s <- probmS( X1%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1])$pr 
p2s <- probmS( X2%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2])$pr 

p12s <- p1s*p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p1s

                      } # intervals


} # independence


list(p12 = p12, p12s = p12s, p1 = p1, p2 = p2)


}



