jc.probs2 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = FALSE){

nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- NULL
CIkt <- tau <- theta <- CItheta <- NULL
CLM.shift <- CLM.shift2 <- 0
dof <- x$dof

p12s  <- matrix(0, 1, 2) 



nC   <- x$nC
BivD <- x$BivD


if(type == "joint"){


if(!missing(newdata)){

nu <- sigma2 <- NA

p1   <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "response")
eta2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata)



###############
if(!is.null(x$VC$K1)){
 
 p1 <- seq(0.001, 0.999, length.out = dim(newdata)[1])  # this is a trick/cheat and will have to be fixed properly

}
###############

if( !is.null(x$X3) && x$margins[2] %in% c(cont2par,cont3par) ){

sigma2 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata), x$margins[2])$vrb

if(x$margins[2] %in% cont2par) eq.th <- 4
if(x$margins[2] %in% cont3par){ eq.nu <- 4; eq.th <- 5}

if(x$margins[2] %in% cont3par) nu <- enu.tr(predict.SemiParBIV(x, eq = eq.nu, newdata = newdata), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = eq.th, newdata = newdata))$teta

}


# new part for discrete
if( !is.null(x$X3) && x$margins[2] %in% cont1par ) theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = 3, newdata = newdata))$teta




if( is.null(x$X3) ){

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}


}



if(missing(newdata)){

p1   <- x$p1
eta2 <- x$eta2

sigma2 <- x$sigma2
nu     <- x$nu 
theta  <- x$theta 

}




if(x$BivD %in% x$BivD2){ 

if(theta >= 0){ nC <- x$fit$nC1; BivD <- x$fit$Cop1 } else{ nC <- x$fit$nC2; BivD <- x$fit$Cop2} 

}





if((cond == 0 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)) || (cond == 1 && x$margins[2] %in% c(x$VC$m2, x$VC$m3))){#*# bin - cont



p0  <- 1 - p1         
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2

p12 <- mm(BiCDF(p0, p2, nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )




if(cond == 1 && y1 == 0) p12 <- p12/p0

                      

if(y1 == 1){                         ### joint prob of y1=1 and y2=x ###

p12  <- p2 - p12

if(cond == 1) p12 <- p12/p1

           }
                      
}#*#





if(cond == 2 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*# bin - cont


p0  <- 1 - p1         
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2


p12 <- copgHsCond(p0, p2, theta, dof = dof, BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2




if(y1 == 1) p12  <- 1 - p12 ### y1=1 | y2=x ###

                   
}#*#









if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){#*# bin - discr

p0   <- 1 - p1
#y2rep <- rep(y2, length(eta2))
ppdf <- distrHsATDiscr(y2, eta2, sigma2, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)
p2   <- ppdf$p2 
pdf2 <- ppdf$pdf2


C1 <- mm(BiCDF(p0, p2,          nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )
C2 <- mm(BiCDF(p0, mm(p2-pdf2, min.pr, max.pr), nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )


A <- ifelse(C1 - C2 < min.pr, min.pr, C1 - C2)

if(y1 == 0){

p12 <- A
if(cond == 1) p12 <- p12/p0
if(cond == 2) p12 <- p12/pdf2

          }

if(y1 == 1){

p12 <- ifelse( pdf2 - A < min.pr, min.pr, pdf2 - A)
if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/pdf2

           }
 
 
 
}#*#








# kendalls' tau

tau   <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = tau.res, t.theta = theta)$tau
theta <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = FALSE,   t.theta = theta)$theta
           
           
if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

#############  
# etas
############# 

if(!missing(newdata)){ X1  <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") 


###############
if(!is.null(x$VC$K1)){
 
 X1 <- matrix(1, dim(newdata)[1], x$X1.d2) # this is a trick/cheat and will have to be fixed properly

}
###############




                       X2s <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       }  


p1s   <- probm( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr 
eta2s <- eta.tr( X2s%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])

#############  
# thetas
#############  

if( is.null(x$X3) ) epds <- bs[,length(x$coefficients)]
  
if( !is.null(x$X3) ){ 


  if(x$VC$margins[2] %in% cont1par){   
  
  
if(!missing(newdata)){ X3s <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}    
  
                epds <- X3s%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])
  
                                   }


  if(x$VC$margins[2] %in% cont2par){   
  
  
  

if (!is.null(x$VC$K1)) {
    
        K1  <- x$VC$K1
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- CLM.shift2 <- 0
}    
  
  
if(!missing(newdata)){ X4s <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
                epds <- X4s%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1+ CLM.shift2):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+ CLM.shift2)])
  
                                   }
         
         
  if(x$VC$margins[2] %in% cont3par){
  
if(!missing(newdata)){ X5s <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X5s <- x$X5s else X5s <- x$X5}    
  
                epds <- X5s%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
                
                                   }
  
  
  
  }

est.RHOb <- teta.tr(x$VC, epds)$teta
   
#############  
# sigmas
#############  

if( x$VC$margins[2] %in% c(cont2par, cont3par)  ){  

      if( is.null(x$X3) )   sigma2.star <- bs[, x$X1.d2 + x$X2.d2 + 1] 
      if( !is.null(x$X3) ) {
      
if(!missing(newdata)){ X3s <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
                
            sigma2.star <- X3s%*%t(bs[,(x$X1.d2+x$X2.d2+1+ CLM.shift2):(x$X1.d2+x$X2.d2+x$X3.d2+ CLM.shift2)]) 
  
                           }
  
sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
}    
    
#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs[, x$X1.d2 + x$X2.d2 + 2] 
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}                    
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
              nu.st <- X4s%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
                      }
  
 nu <- enu.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


#if( is.null(x$X3) ){
#
#est.RHOb <- matrix(rep(est.RHOb, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
#if( x$VC$margins[2] %in% c(cont2par, cont3par)  ) sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
#if( x$VC$margins[2] %in% cont3par )               nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
#
#                   }




if((cond == 0 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)) || (cond == 1 && x$margins[2] %in% c(x$VC$m2, x$VC$m3))){#*#

p0s  <- 1 - p1s         
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2 

if(x$VC$BivD %in% c("N","T")) p12s <- mm(BiCDF(p0s, p2s, nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  ) else{

p12s <- mm(BiCDF(p0s, p2s, nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  ) 


}

if(cond == 1 && y1 == 0) p12s <- p12s/p0s


if(y1 == 1){                         
  
p12s  <- p2s - p12s

if(cond == 1) p12s <- p12s/p1s

           }
 
}#*# 
 
 
 
if(cond == 2 && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*# bin - cont

p0s  <- 1 - p1s         
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2


p12s <- copgHsCond(p0s, p2s, est.RHOb, dof = dof, BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2

if(y1 == 1) p12s <- 1 - p12s ### y1=1 | y2=x ###


}




 
 
 
 
if(x$margins[2] %in% c(x$VC$m1d,x$VC$m2d)){#*#

p0s   <- 1 - p1s 


# y2rep <- rep(y2, length(eta2s))


if( x$margins[2] != c("BB") ) y2rep <- rep(y2, length(eta2s))
if( x$margins[2] == c("BB") ) y2rep <- matrix(rep(y2, each = length(eta2s)), ncol = 2, byrow = FALSE)


ppdf  <- distrHsATDiscr(y2rep, eta2s, sigma2, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2) # not sure about y2m
p2s   <- ppdf$p2 
pdf2s <- ppdf$pdf2

if(x$VC$BivD %in% c("N","T")) C1s <- mm(BiCDF(p0s, p2s, nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  ) else{


C1s <- mm(BiCDF(p0s, p2s, nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  )



}



if(x$VC$BivD %in% c("N","T")) C2s <- mm(BiCDF(p0s, mm(p2s-pdf2s, min.pr, max.pr), nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  ) else{



C2s <-  mm(BiCDF(p0s, mm(p2s-pdf2s, min.pr, max.pr), nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  )


}




As <- ifelse(C1s - C2s < min.pr, min.pr, C1s - C2s)




if(y1 == 0){
p12s <- As
if(cond == 1) p12s <- p12s/p0s
if(cond == 2) p12s <- p12s/pdf2s
          }

if(y1 == 1){
p12s <- ifelse( pdf2s - As < min.pr, min.pr, pdf2s - As)
if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/pdf2s
           }
       
}#*# 
 
   
 

ass.msR <- theta2tau(x$VC$BivD, x$VC$nCa, est.RHOb, tau.res = tau.res, t.theta = theta)
if(tau.res == TRUE) taus <- ass.msR$tau
thetas <- ass.msR$theta

                               
                                   
                                                  
                                                  
if(tau.res == TRUE){ taus <- matrix(taus, 1, n.sim); CIkt <- rowQuantiles(taus, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)}


thetas <- matrix(thetas, 1, n.sim)
CItheta <- rowQuantiles(thetas, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)



} # int


     
} # biv





###########################################
###########################################
###########################################


if(type == "independence"){ # this will not work for ordinal



if(!missing(newdata)){

nu     <- sigma2 <- NA
p1     <- probm(predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
eta2   <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gamlss$coefficients[1:x$X2.d2]



if( !(x$VC$margins[2] %in% cont1par) ){

if( !is.null(x$X3) ){

sigma2 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss$coefficients[(x$X2.d2+1):(x$X2.d2+x$X3.d2)], x$margins[2])$vrb
if(x$margins[2] %in% cont3par) nu <- enu.tr(predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")%*%x$gamlss$coefficients[(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)], x$margins[2])$vrb

                    }

if( is.null(x$X3) ){

if( x$VC$margins[2] %in% c(cont2par,cont3par) ) sigma2 <- x$gamlss$sigma2
if( x$VC$margins[2] %in% c(cont3par) )          nu     <- x$gamlss$nu 

                   }


                                 }


}



if(missing(newdata)){

p1   <- probm(predict.SemiParBIV(x, eq = 1, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
eta2 <- x$gamlss$eta1

sigma2 <- x$gamlss$sigma2
nu     <- x$gamlss$nu 

}



if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*#

if(y1 == 0 || y1 == 1){  

p0  <- 1 - p1                     
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2
p12 <- p0*p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p0

                      }

if(y1 == 1){                                          

p12 <- p1*p2 

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1

           }
           
           
}  #*#         







if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#

p0    <- 1 - p1

#y2rep <- rep(y2, length(eta2))


p2pdf2  <- distrHsATDiscr(y2, eta2, sigma2, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2) 

pdf2  <- p2pdf2$pdf2 
p2    <- p2pdf2$p2  

if(y1 == 0){
p12 <- p0*pdf2
if(cond == 1) p12 <- pdf2
if(cond == 2) p12 <- p0
          }

if(y1 == 1){
p12 <- p1*pdf2
if(cond == 1) p12 <- pdf2
if(cond == 2) p12 <- p1
           }
       
}#*#  







      
if(intervals == TRUE){


bs1 <- rMVN(n.sim, mean = x$gam1$coefficients, sigma=x$gam1$Vp)
bs2 <- rMVN(n.sim, mean = x$gamlss$coefficients, sigma=x$gamlss$Vb)


#############  
# etas
#############  

if(!missing(newdata)){ X1  <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") # this will not work for ordinal
                       X2s <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       } 
  
p1s   <- probm(   X1%*%t(bs1[,1:x$X1.d2]) , x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr   
eta2s <- eta.tr( X2s%*%t(bs2[,1:x$X2.d2]) , x$VC$margins[2])

#############  
# sigmas
#############  

if( x$VC$margins[2] %in% c(cont2par,cont3par) ){ 

      if( is.null(x$X3) )   sigma2.star <- bs2[, x$X2.d2 + 1] 
      
      if( !is.null(x$X3) ){
      
if(!missing(newdata)){ X3s <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")}                     
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
      
                      sigma2.star <- X3s%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X3.d2)]) 

                          }

sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
}    
    
#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs2[,  x$X2.d2 + 2] 
  
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}                      
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}     
  
             nu.st <- X4s%*%t(bs2[,(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
                      }
                      
 nu <- enu.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


#if( is.null(x$X3) ){
#
#if( x$VC$margins[2] %in% c(cont2par,cont3par) ) sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
#if( x$VC$margins[2] %in% cont3par ) nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
#
#                   }
#


if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*#

if(y1 == 0 || y1 == 1){                              
p0s  <- 1 - p1s                    
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2
p12s <- p0s*p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p0s

                      }

if(y1 == 1){                                          

p12s <- p1s*p2s 

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p1s

           }

   
 }#*#   
   
   
   
 

if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#

p0s   <- 1 - p1s


#y2rep <- rep(y2, length(eta2s))


if( x$margins[2] != c("BB") ) y2rep <- rep(y2, length(eta2s))
if( x$margins[2] == c("BB") ) y2rep <- matrix(rep(y2, each = length(eta2s)), ncol = 2, byrow = FALSE)


p2spdf2s <- distrHsATDiscr(y2rep, eta2s, sigma2, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)

pdf2s <- p2spdf2s$pdf2 
p2s   <- p2spdf2s$p2         

if(y1 == 0){
p12s <- p0s*pdf2s
if(cond == 1) p12s <- pdf2s
if(cond == 2) p12s <- p0s
          }

if(y1 == 1){
p12s <- p1s*pdf2s
if(cond == 1) p12s <- pdf2s
if(cond == 2) p12s <- p1s
           }
       
}#*#  
 
 
 
}# int


} # indep


list(p12 = p12, p12s = matrix(p12s, 1, length(p12s)), p1 = p1, p2 = p2, p3 = NULL, CItau = CIkt, tau = tau, CItheta = CItheta, theta = theta)


}

