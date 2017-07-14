jc.probs2 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link){

epsilon <- 0.0000001 
nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- NULL
dof <- x$dof

p12s <- NULL


if(type == "bivariate"){


if(!missing(newdata)){

nu <- sigma2 <- NA

p1   <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "response")
eta2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata)


if( !is.null(x$X3) && x$margins[2] %in% c(cont2par,cont3par) ){

sigma2 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata), x$margins[2])$vrb

if(x$margins[2] %in% cont2par) eq.th <- 4
if(x$margins[2] %in% cont3par){ eq.nu <- 4; eq.th <- 5}

if(x$margins[2] %in% cont3par) nu <- esp.tr(predict.SemiParBIV(x, eq = eq.nu, newdata = newdata), x$margins[2])$vrb

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




if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*#


if(y1 == 0 || y1 == 1){              ### joint prob of y1=0 and y2=x ###

p0  <- 1 - p1         
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2])$p2





if(x$BivD %in% x$BivD2){

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- NA

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof)
if(length(theta) == 1) p12[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof)
                          }
                       
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof)
if(length(theta) == 1) p12[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof)
                          }                       

}



if(!(x$BivD %in% x$BivD2)) p12 <- BiCDF(p0, p2, x$nC, theta, dof)









if(cond == 1) p12 <- p12/p0
if(cond == 2) p12 <- p12/p2

                      }

if(y1 == 1){                         ### joint prob of y1=1 and y2=x ###

p12  <- p2 - p12

if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2

           }
           
           
}#*#





if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#

p0   <- 1 - p1
ppdf <- distrHsATDiscr(y2, eta2, sigma2, nu = 1, x$margins[2], x$VC$y2m)
p2   <- ppdf$p2 
pdf2 <- ppdf$pdf2






if(x$BivD %in% x$BivD2){

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

C1 <- C2 <- NA

if( length(x$teta1) != 0){

if(length(theta) > 1){
C1[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], p2[x$teta.ind1],          nC1, theta[x$teta.ind1], dof)
C2[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], mm(p2[x$teta.ind1]-pdf2[x$teta.ind1]), nC1, theta[x$teta.ind1], dof)
                     }

if(length(theta) == 1){
C1[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], p2[x$teta.ind1],          nC1, theta, dof)
C2[x$teta.ind1] <- BiCDF(p0[x$teta.ind1], mm(p2[x$teta.ind1]-pdf2[x$teta.ind1]), nC1, theta, dof)
                       }

                          }
                       
                       
                       
if( length(x$teta2) != 0){


if(length(theta) > 1){
C1[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], p2[x$teta.ind2],          nC2, theta[x$teta.ind2], dof)
C2[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], mm(p2[x$teta.ind2]-pdf2[x$teta.ind2]), nC2, theta[x$teta.ind2], dof)
                     }

if(length(theta) == 1){
C1[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], p2[x$teta.ind2],          nC2, theta, dof)
C2[x$teta.ind2] <- BiCDF(p0[x$teta.ind2], mm(p2[x$teta.ind2]-pdf2[x$teta.ind2]), nC2, theta, dof)
                       }




                          }                       

}








if(!(x$BivD %in% x$BivD2)){

C1 <- BiCDF(p0, p2,          x$nC, theta, dof)
C2 <- BiCDF(p0, mm(p2-pdf2), x$nC, theta, dof)

}




A <- ifelse(C1 - C2 < epsilon, epsilon, C1 - C2)

if(y1 == 0){
p12 <- A
if(cond == 1) p12 <- p12/p0
if(cond == 2) p12 <- p12/p2
          }

if(y1 == 1){
p12 <- ifelse( pdf2 - A < epsilon, epsilon, pdf2 - A)
if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2
           }
       
}#*#












           
           
if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

#############  
# etas
############# 

if(!missing(newdata)){ X1  <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2s <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       }  


p1s   <- probm( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1])$pr 
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
  
  
if(!missing(newdata)){ X4s <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
                epds <- X4s%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)])
  
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

if( x$VC$margins[2] %in% cont2par ){  

      if( is.null(x$X3) )   sigma2.star <- bs[, x$X1.d2 + x$X2.d2 + 1] 
      if( !is.null(x$X3) ) {
      
if(!missing(newdata)){ X3s <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")}
                       
if( missing(newdata)){ if(x$VC$ccss == "yes") X3s <- x$X3s else X3s <- x$X3}        
                
            sigma2.star <- X3s%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 
  
                           }
  
sigma2 <- esp.tr(sigma2.star, x$VC$margins[2])$vrb   
    
}    
    
#############  
# NUs
#############    
  
if( x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3)  ) nu.st <- bs[, x$X1.d2 + x$X2.d2 + 2] # t(as.matrix(bs[,  x$X1.d2 + x$X2.d2 + 2]))
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}                    
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}    
  
              nu.st <- X4s%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
                      }
  
 nu <- esp.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


if( is.null(x$X3) ){

est.RHOb <- matrix(rep(est.RHOb, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }




if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*#

if(y1 == 0 || y1 == 1){             
p0s  <- 1 - p1s         
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2])$p2

if(x$VC$BivD %in% c("N","T")) p12s <- matrix(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), dim(p1s)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- BiCDF(p0s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- BiCDF(p0s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
                      
                        }

if(!(x$BivD %in% x$BivD2)) p12s <- BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE)




}






if(cond == 1) p12s <- p12s/p0s
if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 1){                         
  
p12s  <- p2s - p12s

if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s

           }
 
}#*# 
 
 
 
 
if(x$margins[2] %in% c(x$VC$m2d, x$VC$m1d)){#*#

p0s   <- 1 - p1s
ppdf  <- distrHsATDiscr(y2, eta2s, sigma2, nu = 1, x$margins[2], x$VC$y2m) # not sure about y2m
p2s   <- ppdf$p2 
pdf2s <- ppdf$pdf2

if(x$VC$BivD %in% c("N","T")) C1s <- matrix(BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE), dim(p1s)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

C1s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])

if( length(x$teta1) != 0) C1s[x$teta.ind1,] <- BiCDF(p0s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) C1s[x$teta.ind2,] <- BiCDF(p0s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
                      
                        }

if(!(x$BivD %in% x$BivD2)) C1s <- BiCDF(p0s, p2s, x$nC, est.RHOb, dof, test = FALSE)



}













if(x$VC$BivD %in% c("N","T")) C2s <- matrix(BiCDF(p0s, mm(p2s-pdf2s), x$nC, est.RHOb, dof, test = FALSE), dim(p1s)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

C2s <- matrix(NA, ncol = n.sim, nrow = dim(p0s)[1])

if( length(x$teta1) != 0) C2s[x$teta.ind1,] <- BiCDF(p0s[x$teta.ind1,], mm(p2s[x$teta.ind1,]-pdf2s[x$teta.ind1,]), nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) C2s[x$teta.ind2,] <- BiCDF(p0s[x$teta.ind2,], mm(p2s[x$teta.ind2,]-pdf2s[x$teta.ind2,]), nC2, -est.RHOb[x$teta.ind2,])
                      
                        }

if(!(x$BivD %in% x$BivD2)) C2s <- BiCDF(p0s, mm(p2s-pdf2s), x$nC, est.RHOb, dof, test = FALSE)


}




As <- ifelse(C1s - C2s < epsilon, epsilon, C1s - C2s)




if(y1 == 0){
p12s <- As
if(cond == 1) p12s <- p12s/p0s
if(cond == 2) p12s <- p12s/p2s
          }

if(y1 == 1){
p12s <- ifelse( pdf2s - As < epsilon, epsilon, pdf2s - As)
if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s
           }
       
}#*# 
 
           

} # int




     
} # biv





###########################################
###########################################
###########################################


if(type == "independence"){



if(!missing(newdata)){

nu     <- sigma2 <- NA
p1     <- probm(predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
eta2   <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gamlss$coefficients[1:x$X2.d2]



if( !(x$VC$margins[2] %in% cont1par) ){

if( !is.null(x$X3) ){

sigma2 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss$coefficients[(x$X2.d2+1):(x$X2.d2+x$X3.d2)], x$margins[2])$vrb
if(x$margins[2] %in% cont3par) nu <- esp.tr(predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")%*%x$gamlss$coefficients[(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)], x$margins[2])$vrb

                    }

if( is.null(x$X3) ){

sigma2 <- x$gamlss$sigma2
nu     <- x$gamlss$nu 

                   }


                                 }


}



if(missing(newdata)){

p1   <- probm(predict.SemiParBIV(x, eq = 1, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
eta2 <- x$gamlss$eta1

sigma2 <- x$gamlss$sigma2
nu     <- x$gamlss$nu 

}



if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*#

if(y1 == 0 || y1 == 1){  

p0  <- 1 - p1                     
p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2])$p2
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

p0   <- 1 - p1
p2 <- distrHsATDiscr(y2, eta2, sigma2, nu = 1, x$margins[2], x$VC$y2m)$pdf2 # note pdf2 here! 

if(y1 == 0){
p12 <- p0*p2
if(cond == 1) p12 <- p12/p0
if(cond == 2) p12 <- p12/p2
          }

if(y1 == 1){
p12 <- p1*p2
if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2
           }
       
}#*#  







      
if(intervals == TRUE){


bs1 <- rMVN(n.sim, mean = x$gam1$coefficients, sigma=x$gam1$Vp)
bs2 <- rMVN(n.sim, mean = x$gamlss$coefficients, sigma=x$gamlss$Vb)


#############  
# etas
#############  

if(!missing(newdata)){ X1  <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2s <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$VC$ccss == "yes") X2s <- x$X2s else X2s <- x$X2  
                       } 
  
p1s   <- probm(   X1%*%t(bs1[,1:x$X1.d2]) , x$VC$margins[1])$pr   
eta2s <- eta.tr( X2s%*%t(bs2[,1:x$X2.d2]) , x$VC$margins[2])

#############  
# sigmas
#############  

if( x$VC$margins[2] %in% cont2par ){ 

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
    
  if( is.null(x$X3)  ) nu.st <- bs2[,  x$X2.d2 + 2] # t(as.matrix(bs2[,  x$X2.d2 + 2]))
  
  if( !is.null(x$X3) ){
  
if(!missing(newdata)){ X4s <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}                      
if( missing(newdata)){ if(x$VC$ccss == "yes") X4s <- x$X4s else X4s <- x$X4}     
  
             nu.st <- X4s%*%t(bs2[,(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
                      }
                      
 nu <- esp.tr(nu.st, x$VC$margins[2])$vrb   
  
} 


#################


if( is.null(x$X3) ){

sigma2   <- matrix(rep(sigma2, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)
nu       <- matrix(rep(nu, each = dim(eta2s)[1]), ncol = n.sim, byrow=FALSE)

                   }



if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)){#*#

if(y1 == 0 || y1 == 1){                              
p0s  <- 1 - p1s                     
p2s  <- distrHsAT(y2, eta2s, sigma2, nu, x$margins[2])$p2
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
pdf2s <- distrHsATDiscr(y2, eta2s, sigma2, nu = 1, x$margins[2], x$VC$y2m)$pdf2 # not sure about y2m

if(y1 == 0){
p12s <- p0s*pdf2s
if(cond == 1) p12s <- p12s/p0s
if(cond == 2) p12s <- p12s/p2s
          }

if(y1 == 1){
p12s <- p1s*pdf2s
if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s
           }
       
}#*#  
 
 
 
}# int


} # indep


list(p12 = p12, p12s = p12s, p1 = p1, p2 = p2)


}



