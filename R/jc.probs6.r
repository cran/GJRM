jc.probs6 <- function(x, y1, y2, y3, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link){

#############################################################################################

epsilon <- 0.0000001 
nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- dof <- p12s <- p111 <- NULL
dof <- x$dof
lf <- length(x$coefficients)
theta12 <- theta13 <- theta23 <- 0
epd12s <- epd13s <- epd23s <- p123s <- 0
CIkt <- tau <- NULL
#############################################################################################
#############################################################################################


X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") 
X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")
X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")


if(!is.null(x$X4)){

X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix") 
X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")
X6 <- predict.SemiParBIV(x, eq = 6, newdata = newdata, type = "lpmatrix")

}





if(type == "joint"){

if(!missing(newdata)){ # !missing(newdata)

b1 <- x$coefficients[1:x$X1.d2]
b2 <- x$coefficients[(x$X1.d2+1):(x$X1.d2+x$X2.d2)]
b3 <- x$coefficients[(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]

p1 <- probm(X1%*%b1, x$margins[1])$pr
p2 <- probm(X2%*%b2, x$margins[2])$pr
p3 <- probm(X3%*%b3, x$margins[3])$pr


if(!is.null(x$X4)){

b4 <- x$coefficients[(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)]
b5 <- x$coefficients[(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)]
b6 <- x$coefficients[(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)]

theta12s <- X4%*%b4
theta13s <- X5%*%b5
theta23s <- X6%*%b6


  Sigma <- list()    
      for(i in 1:dim(X1)[1]){
        Sigma[[i]] <- PosDefCor(Sigma = NULL, Chol = TRUE, theta12.st = theta12s[i], theta13.st = theta13s[i], theta23.st = theta23s[i]) 
        theta12[i] <- Sigma[[i]][1,2]
        theta13[i] <- Sigma[[i]][1,3]
        theta23[i] <- Sigma[[i]][2,3]   
                      }
}



if(is.null(x$X4)){ 
  theta12 <- x$theta12
  theta13 <- x$theta13
  theta23 <- x$theta23
                 }

                     }# !missing(newdata)
                     
                     

if(missing(newdata)){

p1 <- x$fit$p1
p2 <- x$fit$p2
p3 <- x$fit$p3
theta12 <- x$theta12
theta13 <- x$theta13
theta23 <- x$theta23

                    }



###############################

p1 <- mm(as.numeric(p1))
p2 <- mm(as.numeric(p2))
p3 <- mm(as.numeric(p3))

mar1 <- qnorm(p1)
mar2 <- qnorm(p2)
mar3 <- qnorm(p3)

theta12 <- mmf(as.numeric(theta12))
theta13 <- mmf(as.numeric(theta13))
theta23 <- mmf(as.numeric(theta23))

############################### all triv probabilities
  
  p11 <- mm( pbinorm( mar1, mar2, cov12 = theta12) )
  p13 <- mm( pbinorm( mar1, mar3, cov12 = theta13) )
  p23 <- mm( pbinorm( mar2, mar3, cov12 = theta23) )
  
  
  if(is.null(x$X4)){
  
     Sigma <-  matrix( c( 1,       theta12, theta13,
                         theta12,        1, theta23,
                         theta13,  theta23,        1), 3 , 3)
  
    for(i in 1:dim(X1)[1]) p111[i] <- mm( pmnorm(x = c(mar1[i], mar2[i], mar3[i]), varcov = Sigma)[1] ) 
    
                    }
  
  
  if(!is.null(x$X4)){
  
    for(i in 1:dim(X1)[1]){
  
           Sigma <-  matrix( c( 1,       theta12[i], theta13[i],
                               theta12[i],        1, theta23[i],
                               theta13[i],  theta23[i],        1), 3 , 3)
  
             p111[i] <- mm( pmnorm(x = c(mar1[i], mar2[i], mar3[i]), varcov = Sigma)[1] ) 
  
                   }
  
  
  
                    }  
  
  
  p011 <- mm(p23 - p111)
  p101 <- mm(p13 - p111)
  p110 <- mm(p11 - p111)
  p100 <- mm(p1 - p11 - p101)
  p010 <- mm(p2 - p11 - p011)
  p001 <- mm(p3 - p23 - p101)
  p000 <- mm(1 - p111 - p011 - p101 - p110 - p001 - p010 - p100)  
  
###############################
  
if(y1 == 1 && y2 == 1 && y3 == 1){ 

  p123 <- p111 
  if(cond == 1) p123 <- p123/p1
  if(cond == 2) p123 <- p123/p2
  if(cond == 3) p123 <- p123/p3

                                  }

if(y1 == 0 && y2 == 1 && y3 == 1){ 

  p123 <- p011 
  if(cond == 1) p123 <- p123/(1-p1)
  if(cond == 2) p123 <- p123/p2
  if(cond == 3) p123 <- p123/p3

                                  }
                                  
if(y1 == 1 && y2 == 0 && y3 == 1){ 

  p123 <- p101 
  if(cond == 1) p123 <- p123/p1
  if(cond == 2) p123 <- p123/(1-p2)
  if(cond == 3) p123 <- p123/p3

                                  }

if(y1 == 1 && y2 == 1 && y3 == 0){ 

  p123 <- p110 
  if(cond == 1) p123 <- p123/p1
  if(cond == 2) p123 <- p123/p2
  if(cond == 3) p123 <- p123/(1-p3)

                                  }
                                  
if(y1 == 1 && y2 == 0 && y3 == 0){ 

  p123 <- p100 
  if(cond == 1) p123 <- p123/p1
  if(cond == 2) p123 <- p123/(1-p2)
  if(cond == 3) p123 <- p123/(1-p3)

                                  }
                                  
if(y1 == 0 && y2 == 1 && y3 == 0){ 

  p123 <- p010 
  if(cond == 1) p123 <- p123/(1-p1)
  if(cond == 2) p123 <- p123/p2
  if(cond == 3) p123 <- p123/(1-p3)

                                  }

if(y1 == 0 && y2 == 0 && y3 == 1){ 

  p123 <- p001 
  if(cond == 1) p123 <- p123/(1-p1)
  if(cond == 2) p123 <- p123/(1-p2)
  if(cond == 3) p123 <- p123/p3

                                  }                                  

if(y1 == 0 && y2 == 0 && y3 == 0){ 

  p123 <- p000 
  if(cond == 1) p123 <- p123/(1-p1)
  if(cond == 2) p123 <- p123/(1-p2)
  if(cond == 3) p123 <- p123/(1-p3)

                                  }



####
####
####


if(intervals == TRUE){


bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

p1s <- probm(X1%*%t(bs[,1:x$X1.d2]),                                     x$margins[1])$pr
p2s <- probm(X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]),                 x$margins[2])$pr
p3s <- probm(X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]), x$margins[3])$pr

mar1s <- qnorm(p1s)
mar2s <- qnorm(p2s)
mar3s <- qnorm(p3s)




if(is.null(x$X4))  epd12s <- epd13s <- epd23s <- NA  
if(!is.null(x$X4)) epd12s <- epd13s <- epd23s <- p111s <- matrix(NA, dim(X1)[1], n.sim)  


if(x$VC$Chol == FALSE){
  
 epd12s <- teta.tr(x$VC, bs[, lf-2])$teta  
 epd13s <- teta.tr(x$VC, bs[, lf-1])$teta  
 epd23s <- teta.tr(x$VC, bs[, lf]  )$teta  
  
}


  
  
if(x$VC$Chol == TRUE){ 
  
  if(is.null(x$X4)){  
 
     for(i in 1:n.sim){ 

      Sigma <- PosDefCor(Sigma = 1, Chol = TRUE, theta12.st = bs[i, lf-2], theta13.st = bs[i, lf-1], theta23.st = bs[i, lf])
      
      epd12s[i] <- Sigma[1,2]
      epd13s[i] <- Sigma[1,3]
      epd23s[i] <- Sigma[2,3]  
      
                      }
                    }    
    
    
   if(!is.null(x$X4)){  
 
eta4s <- X4%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)])
eta5s <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
eta6s <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)]) 
 
      
    for(i in 1:dim(X1)[1]){    

       for(j in 1:n.sim){
       
         Sigma <- PosDefCor(Sigma = 1, Chol = TRUE, theta12.st = eta4s[i, j], theta13.st = eta5s[i, j], theta23.st = eta6s[i, j])
     
         epd12s[i,j] <- Sigma[1,2]
         epd13s[i,j] <- Sigma[1,3]
         epd23s[i,j] <- Sigma[2,3] 
      
                        }
                    }     
                      }       
}
        



if( is.null(x$X4) ){

           epd12s <- matrix(rep(epd12s, each = dim(p1s)[1]), ncol = n.sim, byrow = FALSE)
           epd13s <- matrix(rep(epd13s, each = dim(p1s)[1]), ncol = n.sim, byrow = FALSE)
           epd23s <- matrix(rep(epd23s, each = dim(p1s)[1]), ncol = n.sim, byrow = FALSE)

                    }





###### p111s

    for(i in 1:dim(X1)[1]){    

       for(j in 1:n.sim){
       
           Sigma <-  matrix( c( 1,           epd12s[i,j], epd13s[i,j],
                                epd12s[i,j],           1, epd23s[i,j],
                                epd13s[i,j], epd23s[i,j],            1), 3 , 3)

           p111s[i, j] <- mm( pmnorm(x = c(mar1s[i,j], mar2s[i,j], mar3s[i,j]), varcov = Sigma)[1] ) 
                      
                         }
                    }                  
######                 

  p11s <- matrix(mm( pbinorm( mar1s, mar2s, cov12 = epd12s) ), dim(p1s)[1], n.sim) 
  p13s <- matrix(mm( pbinorm( mar1s, mar3s, cov12 = epd13s) ), dim(p1s)[1], n.sim) 
  p23s <- matrix(mm( pbinorm( mar2s, mar3s, cov12 = epd23s) ), dim(p1s)[1], n.sim) 


  p011s <- mm(p23s - p111s)
  p101s <- mm(p13s - p111s)
  p110s <- mm(p11s - p111s)
  p100s <- mm(p1s - p11s - p101s)
  p010s <- mm(p2s - p11s - p011s)
  p001s <- mm(p3s - p23s - p101s)
  p000s <- mm(1 - p111s - p011s - p101s - p110s - p001s - p010s - p100s)  
  
###############################
  
if(y1 == 1 && y2 == 1 && y3 == 1){ 

  p123s <- p111s 
  if(cond == 1) p123s <- p123s/p1s
  if(cond == 2) p123s <- p123s/p2s
  if(cond == 3) p123s <- p123s/p3s

                                  }

if(y1 == 0 && y2 == 1 && y3 == 1){ 

  p123s <- p011s 
  if(cond == 1) p123s <- p123s/(1-p1s)
  if(cond == 2) p123s <- p123s/p2s
  if(cond == 3) p123s <- p123s/p3s

                                  }
                                  
if(y1 == 1 && y2 == 0 && y3 == 1){ 

  p123s <- p101s 
  if(cond == 1) p123s <- p123s/p1s
  if(cond == 2) p123s <- p123s/(1-p2s)
  if(cond == 3) p123s <- p123s/p3s

                                  }

if(y1 == 1 && y2 == 1 && y3 == 0){ 

  p123s <- p110s 
  if(cond == 1) p123s <- p123s/p1s
  if(cond == 2) p123s <- p123s/p2s
  if(cond == 3) p123s <- p123s/(1-p3s)

                                  }
                                  
if(y1 == 1 && y2 == 0 && y3 == 0){ 

  p123s <- p100s 
  if(cond == 1) p123s <- p123s/p1s
  if(cond == 2) p123s <- p123s/(1-p2s)
  if(cond == 3) p123s <- p123s/(1-p3s)

                                  }
                                  
if(y1 == 0 && y2 == 1 && y3 == 0){ 

  p123s <- p010s 
  if(cond == 1) p123s <- p123s/(1-p1s)
  if(cond == 2) p123s <- p123s/p2s
  if(cond == 3) p123s <- p123s/(1-p3s)

                                  }

if(y1 == 0 && y2 == 0 && y3 == 1){ 

  p123s <- p001s 
  if(cond == 1) p123s <- p123s/(1-p1s)
  if(cond == 2) p123s <- p123s/(1-p2s)
  if(cond == 3) p123s <- p123s/p3s

                                  }                                  

if(y1 == 0 && y2 == 0 && y3 == 0){ 

  p123s <- p000s 
  if(cond == 1) p123s <- p123s/(1-p1s)
  if(cond == 2) p123s <- p123s/(1-p2s)
  if(cond == 3) p123s <- p123s/(1-p3s)

                                  }
   
}

}


#############################################################################################
#############################################################################################
#############################################################################################




if(type == "independence"){

b1 <- x$gam1$coefficients
b2 <- x$gam2$coefficients
b3 <- x$gam3$coefficients

p1 <- probm(X1%*%b1, x$margins[1])$pr
p2 <- probm(X2%*%b2, x$margins[2])$pr
p3 <- probm(X3%*%b3, x$margins[3])$pr

###############################

p1 <- mm(as.numeric(p1))
p2 <- mm(as.numeric(p2))
p3 <- mm(as.numeric(p3))

###############################
  
if(y1 == 1 && y2 == 1 && y3 == 1){ 

  p123 <- p1*p2*p3 
  if(cond == 1) p123 <- p2*p3
  if(cond == 2) p123 <- p1*p3
  if(cond == 3) p123 <- p1*p2

                                  }

if(y1 == 0 && y2 == 1 && y3 == 1){ 

  p123 <- (1-p1)*p2*p3  
  if(cond == 1) p123 <- p2*p3 
  if(cond == 2) p123 <- (1-p1)*p3
  if(cond == 3) p123 <- (1-p1)*p2

                                  }
                                  
if(y1 == 1 && y2 == 0 && y3 == 1){ 

  p123 <- p1*(1-p2)*p3 
  if(cond == 1) p123 <- (1-p2)*p3 
  if(cond == 2) p123 <- p1*p3 
  if(cond == 3) p123 <- p1*(1-p2) 

                                  }

if(y1 == 1 && y2 == 1 && y3 == 0){ 

  p123 <- p1*p2*(1-p3) 
  if(cond == 1) p123 <- p2*(1-p3)  
  if(cond == 2) p123 <- p1*(1-p3)  
  if(cond == 3) p123 <- p1*p2  

                                  }
                                  
if(y1 == 1 && y2 == 0 && y3 == 0){ 

  p123 <- p1*(1-p2)*(1-p3) 
  if(cond == 1) p123 <- (1-p2)*(1-p3) 
  if(cond == 2) p123 <- p1*(1-p3) 
  if(cond == 3) p123 <- p1*(1-p2) 

                                  }
                                  
if(y1 == 0 && y2 == 1 && y3 == 0){ 

  p123 <- (1-p1)*p2*(1-p3)  
  if(cond == 1) p123 <- p2*(1-p3) 
  if(cond == 2) p123 <- (1-p1)*(1-p3) 
  if(cond == 3) p123 <- (1-p1)*p2 

                                  }

if(y1 == 0 && y2 == 0 && y3 == 1){ 

  p123 <- (1-p1)*(1-p2)*p3 
  if(cond == 1) p123 <- (1-p2)*p3 
  if(cond == 2) p123 <- (1-p1)*p3 
  if(cond == 3) p123 <- (1-p1)*(1-p2) 

                                  }                                  

if(y1 == 0 && y2 == 0 && y3 == 0){ 

  p123 <- (1-p1)*(1-p2)*(1-p3)  
  if(cond == 1) p123 <- (1-p2)*(1-p3) 
  if(cond == 2) p123 <- (1-p1)*(1-p3) 
  if(cond == 3) p123 <- (1-p1)*(1-p2) 

                                  }



####
####
####


if(intervals == TRUE){


bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

p1s <- probm(X1%*%t(bs[,1:x$X1.d2]),                                     x$margins[1])$pr
p2s <- probm(X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]),                 x$margins[2])$pr
p3s <- probm(X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]), x$margins[3])$pr


  
###############################
  
  if(y1 == 1 && y2 == 1 && y3 == 1){ 
  
    p123s <- p1s*p2s*p3s 
    if(cond == 1) p123s <- p2s*p3s
    if(cond == 2) p123s <- p1s*p3s
    if(cond == 3) p123s <- p1s*p2s
  
                                    }
  
  if(y1 == 0 && y2 == 1 && y3 == 1){ 
  
    p123s <- (1-p1s)*p2s*p3s  
    if(cond == 1) p123s <- p2s*p3s 
    if(cond == 2) p123s <- (1-p1s)*p3s
    if(cond == 3) p123s <- (1-p1s)*p2s
  
                                    }
                                    
  if(y1 == 1 && y2 == 0 && y3 == 1){ 
  
    p123s <- p1s*(1-p2s)*p3s 
    if(cond == 1) p123s <- (1-p2s)*p3s 
    if(cond == 2) p123s <- p1s*p3s 
    if(cond == 3) p123s <- p1s*(1-p2s) 
  
                                    }
  
  if(y1 == 1 && y2 == 1 && y3 == 0){ 
  
    p123s <- p1s*p2s*(1-p3s) 
    if(cond == 1) p123s <- p2s*(1-p3s)  
    if(cond == 2) p123s <- p1s*(1-p3s)  
    if(cond == 3) p123s <- p1s*p2s  
  
                                    }
                                    
  if(y1 == 1 && y2 == 0 && y3 == 0){ 
  
    p123s <- p1s*(1-p2s)*(1-p3s) 
    if(cond == 1) p123s <- (1-p2s)*(1-p3s) 
    if(cond == 2) p123s <- p1s*(1-p3s) 
    if(cond == 3) p123s <- p1s*(1-p2s) 
  
                                    }
                                    
  if(y1 == 0 && y2 == 1 && y3 == 0){ 
  
    p123s <- (1-p1s)*p2s*(1-p3s)  
    if(cond == 1) p123s <- p2s*(1-p3s) 
    if(cond == 2) p123s <- (1-p1s)*(1-p3s) 
    if(cond == 3) p123s <- (1-p1s)*p2s 
  
                                    }
  
  if(y1 == 0 && y2 == 0 && y3 == 1){ 
  
    p123s <- (1-p1s)*(1-p2s)*p3s 
    if(cond == 1) p123s <- (1-p2s)*p3s 
    if(cond == 2) p123s <- (1-p1s)*p3s 
    if(cond == 3) p123s <- (1-p1s)*(1-p2s) 
  
                                    }                                  
  
  if(y1 == 0 && y2 == 0 && y3 == 0){ 
  
    p123s <- (1-p1s)*(1-p2s)*(1-p3s)  
    if(cond == 1) p123s <- (1-p2s)*(1-p3s) 
    if(cond == 2) p123s <- (1-p1s)*(1-p3s) 
    if(cond == 3) p123s <- (1-p1s)*(1-p2s) 
  
                                    }




} # indep


}
    
    
    
  

list(p12 = p123, p12s = p123s, p1 = p1, p2 = p2, p3 = p3, 
     theta12 = theta12, theta13 = theta13, theta23 = theta23, 
     theta12s = epd12s, theta13s = epd13s, theta23s = epd23s)





}



