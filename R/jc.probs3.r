jc.probs3 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link){


#############################################################################################

epsilon <- 0.0000001 
nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- dof <- p12s <- NULL
dof <- x$dof
#############################################################################################
#############################################################################################

if(type == "bivariate"){


if(!missing(newdata)){

p1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "response")
p2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "response")

if(!is.null(x$X3)) theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = 3, newdata = newdata))$teta
if(is.null(x$X3))  theta  <- x$theta 

                     }

if(missing(newdata)){

p1 <- x$p1
p2 <- x$p2
theta <- x$theta 

                    }








if(x$BivD %in% x$BivD2){

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- NA

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1])
if(length(theta) == 1) p12[x$teta.ind1] <- BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta)
                          }
                       
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2])
if(length(theta) == 1) p12[x$teta.ind2] <- BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta)
                          }                       

}


if(!(x$BivD %in% x$BivD2)) p12 <- BiCDF(p1, p2, x$nC, theta, dof) 







if(y1 == 1 && y2 == 1){ 

  p12 <- p12 
  if(cond == 1) p12 <- p12/p1
  if(cond == 2) p12 <- p12/p2

                       }

if(y1 == 1 && y2 == 0){ 

  p12 <- pmax(p1 - p12, epsilon)      
  if(cond == 1) p12 <- p12/p1
  if(cond == 2) p12 <- p12/(1-p2)

                        }

if(y1 == 0 && y2 == 1){ 

  p12 <- pmax(p2 - p12, epsilon)      
  if(cond == 1) p12 <- p12/(1-p1)
  if(cond == 2) p12 <- p12/p2

                      }

if(y1 == 0 && y2 == 0){ 

  p12 <- pmax(1 - p12 - pmax(p1 - p12, epsilon) - pmax(p2 - p12, epsilon) , epsilon)      
  if(cond == 1) p12 <- p12/(1-p1)
  if(cond == 2) p12 <- p12/(1-p2)

                      }



####

if(intervals == TRUE){

bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  

#############  
# etas
############# 


if(!missing(newdata)){ X1  <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2s <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       if(x$Model == "BSS") X2s <- x$X2s else X2s <- x$X2  
                       } 


p1s <- probm( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1])$pr 
p2s <- probm(X2s%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])$pr

#############  
# thetas # this may have to be changed when we have dof
#############  

if(x$Model != "BPO0"){

if(is.null(x$X3))  epds <- bs[,length(x$coefficients)]


if(!is.null(x$X3)){

if(!missing(newdata)){ X3s <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix") }                  
if( missing(newdata)){ if(x$Model == "BSS") X3s <- x$X3s else X3s <- x$X3 } 

epds <- X3s%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])

                   }

            
est.RHOb <- teta.tr(x$VC, epds)$teta

}

if(x$Model == "BPO0") est.RHOb <- rep(0, n.sim )
  
if( is.null(x$X3) ) est.RHOb <- matrix(rep(est.RHOb, each = dim(p1s)[1]), ncol = n.sim, byrow=FALSE)


#############

if(x$VC$BivD %in% c("N","T")) p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, dof, test = FALSE), dim(p1s)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
                      
}


if(!(x$BivD %in% x$BivD2)) p12s <- BiCDF(p1s, p2s, x$nC, est.RHOb, dof, test = FALSE)



}







if(y1 == 1 && y2 == 1){ 

  if(cond == 1) p12s <- p12s/p1s
  if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 1 && y2 == 0){ 

  p12s <- pmax(p1s - p12s, epsilon)     
  if(cond == 1) p12s <- p12s/p1s
  if(cond == 2) p12s <- p12s/(1-p2s)

                      }

if(y1 == 0 && y2 == 1){ 

  p12s <- pmax(p2s - p12s, epsilon)      
  if(cond == 1) p12s <- p12s/(1-p1s)
  if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 0 && y2 == 0){ 

  p12s <- pmax(1 - p12s - pmax(p1s - p12s, epsilon) - pmax(p2s - p12s, epsilon) , epsilon)      
  if(cond == 1) p12s <- p12s/(1-p1s)
  if(cond == 2) p12s <- p12s/(1-p2s)

                       }

          

}




     
}


#############################################################################################
#############################################################################################
#############################################################################################




if(type == "independence"){



if(!missing(newdata)){

p1 <- probm( predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
p2 <- probm( predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gam2$coefficients, x$VC$margins[2])$pr

}



if(missing(newdata)){

p1 <- probm( predict.SemiParBIV(x, eq = 1, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1])$pr
p2 <- probm( predict.SemiParBIV(x, eq = 2, type = "lpmatrix")%*%x$gam2$coefficients, x$VC$margins[2])$pr


}



if(y1 == 1 && y2 == 1){ 

  p12 <- p1*p2     
  if(cond == 1) p12 <- p2
  if(cond == 2) p12 <- p1

                      }

if(y1 == 1 && y2 == 0){ 

  p12 <- p1*(1-p2)     
  if(cond == 1) p12 <- 1-p2
  if(cond == 2) p12 <- p1

                      }

if(y1 == 0 && y2 == 1){ 

  p12 <- (1-p1)*p2     
  if(cond == 1) p12 <- p2
  if(cond == 2) p12 <- 1-p1

                      } 

if(y1 == 0 && y2 == 0){ 

  p12 <- (1-p1)*(1-p2)      
  if(cond == 1) p12 <- (1-p2)
  if(cond == 2) p12 <- (1-p1)

                       }




    
    
if(intervals == TRUE){


bs1 <- rMVN(n.sim, mean = x$gam1$coefficients, sigma=x$gam1$Vp)
bs2 <- rMVN(n.sim, mean = x$gam2$coefficients, sigma=x$gam2$Vp)


if(!missing(newdata)){

p1s <- probm( predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1])$pr
p2s <- probm( predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2])$pr

}



if(missing(newdata)){

p1s <- probm( predict.SemiParBIV(x, eq = 1, type = "lpmatrix")%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1])$pr
p2s <- probm( predict.SemiParBIV(x, eq = 2, type = "lpmatrix")%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2])$pr


}



if(y1 == 1 && y2 == 1){ 

  p12s <- p1s*p2s     
  if(cond == 1) p12s <- p2s
  if(cond == 2) p12s <- p1s

                      }

if(y1 == 1 && y2 == 0){ 

  p12s <- p1s*(1-p2s)     
  if(cond == 1) p12s <- 1-p2s
  if(cond == 2) p12s <- p1s

                      }

if(y1 == 0 && y2 == 1){ 

  p12s <- (1-p1s)*p2s     
  if(cond == 1) p12s <- p2s
  if(cond == 2) p12s <- 1-p1s

                      }

if(y1 == 0 && y2 == 0){ 

  p12s <- (1-p1s)*(1-p2s)      
  if(cond == 1) p12s <- (1-p2s)
  if(cond == 2) p12s <- (1-p1s)

                      }



               } # inter

} # indep
    
    
    
  

list(p12 = p12, p12s = p12s, p1 = p1, p2 = p2)





}



