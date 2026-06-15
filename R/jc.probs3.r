jc.probs3 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = FALSE){


#############################################################################################

nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- dof <- p12s <- CIkt <- tau <- theta <- CItheta <- NULL
dof <- x$dof

p12s  <- matrix(0, 1, 2) 

#############################################################################################
#############################################################################################



nC   <- x$nC
BivD <- x$BivD




if(type == "joint"){


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




p1 <- as.numeric(p1)
p2 <- as.numeric(p2)
theta <- as.numeric(theta) 





if(x$BivD %in% x$BivD2){ 

if(theta >= 0){ nC <- x$fit$nC1; BivD <- x$fit$Cop1 } else{ nC <- x$fit$nC2; BivD <- x$fit$Cop2} 

}




p12 <- mm(BiCDF(p1, p2, nC, theta, dof), min.pr = min.pr, max.pr = max.pr  ) 







if(y1 == 1 && y2 == 1){ 

  p12 <- p12 
  if(cond == 1) p12 <- p12/p1
  if(cond == 2) p12 <- p12/p2

                       }

if(y1 == 1 && y2 == 0){ 

  p12 <- pmax(p1 - p12, min.pr)      
  if(cond == 1) p12 <- p12/p1
  if(cond == 2) p12 <- p12/(1-p2)

                        }

if(y1 == 0 && y2 == 1){ 

  p12 <- pmax(p2 - p12, min.pr)      
  if(cond == 1) p12 <- p12/(1-p1)
  if(cond == 2) p12 <- p12/p2

                      }

if(y1 == 0 && y2 == 0){ 

  p12 <- pmax(1 - p12 - pmax(p1 - p12, min.pr) - pmax(p2 - p12, min.pr) , min.pr)      
  if(cond == 1) p12 <- p12/(1-p1)
  if(cond == 2) p12 <- p12/(1-p2)

                      }


####

# kendalls' tau

tau   <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = tau.res, t.theta = theta)$tau
theta <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = FALSE,   t.theta = theta)$theta
     

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


p1s <- probm( X1%*%t(bs[,1:x$X1.d2]), x$VC$margins[1],                      min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr 
p2s <- probm(X2s%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

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
  


#############


p12s <- mm(BiCDF(p1s, p2s, nC, est.RHOb, dof, test = FALSE), min.pr = min.pr, max.pr = max.pr  )




if(y1 == 1 && y2 == 1){ 

  if(cond == 1) p12s <- p12s/p1s
  if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 1 && y2 == 0){ 

  p12s <- pmax(p1s - p12s, min.pr)     
  if(cond == 1) p12s <- p12s/p1s
  if(cond == 2) p12s <- p12s/(1-p2s)

                      }

if(y1 == 0 && y2 == 1){ 

  p12s <- pmax(p2s - p12s, min.pr)      
  if(cond == 1) p12s <- p12s/(1-p1s)
  if(cond == 2) p12s <- p12s/p2s

                      }

if(y1 == 0 && y2 == 0){ 

  p12s <- pmax(1 - p12s - pmax(p1s - p12s, min.pr) - pmax(p2s - p12s, min.pr) , min.pr)      
  if(cond == 1) p12s <- p12s/(1-p1s)
  if(cond == 2) p12s <- p12s/(1-p2s)

                       }







  
ass.msR <- theta2tau(x$VC$BivD, x$VC$nCa, est.RHOb, tau.res = tau.res, t.theta = theta)
if(tau.res == TRUE) taus <- ass.msR$tau
thetas <- ass.msR$theta


                                                  
if(tau.res == TRUE){taus <- matrix(taus, 1, n.sim); CIkt <- rowQuantiles(taus, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)}


thetas <- matrix(thetas, 1, n.sim)
CItheta <- rowQuantiles(thetas, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

          

}

    
}


#############################################################################################
#############################################################################################
#############################################################################################




if(type == "independence"){



if(!missing(newdata)){

p1 <- probm( predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
p2 <- probm( predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gam2$coefficients, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

}



if(missing(newdata)){

p1 <- probm( predict.SemiParBIV(x, eq = 1, type = "lpmatrix")%*%x$gam1$coefficients, x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
p2 <- probm( predict.SemiParBIV(x, eq = 2, type = "lpmatrix")%*%x$gam2$coefficients, x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr


}



p1 <- as.numeric(p1)
p2 <- as.numeric(p2)


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

p1s <- probm( predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
p2s <- probm( predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr

}



if(missing(newdata)){

p1s <- probm( predict.SemiParBIV(x, eq = 1, type = "lpmatrix")%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr
p2s <- probm( predict.SemiParBIV(x, eq = 2, type = "lpmatrix")%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr)$pr


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
    
    
    
  

list(p12 = p12, p12s = matrix(p12s, 1, length(p12s)), p1 = p1, p2 = p2, p3 = NULL, CItau = CIkt, tau = tau, CItheta = CItheta, theta = theta)





}



