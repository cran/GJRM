pred.mvt <- function(x, eq, fun = "mean", n.sim = 100, prob.lev = 0.05, smooth.no = NULL, ...){


 if( !(fun %in% c("mean", "variance", "tau")) ) stop("The value for fun can either be mean, variance or tau.")

 if( fun == "tau" && x$VC$univ.gamls == TRUE) stop("The value for fun can either be mean or variance.")
 
 
 
 
 if(fun != "tau"){
 
 
 
 
 if(x$VC$univ.gamls == FALSE){############### 
 
 if(missing(eq) && x$VC$univ.gamls == FALSE) stop("You must provide the equation number.")
 if(!(eq %in% c(1, 2)) && x$VC$univ.gamls == FALSE) stop("The value for eq can be either 1 or 2.")

 if(!(x$margins[1] %in% c(x$VC$m2, x$VC$m2d, x$VC$bl, x$VC$m1d)) && !(x$margins[2] %in% c(x$VC$m2, x$VC$m3, x$VC$m2d)) ) stop("This is currently implemented for margins with two parameters.\nGet in touch to check progress on the other cases.") 

 if(eq == 1) mar <- x$margins[1]
 if(eq == 2) mar <- x$margins[2] 
 
  if( mar %in% c("DGP") ) stop("Function not ready yet for this case.")

 
 
 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)
 
 
 
if( x$margins[1] %in% c(x$VC$m2, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m2d) ){ ###############
 
  
 if(eq == 1){ind1 <- (1:x$X1.d2);                       ind2 <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) } 
 if(eq == 2){ind1 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2); ind2 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2) }  

if(eq == 1){

 Xfit1 <- predict(x, eq = 1, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 1, type = "link", ...)
 
 Xfit2 <- predict(x, eq = 3, type = "lpmatrix", ...) 
 fit2  <- predict(x, eq = 3, type = "link", ...) 
 
           }

if(eq == 2){

 Xfit1 <- predict(x, eq = 2, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 2, type = "link", ...)

 Xfit2 <- predict(x, eq = 4, type = "lpmatrix", ...)  
 fit2  <- predict(x, eq = 4, type = "link", ...)  
 
           }


if(!is.null(smooth.no)){


if(eq == 1){

Xfit1[, -c(x$gam1$smooth[[smooth.no]]$first.para:x$gam1$smooth[[smooth.no]]$last.para)] <- 0 
Xfit2[, -c(x$gam3$smooth[[smooth.no]]$first.para:x$gam3$smooth[[smooth.no]]$last.para)] <- 0 

           }

if(eq == 2){

Xfit1[, -c(x$gam2$smooth[[smooth.no]]$first.para:x$gam2$smooth[[smooth.no]]$last.para)] <- 0 
Xfit2[, -c(x$gam4$smooth[[smooth.no]]$first.para:x$gam4$smooth[[smooth.no]]$last.para)] <- 0 

           }

fit1 <- Xfit1%*%x$coefficients[ind1]
fit2 <- Xfit2%*%x$coefficients[ind2]

}




fit1Sim <- Xfit1%*%t(bs[, ind1])
fit2Sim <- Xfit2%*%t(bs[, ind2])





}############




if( x$margins[1] %in% c(x$VC$bl, x$VC$m1d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3, x$VC$m2d) ){ #############
 


  
if (!is.null(x$VC$K1)) {
    
        K1  <- x$VC$K1
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- CLM.shift2 <- 0
}  

 
 
 
  
 if(eq == 1){ind1 <- (1:x$X1.d2) } # for ordinal this is not correct
 
 if(eq == 2){ind1 <- (x$X1.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + CLM.shift2); 
 
     if(is.null(x$X3)) ind2 <- x$X1.d2 + x$X2.d2 + 1 + CLM.shift2 else ind2 <- (x$X1.d2 + x$X2.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + CLM.shift2) 
     if(is.null(x$X4)) ind3 <- x$X1.d2 + x$X2.d2 + x$X3.d2 + 1 + CLM.shift2 else ind3 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + CLM.shift2)
     
             
             }  

if(eq == 1){
 Xfit1 <- predict(x, eq = 1, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 1, type = "link", ...)
           }

if(eq == 2){


 Xfit1 <- predict(x, eq = 2, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 2, type = "link", ...)
 
 if(!is.null(x$X3)){
 
   Xfit2 <- predict(x, eq = 3, type = "lpmatrix", ...) 
   fit2  <- predict(x, eq = 3, type = "link", ...)  
   
                   }
                   
 if(!is.null(x$X4)){
 
   Xfit3 <- predict(x, eq = 4, type = "lpmatrix", ...) 
   fit3  <- predict(x, eq = 4, type = "link", ...)  
   
                   }                   
 
 if(is.null(x$X3)){
 
   Xfit2 <- matrix(1, nrow = length(fit1), ncol = 1) 
   fit2  <- x$coefficients["sigma2.star"]
   
                   } 
                   
 if(is.null(x$X4)){
 
   Xfit3 <- matrix(1, nrow = length(fit1), ncol = 1) 
   fit3  <- x$coefficients["nu.star"]
   
                   }                    
 
           }
           
           
fit1Sim <- Xfit1%*%t(bs[, ind1])

if(eq == 2){ 

       if(!is.null(x$X3)) fit2Sim <- Xfit2%*%t(bs[, ind2]) else fit2Sim <- bs[, ind2]
       if(!is.null(x$X4)) fit3Sim <- Xfit3%*%t(bs[, ind3]) else fit3Sim <- bs[, ind3]
       
           }


}##################







}##############################################






if(x$VC$univ.gamls == TRUE){############### 



mar <- x$margins[1]

if( mar %in% c("DGP") ) stop("Function not ready yet for this case.")

 
bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)

 
 

if(mar %in% x$VC$m2){

ind1 <- (1:x$X1.d2)
if( !is.null(x$X2) ) ind2 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2) else ind2 <- x$X1.d2 + 1    

 Xfit1 <- predict(x, eq = 1, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 1, type = "link", ...)
 fit1Sim <- Xfit1%*%t(bs[, ind1])

 
 
 if( !is.null(x$X2) ){
   Xfit2 <- predict(x, eq = 2, type = "lpmatrix", ...) 
   fit2  <- predict(x, eq = 2, type = "link", ...) 
   fit2Sim <- Xfit2%*%t(bs[, ind2])

   
                     }else{
                     
   Xfit2 <- matrix(1, nrow = length(fit1), ncol = 1) 
   fit2  <- x$coefficients["sigma2.star"] 
   fit2Sim <- bs[, ind2]
   
                     
                          }

}



if(mar %in% x$VC$m3){

ind1 <- (1:x$X1.d2)
if( !is.null(x$X2) ) ind2 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2) else ind2 <- x$X1.d2 + 1     
if( !is.null(x$X3) ) ind3 <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) else ind3 <- x$X1.d2 + 2     

 Xfit1 <- predict(x, eq = 1, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 1, type = "link", ...)
 fit1Sim <- Xfit1%*%t(bs[, ind1])

 if( !is.null(x$X2) ){
 
 Xfit2 <- predict(x, eq = 2, type = "lpmatrix", ...) 
 fit2  <- predict(x, eq = 2, type = "link", ...) 
 fit2Sim <- Xfit2%*%t(bs[, ind2])
 
 Xfit3 <- predict(x, eq = 3, type = "lpmatrix", ...) 
 fit3  <- predict(x, eq = 3, type = "link", ...)  
 fit3Sim <- Xfit3%*%t(bs[, ind3])


                      }else{
                      
 Xfit2 <- matrix(1, nrow = length(fit1), ncol = 1) 
 fit2  <- x$coefficients["sigma2.star"]
 fit2Sim <- bs[, ind2]
 
 
 Xfit3 <- matrix(1, nrow = length(fit1), ncol = 1) 
 fit3  <- x$coefficients["nu.star"]                     
 fit3Sim <- bs[, ind3]
                      
                      
                      
                           }
                           
}



 


}########################################













 if(mar %in% c("SM")){
  
     if(fun == "mean"){fit    <- exp(fit1)/gamma(exp(fit3))*gamma( 1+1/sqrt(exp(fit2)) )*gamma( -1/sqrt(exp(fit2))+exp(fit3) )             
                       fitSim <- exp(fit1Sim)/gamma(exp(fit3Sim))*gamma( 1+1/sqrt(exp(fit2Sim)) )*gamma( -1/sqrt(exp(fit2Sim))+exp(fit3Sim) )                
                      } 
     
     if(fun == "variance"){
                           fit    <- exp(fit1)^2*( gamma(1+2/sqrt(exp(fit2)))*gamma(exp(fit3))*gamma(-2/sqrt(exp(fit2))+exp(fit3))-gamma(1+1/sqrt(exp(fit2)))^2*gamma(-1/sqrt(exp(fit2))+exp(fit3))^2 )
                           fitSim <- exp(fit1Sim)^2*( gamma(1+2/sqrt(exp(fit2Sim)))*gamma(exp(fit3Sim))*gamma(-2/sqrt(exp(fit2Sim))+exp(fit3Sim))-gamma(1+1/sqrt(exp(fit2Sim)))^2*gamma(-1/sqrt(exp(fit2Sim))+exp(fit3Sim))^2  )
                          }
                                                  
 } 




 if(mar %in% c("GP")){
 
    if(fun == "mean")     {fit <- exp(fit1) / (  1 - sqrt(fit2) )                          ; fitSim <- exp(fit1Sim) / (1 - sqrt(fit2Sim))                } 
    if(fun == "variance") {fit <- exp(fit1)^2/ (1 + sqrt(fit2))^2 * (1 - 2*sqrt(fit2)); fitSim <- exp(fit1Sim)^2/ (1 + sqrt(fit2Sim))^2 * (1 - 2*sqrt(fit2Sim)) }
                                     
 } # no index included here




 if(mar %in% c("DGP")){
 
    if(fun == "mean")     {fit <- 1;             fitSim <- 1                } 
    if(fun == "variance") {fit <- 1 ; fitSim <- 1 }
                                     
 } # to do, we need truncation here






 if(mar %in% c("BE")){
 
    if(fun == "mean")     {fit <- exp(fit1);             fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- exp(fit1)*(1-exp(fit1))*exp(fit2); fitSim <- exp(fit1Sim)*(1-exp(fit1Sim))*exp(fit2Sim) }
                                     
 } 
 
 
 if(mar %in% c("FISK")){
 

 
    if(fun == "mean"){
    
    sigma <- sqrt(exp(fit2))
    sigma <- ifelse(sigma <= 1, 1.0000001, sigma)    
 
    sigmaSim <- sqrt(exp(fit2Sim))
    sigmaSim <- ifelse(sigmaSim <= 1, 1.0000001, sigmaSim)   
 
    
    fit <- exp(fit1)*pi/sigma/sin(pi/sigma); fitSim <- exp(fit1Sim)*pi/sigmaSim/sin(pi/sigmaSim)                
    
    
    
    } 
    
    
    
    if(fun == "variance"){
    
    sigma <- sqrt(exp(fit2))
    sigma <- ifelse(sigma <= 2, 2.0000001, sigma)    
 
    sigmaSim <- sqrt(exp(fit2Sim))
    sigmaSim <- ifelse(sigmaSim <= 2, 2.0000001, sigmaSim)       
    
    
    fit <- exp(fit1)^2*( 2*pi/sigma/sin(2*pi/sigma)-(pi/sigma)^2/sin(pi/sigma)^2 ); fitSim <- exp(fit1Sim)^2*( 2*pi/sigmaSim/sin(2*pi/sigmaSim)-(pi/sigmaSim)^2/sin(pi/sigmaSim)^2 )  
    
    
    
    }
                                     
 }  
 




 if(mar %in% c("GU")){
 
    if(fun == "mean")     {fit <- fit1 - 0.57722*sqrt(exp(fit2));  fitSim <- fit1Sim - 0.57722*sqrt(exp(fit2Sim))            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)/6; fitSim <- pi^2*exp(fit2Sim)/6}
                                     
 } 
 
 
 
 if(mar %in% c("rGU")){
 
    if(fun == "mean")     {fit <- fit1 + 0.57722*sqrt(exp(fit2));  fitSim <- fit1Sim + 0.57722*sqrt(exp(fit2Sim))            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)/6; fitSim <- pi^2*exp(fit2Sim)/6}
                                     
 }  




 if(mar %in% c("LO")){
 
    if(fun == "mean")     {fit <- fit1;             fitSim <- fit1Sim            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)/3; fitSim <- pi^2*exp(fit2Sim)/3}
                                     
 } 
 
 
 if(mar %in% c("N")){
 
    if(fun == "mean")     {fit <- fit1;      fitSim <- fit1Sim            } 
    if(fun == "variance") {fit <- exp(fit2); fitSim <- exp(fit2Sim)       }
                                     
 } 
 
 if(mar %in% c("N2")){
 
    if(fun == "mean")     {fit <- fit1;      fitSim <- fit1Sim            } 
    if(fun == "variance") {fit <- sqrt(exp(fit2)); fitSim <- sqrt(exp(fit2Sim))       }
                                     
 }  
 
 
 

 if(mar %in% c("LN")){
 
    if(fun == "mean")     {fit <- exp(fit1)*sqrt(exp(exp(fit2)));                    fitSim <- exp(fit1Sim)*sqrt(exp(exp(fit2Sim)))                       } 
    if(fun == "variance") {fit <- exp(exp(fit2))*( exp(exp(fit2)) - 1 )*exp(2*fit1); fitSim <- exp(exp(fit2Sim))*( exp(exp(fit2Sim)) - 1 )*exp(2*fit1Sim) }
                                    
 } 
 

 if(mar %in% c("iG")){
 
    if(fun == "mean")     {fit <- exp(fit1);             fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- exp(fit1)^3*exp(fit2); fitSim <- exp(fit1Sim)^3*exp(fit2Sim) }
                                    
 }  
 
 
 
  if(mar %in% c("GA")){
  
     if(fun == "mean")     {fit <- exp(fit1);             fitSim <- exp(fit1Sim)                } 
     if(fun == "variance") {fit <- exp(fit1)^2*exp(fit2); fitSim <- exp(fit1Sim)^2*exp(fit2Sim) }
                                     
 } 
 
 
  if(mar %in% c("WEI")){
  
     if(fun == "mean")     {fit <- exp(fit1)*gamma(1+1/sqrt(exp(fit2))); fitSim <- exp(fit1Sim)*gamma(1+1/sqrt(exp(fit2Sim)))                } 
     if(fun == "variance") {fit <- exp(fit1)^2*( gamma(1+2/sqrt(exp(fit2))) - gamma( 1+1/sqrt(exp(fit2)) )^2  )  ; fitSim <-  exp(fit1Sim)^2*( gamma(1+2/sqrt(exp(fit2Sim))) - gamma( 1+1/sqrt(exp(fit2Sim)) )^2  )   }
                                     
 }  
 
 
 
 
 
  if(mar %in% c("DAGUM")){
  
     if(fun == "mean"){
     
         sigma <- sqrt(exp(fit2))
         sigma <- ifelse(sigma <= 1, 1.0000001, sigma)    
      
         sigmaSim <- sqrt(exp(fit2Sim))
    sigmaSim <- ifelse(sigmaSim <= 1, 1.0000001, sigmaSim)  
    
     
     
     fit    <- -(exp(fit1)/sigma)*gamma(-1/sigma)*gamma(1/sigma+exp(fit3))/gamma(exp(fit3))             
                       fitSim <- -(exp(fit1Sim)/sigmaSim)*gamma(-1/sigmaSim)*gamma(1/sigmaSim+exp(fit3Sim))/gamma(exp(fit3Sim))                
                      } 
     
     if(fun == "variance"){
     
     
         sigma <- sqrt(exp(fit2))
         sigma <- ifelse(sigma <= 2, 2.0000001, sigma)    
      
         sigmaSim <- sqrt(exp(fit2Sim))
    sigmaSim <- ifelse(sigmaSim <= 2, 2.0000001, sigmaSim)  
     
     
                           fit    <- -(exp(fit1)/sigma)^2*( 2*sigma*gamma(-2/sigma)*gamma(2/sigma + exp(fit3))/gamma(exp(fit3)) + ( gamma(-1/sigma)*gamma(1/sigma + exp(fit3))/gamma(exp(fit3))  )^2   )
                           fitSim <- -(exp(fit1Sim)/sigmaSim)^2*( 2*sigmaSim*gamma(-2/sigmaSim)*gamma(2/sigmaSim + exp(fit3Sim))/gamma(exp(fit3Sim)) + ( gamma(-1/sigmaSim)*gamma(1/sigmaSim + exp(fit3Sim))/gamma(exp(fit3Sim))  )^2   ) 
                          }
                           
                           
 } 
 
 
 
 
 
 if(mar %in% c("PO")){ # p. 59
 
    if(fun == "mean" || fun == "variance")     {fit <- exp(fit1); fitSim <- exp(fit1Sim) } 
                                     
 } 
 
 
 
 if(mar %in% c("NBI", "PIG")){ # p. 63, 65
 
    if(fun == "mean")     {fit <- exp(fit1);                               fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- exp(fit1) + sqrt(exp(fit2))*exp(fit1)^2; fitSim <- exp(fit1Sim) + sqrt(exp(fit2Sim))*exp(fit1Sim)^2 }
                                     
 } 
 
 
 if(mar %in% c("NBII")){ # p. 63
 
    if(fun == "mean")     {fit <- exp(fit1);                         fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- ( 1 + sqrt(exp(fit2)) )*exp(fit1); fitSim <- ( 1 + sqrt(exp(fit2Sim)) )*exp(fit1Sim) }
                                     
 }  
 
 if(mar %in% c("ZTP")){ # paper Cza
 
    if(fun == "mean")     {fit <- exp(fit1)/( 1 - exp(-exp(fit1)) ); fitSim <- exp(fit1Sim)/( 1 - exp(-exp(fit1Sim)) )                } 
    if(fun == "variance") {fit <- ( exp(fit1)*( 1 - exp(-exp(fit1))*(exp(fit1) + 1)) )/( 1 - exp(-exp(fit1)) )^2  ; fitSim <- ( exp(fit1Sim)*( 1 - exp(-exp(fit1Sim))*(exp(fit1Sim) + 1)) )/( 1 - exp(-exp(fit1Sim)) )^2 }
                                     
 }   
  
 
 
 
 
 
 
 
  
 
CIpred <- rowQuantiles(fitSim, probs = c(prob.lev/2, 1-prob.lev/2), na.rm = TRUE)


}



if(fun == "tau"){


res <- jc.probs(x, mean(x$y1), mean(x$y2), intervals = TRUE, n.sim = n.sim, prob.lev = prob.lev, ...)

fit    <- res[, 6]
CIpred <- res[, 7:8]


}


list(pred = fit, CIpred = CIpred)

}


