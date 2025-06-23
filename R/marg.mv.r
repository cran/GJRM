marg.mv <- function(x, eq, newdata, fun = "mean", n.sim = 100, prob.lev = 0.05, bin.model = NULL){


 pb <- pbS <- 1

 if( !(fun %in% c("mean", "variance")) ) stop("The value for fun can be either mean or variance.")
 
 if( !is.null(x$ordinal) && x$margins[1] %in% c("logit", "probit")  && x$margins[2] %in% c("logit", "probit") && fun %in% c("mean", "variance")  ) stop("This function is not suitable for the fitted model.")

 if( !is.null(x$ordinal) && x$margins[1] %in% c("logit", "probit")  && !(x$margins[2] %in% c("logit", "probit")) && eq == 1  ) stop("Mean and variance are available for the second margin only.")

 if(missing(newdata))                               stop("You must provide a data frame.")
 
 if(!missing(newdata)){ 
                        if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                        if(dim(newdata)[1] > 1)     stop("This calculation is currently supported for single row data frames.")
 
                      }

 
 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)
 CLM.shift <- 0
 


 if(!is.null(bin.model)){
 
 if(class(bin.model)[1] != "gam")           stop("The binary model has to be fitted using gam() from mgcv.")
 if(family(bin.model)$family != "binomial") stop("A binary model is required here.")


 Xb  <- predict(bin.model, newdata = newdata, type = "lpmatrix")
 bss <- rMVN(n.sim, mean = bin.model$coefficients, sigma = bin.model$Vp)
  
 pb  <- predict(bin.model, newdata = newdata, type = "response")
 
 if( family(bin.model)$link == "probit")  pbS <- pnorm(Xb%*%t(bss)) 
 if( family(bin.model)$link == "logit")   pbS <- plogis(Xb%*%t(bss)) 
 if( family(bin.model)$link == "cloglog") pbS <- 1-exp(-exp(Xb%*%t(bss)))
 
 } 
 
 
 
 if(x$VC$univ.gamls == FALSE){############### 
 
 if(missing(eq)) stop("You must provide the equation number.")
 if( !(eq %in% c(1, 2)) ) stop("The value for eq can be either 1 or 2.")

 if(eq == 1){ mar <- x$margins[1]; ltt <- x$left.trunc1} 
 if(eq == 2){ mar <- x$margins[2]; ltt <- x$left.trunc2} 
 
 if( mar %in% c("GP","GPII","GPo","DGP","DGPII","DGP0") ) stop("Function not ready yet for the chosen distribution(s).")
 #if( mar %in% c("logit","probit","cloglog") )             stop("Function not suitable for binary margin(s).")

 
 

if( x$margins[1] %in% c(x$VC$m2, x$VC$m2d, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m2d, x$VC$m3) ){ ############### 2/3-2/3 params
 
  
 if(eq == 1){  
              ind1 <- (1:x$X1.d2)
              
              if(!is.null(x$X3)) ind2 <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) else ind2 <- x$X1.d2 + x$X2.d2 + 1  
            
              if( x$margins[1] %in% c(x$VC$m3) ){ if(!is.null(x$X5)) ind3 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2) else ind3 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1  }
            
            }
 
 
 
 
 
 if(eq == 2){
              ind1 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
              
              if(!is.null(x$X4)) ind2 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2) else ind2 <- x$X1.d2 + x$X2.d2 + 1 + 1   
            
              if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3) ){ if(!is.null(x$X6)) ind3 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2) else ind3 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1  }

              if( x$margins[1] %in% c(x$VC$m2,x$VC$m2d) && x$margins[2] %in% c(x$VC$m3) ){ if(!is.null(x$X5)) ind3 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2) else ind3 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1  }
            
            }





if(eq == 1){

 Xfit1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
 fit1  <- predict(x, eq = 1, newdata = newdata, type = "link")
 
 if(!is.null(x$X3)){ 
 
                     Xfit2 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix") 
                     fit2  <- predict(x, eq = 3, newdata = newdata, type = "link") 
 
             if( x$margins[1] %in% c(x$VC$m3)){
             
                     Xfit3 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix") 
                     fit3  <- predict(x, eq = 5, newdata = newdata, type = "link") 
             
                                              }
                    }
                    
 if(is.null(x$X3)){ 
 
                     Xfit2 <- matrix(1, nrow = length(fit1), ncol = 1)  
                     fit2  <- x$coefficients[x$X1.d2 + x$X2.d2 + 1] 
 
             if( x$margins[1] %in% c(x$VC$m3)){
             
                     Xfit3 <- matrix(1, nrow = length(fit1), ncol = 1)  
                     fit3  <- x$coefficients[x$X1.d2 + x$X2.d2 + 1 + 1 + 1] 
             
                                              }
                    }                    
 
 
fit1Sim <- Xfit1%*%t(bs[, ind1])
fit2Sim <- Xfit2%*%t(bs[, ind2])
if( x$margins[1] %in% c(x$VC$m3)) fit3Sim <- Xfit3%*%t(bs[, ind3])


           }





if(eq == 2){

 Xfit1 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix")
 fit1  <- predict(x, eq = 2, newdata = newdata, type = "link")
 
 if(!is.null(x$X4)){ 
 
                     Xfit2 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix") 
                     fit2  <- predict(x, eq = 4, newdata = newdata, type = "link") 
 
             if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3)){
             
                     Xfit3 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix") 
                     fit3  <- predict(x, eq = 6, newdata = newdata, type = "link") 
             
                                              }

             if( x$margins[1] %in% c(x$VC$m2,x$VC$m2d) && x$margins[2] %in% c(x$VC$m3)){
             
                     Xfit3 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix") 
                     fit3  <- predict(x, eq = 5, newdata = newdata, type = "link") 
             
                                              }
                                              
                                              
                    }
                    
 if(is.null(x$X4)){ 
 
                     Xfit2 <- matrix(1, nrow = length(fit1), ncol = 1)  
                     fit2  <- x$coefficients[x$X1.d2 + x$X2.d2 + 1 + 1] 
 
             if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3)){
             
                     Xfit3 <- matrix(1, nrow = length(fit1), ncol = 1)  
                     fit3  <- x$coefficients[x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1] 
             
                                              }
                                              
             if( x$margins[1] %in% c(x$VC$m2,x$VC$m2d) && x$margins[2] %in% c(x$VC$m3)){
             
                     Xfit3 <- matrix(1, nrow = length(fit1), ncol = 1)  
                     fit3  <- x$coefficients[x$X1.d2 + x$X2.d2 + 1 + 1 + 1] 
             
                                              }                                              
                                              
                    }                    
 
 
fit1Sim <- Xfit1%*%t(bs[, ind1])
fit2Sim <- Xfit2%*%t(bs[, ind2])
if( x$margins[2] %in% c(x$VC$m3)) fit3Sim <- Xfit3%*%t(bs[, ind3])
 
 

           }











}############




if( x$margins[1] %in% c(x$VC$bl, x$VC$m1d) && x$margins[2] %in% c(x$VC$bl, x$VC$m1d, x$VC$m2, x$VC$m3, x$VC$m2d) ){ #############  1-1/2/3 params
 


  
if(!is.null(x$VC$K1)) CLM.shift <- x$VC$K1 - 1  

  
 if(eq == 1){ind1 <- ((1 + CLM.shift) : (x$X1.d2 + CLM.shift)) }     
 
 if(eq == 2){ind1 <- (x$X1.d2 + 1 + CLM.shift):(x$X1.d2 + x$X2.d2 + CLM.shift) 
 
     
    if( !(x$margins[2] %in% c(x$VC$bl, x$VC$m1d)) ){
     if(is.null(x$X3)) ind2 <- x$X1.d2 + x$X2.d2 + 1 + CLM.shift           else ind2 <- (x$X1.d2 + x$X2.d2 + 1 + CLM.shift):(x$X1.d2 + x$X2.d2 + x$X3.d2 + CLM.shift) 
     if(is.null(x$X4)) ind3 <- x$X1.d2 + x$X2.d2 + 1 + 1 + CLM.shift else       ind3 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1 + CLM.shift):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + CLM.shift)
                                           }
             
             }  

if(eq == 1){
 Xfit1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
 fit1  <- predict(x, eq = 1, newdata = newdata, type = "link")
           }

if(eq == 2){


 Xfit1 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix")
 fit1  <- predict(x, eq = 2, newdata = newdata, type = "link")
 
 
 if( !(x$margins[2] %in% c(x$VC$bl, x$VC$m1d))) { ###
 
 if(!is.null(x$X3)){
 
   Xfit2 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix") 
   fit2  <- predict(x, eq = 3, newdata = newdata, type = "link")  
   
                   }
                   
 if(!is.null(x$X4)){
 
   Xfit3 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix") 
   fit3  <- predict(x, eq = 4, newdata = newdata, type = "link")  
   
                   }                   
 
 if(is.null(x$X3)){
 
   Xfit2 <- matrix(1, nrow = length(fit1), ncol = 1) 
   fit2  <- x$coefficients[ind2]
   
                   } 
                   
 if(is.null(x$X4)){
 
   Xfit3 <- matrix(1, nrow = length(fit1), ncol = 1) 
   fit3  <- x$coefficients[ind3]
   
                   }                    
                                      } ###
           }
           
           
fit1Sim <- Xfit1%*%t(bs[, ind1])

if(eq == 2){ 

     if( !(x$margins[2] %in% c(x$VC$bl, x$VC$m1d))) {

       if(!is.null(x$X3)) fit2Sim <- Xfit2%*%t(bs[, ind2]) else fit2Sim <- bs[, ind2]
       if(!is.null(x$X4)) fit3Sim <- Xfit3%*%t(bs[, ind3]) else fit3Sim <- bs[, ind3]
       
                                           }
       
           }


}##################









}##############################################






if(x$VC$univ.gamls == TRUE){############### 


ltt <- x$left.trunc 

mar <- x$margins[1]

if( mar %in% c("GP","GPII","GPo","DGP","DGPII") ) stop("Function not ready yet for the chosen distribution.")
 
 
                     ind1 <- (1:x$X1.d2)
if( !is.null(x$X2) ) ind2 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)                    
if(  is.null(x$X2) ) ind2 <- x$X1.d2 + 1   
if( !is.null(x$X3) ) ind3 <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)
if(  is.null(x$X3) ) ind3 <- x$X1.d2 + 1 + 1    


 Xfit1   <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
 fit1    <- predict(x, eq = 1, newdata = newdata, type = "link")
 fit1Sim <- Xfit1%*%t(bs[, ind1])

 
 
 if( !is.null(x$X2) ){
   Xfit2   <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix") 
   fit2    <- predict(x, eq = 2, newdata = newdata, type = "link") 
   fit2Sim <- Xfit2%*%t(bs[, ind2])

                     }
                     
 if( is.null(x$X2) && mar %in% c(x$VC$m2, x$VC$m2d, x$VC$m3) ){                     
                     
   fit2    <- x$coefficients[ind2] 
   fit2Sim <- bs[, ind2]                     
                                                      }

 if( !is.null(x$X3) ){

   Xfit3   <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix") 
   fit3    <- predict(x, eq = 3, newdata = newdata, type = "link")  
   fit3Sim <- Xfit3%*%t(bs[, ind3]) 
 
                     }

 if( is.null(x$X3) && mar %in% c(x$VC$m3) ){

   fit3    <- x$coefficients[ind3]                     
   fit3Sim <- bs[, ind3]

                                           }
                          
} #################################################


 if(mar %in% c(x$VC$bl)){ 
 
    if(mar == "probit" ){ t.fit <- pnorm(fit1);       t.fitS <- pnorm(fit1Sim)      }  
    if(mar == "logit"  ){ t.fit <- plogis(fit1);      t.fitS <- plogis(fit1Sim)     }
    if(mar == "cloglog"){ t.fit <- 1-exp(-exp(fit1)); t.fitS <- 1-exp(-exp(fit1Sim))}
 
    if(fun == "mean")     { fit <- t.fit;           fitSim <- t.fitS             } 
    if(fun == "variance") { fit <- t.fit*(1-t.fit); fitSim <- t.fitS*(1-t.fitS)  }
                                                                         
 } 
 
 
 




 if(mar %in% c("SM")){
  
     if(fun == "mean"){fit    <- exp(fit1)/gamma(exp(fit3))*gamma( 1+1/exp(fit2) )*gamma( -1/exp(fit2)+exp(fit3) )             
                       fitSim <- exp(fit1Sim)/gamma(exp(fit3Sim))*gamma( 1+1/exp(fit2Sim) )*gamma( -1/exp(fit2Sim) + exp(fit3Sim) )   
                       
                       sn.cond  <- as.logical(exp(fit2)*exp(fit3) > 1       )
                       sn.condS <- as.logical(exp(fit2Sim)*exp(fit3Sim) > 1 )
                       
                       if(sn.cond == FALSE)                        stop("Expectation not defined.")
                       if(sum(ifelse(sn.condS == TRUE, 1, 0)) < 3) stop("Expectation not defined.")
                       
                       fitSim <- fitSim[which(sn.condS, TRUE)]

                      } 
     
     if(fun == "variance"){
                           fit    <- exp(fit1)^2*( gamma(1+2/exp(fit2))*gamma(exp(fit3))*gamma(-2/exp(fit2)+exp(fit3))-gamma(1+1/exp(fit2))^2*gamma(-1/exp(fit2)+exp(fit3))^2 )
                           fitSim <- exp(fit1Sim)^2*( gamma(1+2/exp(fit2Sim))*gamma(exp(fit3Sim))*gamma(-2/exp(fit2Sim)+exp(fit3Sim))-gamma(1+1/exp(fit2Sim))^2*gamma(-1/exp(fit2Sim)+exp(fit3Sim))^2  )
                          
                       sn.cond  <- as.logical(exp(fit2)*exp(fit3) > 2       )
                       sn.condS <- as.logical(exp(fit2Sim)*exp(fit3Sim) > 2 ) 
                       
                       if(sn.cond == FALSE)                        stop("Variance not defined.")
                       if(sum(ifelse(sn.condS == TRUE, 1, 0)) < 3) stop("Variance not defined.")
                       
                       fitSim <- fitSim[which(sn.condS, TRUE)]                        
                          
                          
                          }
                                                  
 } 



 if(mar %in% c("BE")){
 
    if(fun == "mean")     {fit <- plogis(fit1);             fitSim <- plogis(fit1Sim)                } 
    if(fun == "variance") {fit <- plogis(fit1)*(1-plogis(fit1))*plogis(fit2)^2; fitSim <- plogis(fit1Sim)*(1-plogis(fit1Sim))*plogis(fit2Sim)^2 }
                                     
 } 
 
 
 if(mar %in% c("FISK")){
 

 
    if(fun == "mean"){
    
    sigma    <- exp(fit2)
    sigmaSim <- exp(fit2Sim)
 
 
    fit    <- exp(fit1)*pi/sigma/sin(pi/sigma)
    fitSim <- exp(fit1Sim)*pi/sigmaSim/sin(pi/sigmaSim)                
 
 
 
    sn.cond  <- as.logical(sigma    > 1 )
    sn.condS <- as.logical(sigmaSim > 1 )
                        
    if(sn.cond == FALSE)                        stop("Expectation not defined.")
    if(sum(ifelse(sn.condS == TRUE, 1, 0)) < 3) stop("Expectation not defined.")
                        
    fitSim <- fitSim[which(sn.condS, TRUE)]
 

    } 
    
    
    
    if(fun == "variance"){
    
    sigma    <- exp(fit2)
    sigmaSim <- exp(fit2Sim)
    
    
    fit <- exp(fit1)^2*( 2*pi/sigma/sin(2*pi/sigma)-(pi/sigma)^2/sin(pi/sigma)^2 )
    fitSim <- exp(fit1Sim)^2*( 2*pi/sigmaSim/sin(2*pi/sigmaSim)-(pi/sigmaSim)^2/sin(pi/sigmaSim)^2 )  
    
    
    sn.cond  <- as.logical(sigma    > 2 )
    sn.condS <- as.logical(sigmaSim > 2 ) 
   
    if(sn.cond == FALSE)                        stop("Variance not defined.")
    if(sum(ifelse(sn.condS == TRUE, 1, 0)) < 3) stop("Variance not defined.")
   
    fitSim <- fitSim[which(sn.condS, TRUE)]      
    
    
    
    }
                                     
 }  
 




 if(mar %in% c("GU")){
 
    if(fun == "mean")     {fit <- fit1 - 0.57722*exp(fit2);  fitSim <- fit1Sim - 0.57722*exp(fit2Sim)            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)^2/6; fitSim <- pi^2*exp(fit2Sim)^2/6}
                                     
 } 
 
 
 
 if(mar %in% c("rGU")){
 
    if(fun == "mean")     {fit <- fit1 + 0.57722*exp(fit2);  fitSim <- fit1Sim + 0.57722*exp(fit2Sim)            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)^2/6; fitSim <- pi^2*exp(fit2Sim)^2/6}
                                     
 }  




 if(mar %in% c("LO")){
 
    if(fun == "mean")     {fit <- fit1;             fitSim <- fit1Sim            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)^2/3; fitSim <- pi^2*exp(fit2Sim)^2/3}
                                     
 } 
 
 
 if(mar %in% c("N")){
 
    if(fun == "mean")     {fit <- fit1;      fitSim <- fit1Sim            } 
    if(fun == "variance") {fit <- exp(fit2)^2; fitSim <- exp(fit2Sim)^2       }
                                     
 } 
 
 
 
 

 if(mar %in% c("LN")){ # fit1 is mu taking values in R, table in M&R CSDA contains a mistake
 
    if(fun == "mean")     {fit <- exp(fit1)*sqrt(exp(exp(fit2)^2));                      fitSim <- exp(fit1Sim)*sqrt(exp(exp(fit2Sim)^2))                         } 
    if(fun == "variance") {fit <- exp(exp(fit2)^2)*( exp(exp(fit2)^2) - 1 )*exp(2*fit1); fitSim <- exp(exp(fit2Sim)^2)*( exp(exp(fit2Sim)^2) - 1 )*exp(2*fit1Sim) }
                                    
 } 
 

 if(mar %in% c("IG")){
 
    if(fun == "mean")     {fit <- exp(fit1);             fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- exp(fit1)^3*exp(fit2)^2; fitSim <- exp(fit1Sim)^3*exp(fit2Sim)^2 }
                                    
 }  
 
 
 
  if(mar %in% c("GA")){
  
     if(fun == "mean")     {fit <- exp(fit1);             fitSim <- exp(fit1Sim)                } 
     if(fun == "variance") {fit <- exp(fit1)^2*exp(fit2)^2; fitSim <- exp(fit1Sim)^2*exp(fit2Sim)^2 }
                                     
 } 
 
 
  if(mar %in% c("WEI")){
  
     if(fun == "mean")     {fit <- exp(fit1)*gamma(1+1/exp(fit2)); fitSim <- exp(fit1Sim)*gamma(1+1/exp(fit2Sim))                } 
     if(fun == "variance") {fit <- exp(fit1)^2*( gamma(1+2/exp(fit2)) - gamma( 1+1/exp(fit2) )^2  )  ; fitSim <-  exp(fit1Sim)^2*( gamma(1+2/exp(fit2Sim)) - gamma( 1+1/exp(fit2Sim) )^2  )   }
                                     
 }  
 
 
 
 
 
  if(mar %in% c("DAGUM")){
  
     if(fun == "mean"){
     
         sigma    <- exp(fit2)
         sigmaSim <- exp(fit2Sim)
    
     
     
         fit      <- -(exp(fit1)/sigma)*gamma(-1/sigma)*gamma(1/sigma+exp(fit3))/gamma(exp(fit3))             
         fitSim   <- -(exp(fit1Sim)/sigmaSim)*gamma(-1/sigmaSim)*gamma(1/sigmaSim+exp(fit3Sim))/gamma(exp(fit3Sim))    
         
 
    sn.cond  <- as.logical(sigma    > 1 )
    sn.condS <- as.logical(sigmaSim > 1 )
                        
    if(sn.cond == FALSE)                        stop("Expectation not defined.")
    if(sum(ifelse(sn.condS == TRUE, 1, 0)) < 3) stop("Expectation not defined.")
                        
    fitSim <- fitSim[which(sn.condS, TRUE)]         
                       
                       
                       
                      } 
     
     if(fun == "variance"){
     
     
         sigma    <- exp(fit2)
         sigmaSim <- exp(fit2Sim)
     
     
         fit    <- -(exp(fit1)/sigma)^2*( 2*sigma*gamma(-2/sigma)*gamma(2/sigma + exp(fit3))/gamma(exp(fit3)) + ( gamma(-1/sigma)*gamma(1/sigma + exp(fit3))/gamma(exp(fit3))  )^2   )
         fitSim <- -(exp(fit1Sim)/sigmaSim)^2*( 2*sigmaSim*gamma(-2/sigmaSim)*gamma(2/sigmaSim + exp(fit3Sim))/gamma(exp(fit3Sim)) + ( gamma(-1/sigmaSim)*gamma(1/sigmaSim + exp(fit3Sim))/gamma(exp(fit3Sim))  )^2   ) 
         
         
    sn.cond  <- as.logical(sigma    > 2 )
    sn.condS <- as.logical(sigmaSim > 2 ) 
   
    if(sn.cond == FALSE)                        stop("Variance not defined.")
    if(sum(ifelse(sn.condS == TRUE, 1, 0)) < 3) stop("Variance not defined.")
   
    fitSim <- fitSim[which(sn.condS, TRUE)]           
          
                          
                          
                          }
                           
                           
 } 
 
 
 
 
 
 if(mar %in% c("P")){ # p. 59
 
    if(fun == "mean" || fun == "variance")     {fit <- exp(fit1); fitSim <- exp(fit1Sim) } 
                                     
 } 
 
 
 
 if(mar %in% c("NBI", "PIG")){ # p. 63, 65
 
    if(fun == "mean")     {fit <- exp(fit1);                               fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- exp(fit1) + exp(fit2)*exp(fit1)^2; fitSim <- exp(fit1Sim) + exp(fit2Sim)*exp(fit1Sim)^2 }
                                     
 } 
 
 
 if(mar %in% c("NBII")){ # p. 63
 
    if(fun == "mean")     {fit <- exp(fit1);                         fitSim <- exp(fit1Sim)                } 
    if(fun == "variance") {fit <- ( 1 + exp(fit2) )*exp(fit1); fitSim <- ( 1 + exp(fit2Sim) )*exp(fit1Sim) }
                                     
 }
 
 
  if(mar %in% c("tNBII")){ 
  
     if(fun == "mean")     {fit <- exp(fit1)/(1 - pNBII(rep(ltt, length(fit1)), mu = exp(fit1), sigma = exp(fit2)));  fitSim <- exp(fit1Sim)/(1 - pNBII(rep(ltt, length(fit1Sim)), mu = exp(fit1Sim), sigma = exp(fit2Sim)))                } 
     
      if(fun == "variance") stop("Variance not implemented yet for this distribution. \nGet in touch to check progress.")
                                      
 } 
 
 
   if(mar %in% c("tNBI")){ 
   
      if(fun == "mean")     {fit <- exp(fit1)/(1 - pNBI(rep(ltt, length(fit1)), mu = exp(fit1), sigma = exp(fit2)));  fitSim <- exp(fit1Sim)/(1 - pNBI(rep(ltt, length(fit1Sim)), mu = exp(fit1Sim), sigma = exp(fit2Sim)))                } 
      
      if(fun == "variance") stop("Variance not implemented yet for this distribution. \nGet in touch to check progress.")
                                       
 } 
 
   if(mar %in% c("tPIG")){ 
   
      if(fun == "mean")     {fit <- exp(fit1)/(1 - pPIG(rep(ltt, length(fit1)), mu = exp(fit1), sigma = exp(fit2)));  fitSim <- exp(fit1Sim)/(1 - pPIG(rep(ltt, length(fit1Sim)), mu = exp(fit1Sim), sigma = exp(fit2Sim)))                } 
      
      if(fun == "variance") stop("Variance not implemented yet for this distribution. \nGet in touch to check progress.")
                                       
 }  
 
 
 
 if(mar %in% c("tP")){ # 
 
    if(fun == "mean")     {fit <- exp(fit1)/ (1 - pPO(rep(ltt, length(fit1)), exp(fit1))); fitSim <- exp(fit1Sim)/ (1 - pPO(rep(ltt, length(fit1Sim)), exp(fit1Sim)))                } 
    #if(fun == "variance") {fit <- (exp(fit1) + exp(fit1)^2)/( 1 - exp(-exp(fit1)) ) - exp(fit1)^2/( 1 - exp(-exp(fit1)) )^2   ; fitSim <- (exp(fit1Sim) + exp(fit1Sim)^2)/( 1 - exp(-exp(fit1Sim)) ) - exp(fit1Sim)^2/( 1 - exp(-exp(fit1Sim)) )^2 }
    if(fun == "variance") stop("Variance not implemented yet for this distribution. \nGet in touch to check progress.")

 
 }   
 
 
 

  if(mar %in% c("tN")){
  
  lttr <- rep(ltt, length(fit1))  
  alpha  <- (lttr - fit1) / exp(fit2) 
  alphaS <- (lttr - fit1Sim) / exp(fit2Sim)
  

     if(fun == "mean"){
     
             fit    <- fit1    + exp(fit2)*(dnorm(alpha)/(1 - pnorm(alpha)))
             fitSim <- fit1Sim + exp(fit2Sim)*(dnorm(alphaS)/(1 - pnorm(alphaS)))                
             
             } 
     
     if(fun == "variance"){
     
             fit    <- exp(fit2)^2 * (1 + (alpha * dnorm(alpha)) / (1 - pnorm(alpha)) - (dnorm(alpha) / (1 - pnorm(alpha)))^2) 
             fitSim <- exp(fit2Sim)^2 * (1 + (alphaS * dnorm(alphaS)) / (1 - pnorm(alphaS)) - (dnorm(alphaS) / (1 - pnorm(alphaS)))^2)   
             
             }
                                     
 }  
 
 
 
# if(mar %in% c("tP")){ # 
# 
#    if(fun == "mean")     {fit <- exp(fit1)/( 1 - exp(-exp(fit1)) ); fitSim <- exp(fit1Sim)/( 1 - exp(-exp(fit1Sim)) )                } 
#    if(fun == "variance") {fit <- (exp(fit1) + exp(fit1)^2)/( 1 - exp(-exp(fit1)) ) - exp(fit1)^2/( 1 - exp(-exp(fit1)) )^2   ; fitSim <- (exp(fit1Sim) + exp(fit1Sim)^2)/( 1 - exp(-exp(fit1Sim)) ) - exp(fit1Sim)^2/( 1 - exp(-exp(fit1Sim)) )^2 }
#                                     
# }  
 
  
 
 
 if(mar %in% c("TW")){ 
 
    if(fun == "mean"){fit <- exp(fit1); fitSim <- exp(fit1Sim) } 
    
    
    if(fun == "variance"){
 
 
       aTW <- 1.001  
       bTW <- 1.999
  
       nu  <- (aTW + bTW*exp(fit3))/(1 + exp(fit3)) 
       nuS <- (aTW + bTW*exp(fit3Sim))/(1 + exp(fit3Sim)) 
   
 
                      sigma <- exp(fit2); sigmaS <- exp(fit2Sim)
 
                      fit <- sigma*exp(fit1)^nu
                      fitSim <- sigmaS*exp(fit1Sim)^nuS  
                              
    
                          }
                                 
 }  
 
 
 
lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "") 
 
 
 res.CI <- as.numeric(quantile(pbS*fitSim, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))
 
 res <- c(res.CI[1], pb*fit, res.CI[2])
 
  if(fun == "mean")     names(res) <- c(lbn, "mean",     ubn)
  if(fun == "variance") names(res) <- c(lbn, "variance", ubn) 
  


out <- list(res = res, sim.mv = pbS*fitSim, prob.lev = prob.lev, fun = fun)
 							 
class(out) <- "marg.mv"

out



}


