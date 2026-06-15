jc.probs9 <- function(x, y1, y2, y3, newdata, intervals, n.sim, prob.lev, min.pr, max.pr, tau.res = FALSE, method = "adaptive"){

res.vCI <- res.mCI <- eta.eq1.3 <- eta.eq1.3S <- eta.eq2.3 <- eta.eq2.3S <- eta.eq3.3 <- eta.eq3.3S <- theta3 <- theta3S <- NULL
eta.eq1.3S <- res.relerr <- res.niter <- nConst <- nConstTemp <- NULL
sI <- NULL

res.mCI  <- matrix(0, 1, 2) 

if (!(method %in% c("adaptive", "quad"))) stop("Invalid 'method'. Must be one of 'adaptive' or 'quad'.")


if(is.null(y1) || is.null(y2) || is.null(y3) ) stop("You must provide values for y1, y2 and y3.")

if( !(y1 >= range(x$m1$y1)[1] && y1 <= range(x$m1$y1)[2]) ) stop("A value for y1 within its observed range must be provided.")
if( !(y2 >= range(x$m2$y1)[1] && y2 <= range(x$m2$y1)[2]) ) stop("A value for y2 within its observed range must be provided.")
if( !(y3 >= range(x$m3$y1)[1] && y3 <= range(x$m3$y1)[2]) ) stop("A value for y3 within its observed range must be provided.")
    

eq.theta <- x$eq.theta
#eq <- x$cond.var   # this is the cond in the factorization, in cond mean function is instead the var I am conditioning on
cop.sel <- x$copula.sel

if(x$cond.var == 1) nC <- x$b23.1$nC
if(x$cond.var == 2) nC <- x$b13.2$nC
if(x$cond.var == 3) nC <- x$b12.3$nC


if(missing(newdata))                               stop("You must provide a data frame.")
if(!missing(newdata)){ 
                       if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                       if(dim(newdata)[1] > 1)     stop("This calculation is supported for single row data frames.")
                     }

bl  <- x$m1$VC$bl
m2  <- x$m1$VC$m2
m3  <- x$m1$VC$m3
m1d <- c("P", "tP")
m2d <- c("NBI", "NBII", "PIG","tNBI", "tNBII", "tPIG")

lt  <- x$left.trunc.v  
margins <- x$margins

if (!all(margins %in% c(m2, m3))) stop("Only continuous marginal distributions currently allowed for.")
if (any(margins %in% c("GP", "GPII", "GPo", "DGP", "DGPII", "DGP0"))) stop("Function not ready yet for the chosen distribution(s).")


y.vec <- c(y1, y2, y3)
uB    <- y.vec[x$cond.var]

if( margins[x$cond.var] %in% c("N","GU","rGU","LO") ){                             lB <- -Inf;                    }
if( margins[x$cond.var] %in% c("tN") ){                                            lB <- lt[x$cond.var];                  }
if( margins[x$cond.var] %in% c("LN","WEI","IG","GA","DAGUM","SM","FISK","TW")  ) { lB <- sqrt(.Machine$double.eps)}
if( margins[x$cond.var] %in% c("BE")  ) {                                          lB <- sqrt(.Machine$double.eps)}




############################      

 Xeq1.1 <- predict(x$m1, eq = 1, newdata = newdata, type = "lpmatrix")
 if(!is.null(x$m1$X2)){
                                                 Xeq1.2 <- predict(x$m1, eq = 2, newdata = newdata, type = "lpmatrix")
                          if(margins[1] %in% m3) Xeq1.3 <- predict(x$m1, eq = 3, newdata = newdata, type = "lpmatrix")
                        }else{Xeq1.2 <- Xeq1.3 <- matrix(1, dim(Xeq1.1)[1], 1)}
 

 Xeq2.1 <- predict(x$m2, eq = 1, newdata = newdata, type = "lpmatrix")
 if(!is.null(x$m2$X2)){
                                                 Xeq2.2 <- predict(x$m2, eq = 2, newdata = newdata, type = "lpmatrix")
                          if(margins[2] %in% m3) Xeq2.3 <- predict(x$m2, eq = 3, newdata = newdata, type = "lpmatrix")
                        }else{Xeq2.2 <- Xeq2.3 <- matrix(1, dim(Xeq2.1)[1], 1)}
 

 Xeq3.1 <- predict(x$m3, eq = 1, newdata = newdata, type = "lpmatrix")
 if(!is.null(x$m3$X2)){
                                                 Xeq3.2 <- predict(x$m3, eq = 2, newdata = newdata, type = "lpmatrix")
                          if(margins[3] %in% m3) Xeq3.3 <- predict(x$m3, eq = 3, newdata = newdata, type = "lpmatrix")
                        }else{Xeq3.2 <- Xeq3.3 <- matrix(1, dim(Xeq3.1)[1], 1)}
  
############################
 
 bs1 <- rMVN(n.sim, mean = x$m1$coefficients, sigma = x$m1$Vb)
 
 eta.eq1.1  <- c(Xeq1.1 %*% x$m1$coefficients[1:x$m1$X1.d2]) 
 eta.eq1.1S <- c(Xeq1.1 %*% t(bs1[, 1:x$m1$X1.d2])) 
 
 if(!is.null(x$m1$X2)){   ind.t      <- (x$m1$X1.d2 + 1):(x$m1$X1.d2 + x$m1$X2.d2) 
                          eta.eq1.2  <- c(esp.tr(Xeq1.2 %*% x$m1$coefficients[ind.t], margins[1])$vrb) # these are transformed already to be used in the distrAT
                          eta.eq1.2S <- c(esp.tr(Xeq1.2 %*% t(bs1[, ind.t]), margins[1])$vrb)
                          
                          if(margins[1] %in% m3 ){
                          
                            ind.t      <- (x$m1$X1.d2 + x$m1$X2.d2 + 1):(x$m1$X1.d2 + x$m1$X2.d2 + x$m1$X3.d2)
                            eta.eq1.3  <- c(enu.tr(Xeq1.3 %*% x$m1$coefficients[ind.t], margins[1])$vrb) 
                            eta.eq1.3S <- c(enu.tr(Xeq1.3 %*% t(bs1[, ind.t]), margins[1])$vrb) 
                          
                                                 }

                        }else{
                        
                              eta.eq1.2  <- c(esp.tr(x$m1$coefficients[x$m1$X1.d2 + 1], margins[1])$vrb)   
                              eta.eq1.2S <- c(esp.tr(bs1[, x$m1$X1.d2 + 1], margins[1])$vrb) 
                              
                              
                              if(margins[1] %in% m3 ){
                              eta.eq1.3  <- c(enu.tr(x$m1$coefficients[x$m1$X1.d2 + 1 + 1], margins[1])$vrb)
                              eta.eq1.3S <- c(enu.tr(bs1[, x$m1$X1.d2 + 1 + 1], margins[1])$vrb)                              
                                                      }
                             }
 
 ##
 
  bs2 <- rMVN(n.sim, mean = x$m2$coefficients, sigma = x$m2$Vb)
 
  eta.eq2.1  <- c(Xeq2.1 %*% x$m2$coefficients[1:x$m2$X1.d2]) 
  eta.eq2.1S <- c(Xeq2.1 %*% t(bs2[, 1:x$m2$X1.d2])) 
  
  if(!is.null(x$m2$X2)){
                           ind.t      <- (x$m2$X1.d2 + 1):(x$m2$X1.d2 + x$m2$X2.d2)
                           eta.eq2.2  <- c(esp.tr(Xeq2.2 %*% x$m2$coefficients[ind.t], margins[2])$vrb) 
                           eta.eq2.2S <- c(esp.tr(Xeq2.2 %*% t(bs2[, ind.t]) , margins[2])$vrb) 
                           
                           if(margins[2] %in% m3 ){
                           
                             ind.t      <- (x$m2$X1.d2 + x$m2$X2.d2 + 1):(x$m2$X1.d2 + x$m2$X2.d2 + x$m2$X3.d2)
                             eta.eq2.3  <- c(enu.tr(Xeq2.3 %*% x$m2$coefficients[ind.t] , margins[2])$vrb)
                             eta.eq2.3S <- c(enu.tr(Xeq2.3 %*% t(bs2[, ind.t]) , margins[2])$vrb)
                           
                                                  }
  
  
                         }else{
                         
                               eta.eq2.2  <- c(esp.tr(x$m2$coefficients[x$m2$X1.d2 + 1] , margins[2])$vrb)
                               eta.eq2.2S <- c(esp.tr(bs2[, x$m2$X1.d2 + 1], margins[2])$vrb)   
                               
                               if(margins[2] %in% m3 ){

                               eta.eq2.3  <- c(enu.tr(x$m2$coefficients[x$m2$X1.d2 + 1 + 1], margins[2])$vrb)
                               eta.eq2.3S <- c(enu.tr(bs2[, x$m2$X1.d2 + 1 + 1], margins[2])$vrb)      
                               
                               }
                         
                              }
  
##

  bs3 <- rMVN(n.sim, mean = x$m3$coefficients, sigma = x$m3$Vb)
 
  eta.eq3.1  <- c(Xeq3.1 %*% x$m3$coefficients[1:x$m3$X1.d2]) 
  eta.eq3.1S <- c(Xeq3.1 %*% t(bs3[, 1:x$m3$X1.d2])) 
  
  if(!is.null(x$m3$X2)){
                           ind.t      <- (x$m3$X1.d2 + 1):(x$m3$X1.d2 + x$m3$X2.d2)
                           eta.eq3.2  <- c(esp.tr(Xeq3.2 %*% x$m3$coefficients[ind.t] , margins[3])$vrb)
                           eta.eq3.2S <- c(esp.tr(Xeq3.2 %*% t(bs3[, ind.t]), margins[3])$vrb) 
                           
                           if(margins[3] %in% m3 ){
                           
                             ind.t      <- (x$m3$X1.d2 + x$m3$X2.d2 + 1):(x$m3$X1.d2 + x$m3$X2.d2 + x$m3$X3.d2)
                             eta.eq3.3  <- c(enu.tr(Xeq3.3 %*% x$m3$coefficients[ind.t], margins[3])$vrb) 
                             eta.eq3.3S <- c(enu.tr(Xeq3.3 %*% t(bs3[, ind.t]) , margins[3])$vrb) 
                           
                                                  }
  
  
                         }else{
                         
                               eta.eq3.2  <- c(esp.tr(x$m3$coefficients[x$m3$X1.d2 + 1] , margins[3])$vrb)  
                               eta.eq3.2S <- c(esp.tr(bs3[, x$m3$X1.d2 + 1], margins[3])$vrb)  

                               if(margins[3] %in% m3 ){

                               eta.eq3.3  <- c(enu.tr(x$m3$coefficients[x$m3$X1.d2 + 1 + 1], margins[3])$vrb)
                               eta.eq3.3S <- c(enu.tr(bs3[, x$m3$X1.d2 + 1 + 1], margins[3])$vrb)
                               
                                                      }
                         
                              }
   
 
 
##



if(x$cond.var == 1){  
                            
  bs12   <- rMVN(n.sim, mean = x$b12$coefficients,   sigma = x$b12$Vb)
  bs13   <- rMVN(n.sim, mean = x$b13$coefficients,   sigma = x$b13$Vb)
  bs23.1 <- rMVN(n.sim, mean = x$b23.1$coefficients, sigma = x$b23.1$Vb)
  
  if(!is.null(x$b12$X5))   Xt1 <- predict(x$b12,   newdata = newdata, type = "lpmatrix") else Xt1 <- matrix(1, dim(Xeq1.1)[1], 1)
  if(!is.null(x$b13$X5))   Xt2 <- predict(x$b13,   newdata = newdata, type = "lpmatrix") else Xt2 <- matrix(1, dim(Xeq1.1)[1], 1)
  if(x$simplified == TRUE){ if(!is.null(x$b23.1$X5)) Xt3 <- predict(x$b23.1, newdata = newdata, type = "lpmatrix") else Xt3 <- matrix(1, dim(Xeq1.1)[1], 1) }

  
  eq.theta1  <- Xt1 %*% x$b12$coefficients
  eq.theta1S <- Xt1 %*% t(bs12)
  theta1     <- c(teta.tr( x$b12$VC, eq.theta1)$teta)
  theta1S    <- c(teta.tr( x$b12$VC, eq.theta1S)$teta)

  if(x$b12$copula %in% x$b12$BivD2){

    theta1S <- sign(theta1)*abs(theta1S)
    cop.sel[1] <- ifelse(theta1 >= 0, x$b12$fit$Cop1, x$b12$fit$Cop2)
    
                                       } 


  eq.theta2  <- Xt2 %*% x$b13$coefficients
  eq.theta2S <- Xt2 %*% t(bs13)
  theta2     <- c(teta.tr( x$b13$VC, eq.theta2)$teta)
  theta2S    <- c(teta.tr( x$b13$VC, eq.theta2S)$teta)


  if(x$b13$copula %in% x$b13$BivD2){

    theta2S <- sign(theta2)*abs(theta2S)
    cop.sel[2] <- ifelse(theta2 >= 0, x$b13$fit$Cop1, x$b13$fit$Cop2)
    
                                       } 



  if(x$simplified == TRUE && is.null(x$b23.1$X5)){
  eq.theta3  <- Xt3 %*% x$b23.1$coefficients
  eq.theta3S <- Xt3 %*% t(bs23.1)
  theta3     <- c(teta.tr( x$b23.1$VC, eq.theta3)$teta)
  theta3S    <- c(teta.tr( x$b23.1$VC, eq.theta3S)$teta)
  
  if(x$b23.1$copula %in% x$b23.1$BivD2){

    nC <- ifelse( theta3 >= 0, x$b23.1$fit$nC1, x$b23.1$fit$nC2) 
    theta3S <- sign(theta3)*abs(theta3S)
    cop.sel[3] <- ifelse(theta3 >= 0, x$b23.1$fit$Cop1, x$b23.1$fit$Cop2)
    
                                       }  
                                                  }
                                  
                 } 



if(x$cond.var == 2){      
                            
  bs12   <- rMVN(n.sim, mean = x$b12$coefficients,   sigma = x$b12$Vb)
  bs23   <- rMVN(n.sim, mean = x$b23$coefficients,   sigma = x$b23$Vb)
  bs13.2 <- rMVN(n.sim, mean = x$b13.2$coefficients, sigma = x$b13.2$Vb)
  
  if(!is.null(x$b12$X5))   Xt1 <- predict(x$b12,   newdata = newdata, type = "lpmatrix") else Xt1 <- matrix(1, dim(Xeq1.1)[1], 1)
  if(!is.null(x$b23$X5))   Xt2 <- predict(x$b23,   newdata = newdata, type = "lpmatrix") else Xt2 <- matrix(1, dim(Xeq1.1)[1], 1)
  if(x$simplified == TRUE){ if(!is.null(x$b13.2$X5)) Xt3 <- predict(x$b13.2, newdata = newdata, type = "lpmatrix") else Xt3 <- matrix(1, dim(Xeq1.1)[1], 1) }
  
  # Xt3 is useless above but left for symmetry when reading the code
  
  eq.theta1  <- Xt1 %*% x$b12$coefficients
  eq.theta1S <- Xt1 %*% t(bs12)
  theta1     <- c(teta.tr( x$b12$VC, eq.theta1)$teta)
  theta1S    <- c(teta.tr( x$b12$VC, eq.theta1S)$teta)

 if(x$b12$copula %in% x$b12$BivD2){

    theta1S <- sign(theta1)*abs(theta1S)
    cop.sel[1] <- ifelse(theta1 >= 0, x$b12$fit$Cop1, x$b12$fit$Cop2)
    
                                       } 

  eq.theta2  <- Xt2 %*% x$b23$coefficients
  eq.theta2S <- Xt2 %*% t(bs23)
  theta2     <- c(teta.tr( x$b23$VC, eq.theta2)$teta)
  theta2S    <- c(teta.tr( x$b23$VC, eq.theta2S)$teta)

 if(x$b23$copula %in% x$b23$BivD2){

    theta2S <- sign(theta2)*abs(theta2S)
    cop.sel[2] <- ifelse(theta2 >= 0, x$b23$fit$Cop1, x$b23$fit$Cop2)
    
                                       } 



  if(x$simplified == TRUE && is.null(x$b13.2$X5)){
  eq.theta3  <- Xt3 %*% x$b13.2$coefficients
  eq.theta3S <- Xt3 %*% t(bs13.2)
  theta3     <- c(teta.tr( x$b13.2$VC, eq.theta3)$teta)
  theta3S    <- c(teta.tr( x$b13.2$VC, eq.theta3S)$teta)
  
    if(x$b13.2$copula %in% x$b13.2$BivD2){
  
      nC <- ifelse( theta3 >= 0, x$b13.2$fit$nC1, x$b13.2$fit$nC2) 
      theta3S <- sign(theta3)*abs(theta3S)
      cop.sel[3] <- ifelse(theta3 >= 0, x$b13.2$fit$Cop1, x$b13.2$fit$Cop2)
      
                                       }
  
  
                                                 }
  
                 } 



if(x$cond.var == 3){     
                         
  bs13   <- rMVN(n.sim, mean = x$b13$coefficients,   sigma = x$b13$Vb)
  bs23   <- rMVN(n.sim, mean = x$b23$coefficients,   sigma = x$b23$Vb)
  bs12.3 <- rMVN(n.sim, mean = x$b12.3$coefficients, sigma = x$b12.3$Vb)
  
  if(!is.null(x$b13$X5))   Xt1 <- predict(x$b13,   newdata = newdata, type = "lpmatrix") else Xt1 <- matrix(1, dim(Xeq1.1)[1], 1)
  if(!is.null(x$b23$X5))   Xt2 <- predict(x$b23,   newdata = newdata, type = "lpmatrix") else Xt2 <- matrix(1, dim(Xeq1.1)[1], 1)
  if(x$simplified == TRUE){ if(!is.null(x$b12.3$X5)) Xt3 <- predict(x$b12.3, newdata = newdata, type = "lpmatrix") else Xt3 <- matrix(1, dim(Xeq1.1)[1], 1) }
  
  eq.theta1  <- Xt1 %*% x$b13$coefficients
  eq.theta1S <- Xt1 %*% t(bs13)
  theta1     <- c(teta.tr( x$b13$VC, eq.theta1)$teta)
  theta1S    <- c(teta.tr( x$b13$VC, eq.theta1S)$teta)

 if(x$b13$copula %in% x$b13$BivD2){

    theta1S <- sign(theta1)*abs(theta1S)
    cop.sel[1] <- ifelse(theta1 >= 0, x$b13$fit$Cop1, x$b13$fit$Cop2)
    
                                       } 


  eq.theta2  <- Xt2 %*% x$b23$coefficients
  eq.theta2S <- Xt2 %*% t(bs23)
  theta2     <- c(teta.tr( x$b23$VC, eq.theta2)$teta)
  theta2S    <- c(teta.tr( x$b23$VC, eq.theta2S)$teta)

 if(x$b23$copula %in% x$b23$BivD2){

    theta2S <- sign(theta2)*abs(theta2S)
    cop.sel[2] <- ifelse(theta2 >= 0, x$b23$fit$Cop1, x$b23$fit$Cop2)
    
                                       } 



  if(x$simplified == TRUE && is.null(x$b12.3$X5)){
  eq.theta3  <- Xt3 %*% x$b12.3$coefficients
  eq.theta3S <- Xt3 %*% t(bs12.3)
  theta3     <- c(teta.tr( x$b12.3$VC, eq.theta3)$teta)
  theta3S    <- c(teta.tr( x$b12.3$VC, eq.theta3S)$teta)
  
    if(x$b12.3$copula %in% x$b12.3$BivD2){
  
      nC <- ifelse( theta3 >= 0, x$b12.3$fit$nC1, x$b12.3$fit$nC2) 
      theta3S <- sign(theta3)*abs(theta3S)
      cop.sel[3] <- ifelse(theta3 >= 0, x$b12.3$fit$Cop1, x$b12.3$fit$Cop2)
      
                                       }  
  
                                                  }
  
                 } 


#################################################################


if(x$cond.var == 1){ beta.theta3 <- x$b23.1$coefficients; BS3 <- bs23.1; obj.cB = x$b23.1}  
if(x$cond.var == 2){ beta.theta3 <- x$b13.2$coefficients; BS3 <- bs13.2; obj.cB = x$b13.2}
if(x$cond.var == 3){ beta.theta3 <- x$b12.3$coefficients; BS3 <- bs12.3; obj.cB = x$b12.3}

if(obj.cB$copula %in% obj.cB$BivD2){

	if( is.null(theta3) ){
   
          if(x$simplified == FALSE){ 
          
			if(x$cond.var == 1) p1 <- distrHsAT(y1, eta2 = eta.eq1.1, sigma2 = eta.eq1.2, nu = eta.eq1.3, margin2 = margins[1], min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, left.trunc = lt[1])$p2
			if(x$cond.var == 2) p2 <- distrHsAT(y2, eta2 = eta.eq2.1, sigma2 = eta.eq2.2, nu = eta.eq2.3, margin2 = margins[2], min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, left.trunc = lt[2])$p2
			if(x$cond.var == 3) p3 <- distrHsAT(y3, eta2 = eta.eq3.1, sigma2 = eta.eq3.2, nu = eta.eq3.3, margin2 = margins[3], min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, left.trunc = lt[3])$p2
          
                        if(x$cond.var == 1) Xpr.t <- predict(obj.cB, newdata = merge(data.frame(.gjrm_p1 = p1), newdata), type = "lpmatrix" ) 
                        if(x$cond.var == 2) Xpr.t <- predict(obj.cB, newdata = merge(data.frame(.gjrm_p2 = p2), newdata), type = "lpmatrix" ) 
                        if(x$cond.var == 3) Xpr.t <- predict(obj.cB, newdata = merge(data.frame(.gjrm_p3 = p3), newdata), type = "lpmatrix" ) 
                       
                                   }
                                   
          if(x$simplified == TRUE) Xpr.t <- predict(obj.cB, newdata = newdata, type = "lpmatrix" ) 
          
          eta.theta3T <- Xpr.t %*% beta.theta3
          theta3T     <- c(teta.tr(obj.cB$VC, eta.theta3T)$teta)  
          
    	  sI <- ifelse(theta3T < 0, -1, 1) # (theta3T >= 0) 1    
    	
	                      }
	
                                   }




              
if(method == "adaptive"){ 

  res.m <- distrExIntegrate(cdf.pcc, lower = lB, upper = uB, y1 = y1, y2 = y2, y3 = y3, eta1 = eta.eq1.1, eta2 = eta.eq2.1, eta3 = eta.eq3.1, 
                            sig1 = eta.eq1.2, sig2 = eta.eq2.2, sig3 = eta.eq3.2, nu1 = eta.eq1.3, nu2 = eta.eq2.3, nu3 = eta.eq3.3, 
                            theta1 = theta1, theta2 = theta2, theta3 = theta3, beta.theta3 = beta.theta3,
                            margins = x$margins, cop.sel = cop.sel, dof.vec = x$dof.vec, obj.c = obj.cB, simplified = x$simplified, eq.theta = eq.theta,
                            newdata = newdata, min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, eq = NULL, 
                            left.trunc.v = x$left.trunc.v, cond.var = x$cond.var, nC = nC, sign.theta = sI)  }else{    
    
  res.I <- quadinf(cdf.pcc, xa = lB, xb = uB, y1 = y1, y2 = y2, y3 = y3, eta1 = eta.eq1.1, eta2 = eta.eq2.1, eta3 = eta.eq3.1, 
                            sig1 = eta.eq1.2, sig2 = eta.eq2.2, sig3 = eta.eq3.2, nu1 = eta.eq1.3, nu2 = eta.eq2.3, nu3 = eta.eq3.3, 
                            theta1 = theta1, theta2 = theta2, theta3 = theta3, beta.theta3 = beta.theta3,
                            margins = x$margins, cop.sel = cop.sel, dof.vec = x$dof.vec, obj.c = obj.cB, simplified = x$simplified, eq.theta = eq.theta,
                            newdata = newdata, min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, eq = NULL, 
                            left.trunc.v = x$left.trunc.v, cond.var = x$cond.var, nC = nC, sign.theta = sI)                                                   
  res.m      <- res.I$Q  
  res.relerr <- res.I$relerr
  res.niter  <- res.I$niter                                                    
                            }
                            
                                                     


   

if(intervals == TRUE){

for(i in 1:n.sim){

        
 if(method == "adaptive"){  
  res.mCITemp <- suppressMessages(try(distrExIntegrate(cdf.pcc, lower = lB, upper = uB, y1 = y1, y2 = y2, y3 = y3, eta1 = eta.eq1.1S[i], eta2 = eta.eq2.1S[i], eta3 = eta.eq3.1S[i], 
                            sig1 = eta.eq1.2S[i], sig2 = eta.eq2.2S[i], sig3 = eta.eq3.2S[i], nu1 = eta.eq1.3S[i], nu2 = eta.eq2.3S[i], nu3 = eta.eq3.3S[i], 
                            theta1 = theta1S[i], theta2 = theta2S[i], theta3 = theta3S[i], beta.theta3 = BS3[i, ],
                            margins = x$margins, cop.sel = cop.sel, dof.vec = x$dof.vec, obj.c = obj.cB, simplified = x$simplified, eq.theta = eq.theta,
                            newdata = newdata, min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, eq = NULL, 
                            left.trunc.v = x$left.trunc.v, cond.var = x$cond.var, nC = nC, sign.theta = sI) , silent = TRUE)) }else{                      
  res.mCITemp <- suppressMessages(try(quadinf(cdf.pcc, xa = lB, xb = uB, y1 = y1, y2 = y2, y3 = y3, eta1 = eta.eq1.1S[i], eta2 = eta.eq2.1S[i], eta3 = eta.eq3.1S[i], 
                            sig1 = eta.eq1.2S[i], sig2 = eta.eq2.2S[i], sig3 = eta.eq3.2S[i], nu1 = eta.eq1.3S[i], nu2 = eta.eq2.3S[i], nu3 = eta.eq3.3S[i], 
                            theta1 = theta1S[i], theta2 = theta2S[i], theta3 = theta3S[i], beta.theta3 = BS3[i, ],
                            margins = x$margins, cop.sel = cop.sel, dof.vec = x$dof.vec, obj.c = obj.cB, simplified = x$simplified, eq.theta = eq.theta,
                            newdata = newdata, min.dn = x$m1$VC$min.dn, min.pr = x$m1$VC$min.pr, max.pr = x$m1$VC$max.pr, eq = NULL, 
                            left.trunc.v = x$left.trunc.v, cond.var = x$cond.var, nC = nC, sign.theta = sI)$Q , silent = TRUE))      
                                                                                                         } 
   if( is.numeric(res.mCITemp) == FALSE ) next 
   res.mCI[i] <- as.numeric(res.mCITemp)
   
   
}### end intervals

}




list(p123 = res.m, p123s = matrix(res.mCI, 1, n.sim)) #, p1 = p1, p2 = p2, p3 = p3, 
     #theta12 = theta12, theta13 = theta13, theta23 = theta23, 
     #theta12s = epd12s, theta13s = epd13s, theta23s = epd23s)




}