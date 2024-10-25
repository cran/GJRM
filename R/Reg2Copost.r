Reg2Copost <- function(SemiParFit, VC, theta, tau.res = FALSE){

thetaT <- theta; tau <- NA


if(length(thetaT) > 1){

if( length(SemiParFit$fit$teta1) != 0) thetaT[SemiParFit$fit$teta.ind1] <-  theta[SemiParFit$fit$teta.ind1] 
if( length(SemiParFit$fit$teta2) != 0) thetaT[SemiParFit$fit$teta.ind2] <- -theta[SemiParFit$fit$teta.ind2] 

                       }else{

if( length(SemiParFit$fit$teta1) != 0) thetaT <-  theta 
if( length(SemiParFit$fit$teta2) != 0) thetaT <- -theta 

                             }



theta <- thetaT

###################################


nCa1 <- VC$cta[which(VC$cta[,1] == SemiParFit$fit$Cop1),2] 
nCa2 <- VC$cta[which(VC$cta[,1] == SemiParFit$fit$Cop2),2]

if( length(SemiParFit$fit$teta1) != 0){

if(length(theta) > 1)  ass.msR1 <- theta2tau(SemiParFit$fit$Cop1, nCa1, theta[SemiParFit$fit$teta.ind1], tau.res = tau.res)
if(length(theta) == 1) ass.msR1 <- theta2tau(SemiParFit$fit$Cop1, nCa1, theta, tau.res = tau.res)


theta1   <- ass.msR1$theta
if(tau.res == TRUE) tau1 <- ass.msR1$tau  

if(length(theta) > 1){
                    theta[SemiParFit$fit$teta.ind1] <- theta1 
                    if(tau.res == TRUE) tau[SemiParFit$fit$teta.ind1] <- tau1
                     }

if(length(theta) == 1){
                    theta <- theta1 
                    if(tau.res == TRUE) tau <- tau1 
                       }
        
}



if( length(SemiParFit$fit$teta2) != 0){

if(length(theta) > 1)  ass.msR2 <- theta2tau(SemiParFit$fit$Cop2, nCa2, theta[SemiParFit$fit$teta.ind2], tau.res = tau.res)
if(length(theta) == 1) ass.msR2 <- theta2tau(SemiParFit$fit$Cop2, nCa2, theta, tau.res = tau.res)

theta2   <- ass.msR2$theta
if(tau.res == TRUE) tau2  <- ass.msR2$tau  

if(length(theta) > 1){ 
                     theta[SemiParFit$fit$teta.ind2] <- theta2 
                     if(tau.res == TRUE) tau[SemiParFit$fit$teta.ind2] <- tau2 
                     }
if(length(theta) == 1){ 
                     theta <- theta2 
                     if(tau.res == TRUE) tau <- tau2
                     } 
                   
                   
}



list(theta = theta, theta.a = mean(theta), tau = tau, tau.a = mean(tau))

}

