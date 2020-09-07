Reg2Copost <- function(SemiParFit, VC, theta){

thetaT <- theta 

if(length(thetaT) > 1){

if( length(SemiParFit$fit$teta1) != 0) thetaT[SemiParFit$fit$teta.ind1] <-  theta[SemiParFit$fit$teta.ind1] 
if( length(SemiParFit$fit$teta2) != 0) thetaT[SemiParFit$fit$teta.ind2] <- -theta[SemiParFit$fit$teta.ind2] 

                       }else{

if( length(SemiParFit$fit$teta1) != 0) thetaT <-  theta 
if( length(SemiParFit$fit$teta2) != 0) thetaT <- -theta 

                             }



theta <- thetaT

###################################

tau <- NA

nCa1 <- VC$cta[which(VC$cta[,1] == SemiParFit$fit$Cop1),2] 
nCa2 <- VC$cta[which(VC$cta[,1] == SemiParFit$fit$Cop2),2]

if( length(SemiParFit$fit$teta1) != 0){

if(length(theta) > 1)  ass.msR1 <- ass.ms(SemiParFit$fit$Cop1, nCa1, theta[SemiParFit$fit$teta.ind1])
if(length(theta) == 1) ass.msR1 <- ass.ms(SemiParFit$fit$Cop1, nCa1, theta)


theta1   <- ass.msR1$theta
tau1     <- ass.msR1$tau  

if(length(theta) > 1){
                    theta[SemiParFit$fit$teta.ind1] <- theta1 
                    tau[SemiParFit$fit$teta.ind1]   <- tau1
                     }

if(length(theta) == 1){
                    theta <- theta1 
                    tau   <- tau1 
                       }
        
}



if( length(SemiParFit$fit$teta2) != 0){

if(length(theta) > 1)  ass.msR2 <- ass.ms(SemiParFit$fit$Cop2, nCa2, theta[SemiParFit$fit$teta.ind2])
if(length(theta) == 1) ass.msR2 <- ass.ms(SemiParFit$fit$Cop2, nCa2, theta)

theta2   <- ass.msR2$theta
tau2     <- ass.msR2$tau  

if(length(theta) > 1){ 
                     theta[SemiParFit$fit$teta.ind2] <- theta2 
                     tau[SemiParFit$fit$teta.ind2] <- tau2 
                     }
if(length(theta) == 1){ 
                     theta <- theta2 
                     tau <- tau2
                     } 
                   
                   
}



#if(length(tau) > 1){tau1 <- theta; tau1[1:length(tau),] <- tau; tau <- tau1}


list(theta = theta, tau = tau, theta.a = mean(theta), tau.a = mean(tau) )

}

