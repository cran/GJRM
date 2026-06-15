theta2tau <- function(BivD, nCa, theta, tau.res = FALSE, t.theta = NULL){

tau <- tau.a <- theta1 <- theta2 <- NA

BivD2 <- c("C0C90","C0C270","C180C90","C180C270", "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270",
           "J0J90","J0J270","J180J90","J180J270", "G0G90","G0G270","G180G90","G180G270")




if(BivD %in% BivD2){

Cop1Cop2R <- Cop1Cop2(BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2

teta.ind1 <- as.logical( t.theta >= 0 ) # for + coefs
teta.ind2 <- teta.ind1 == FALSE         # for - coefs

if( length(teta.ind1) == 1 & length(teta.ind2) == 1 & length(theta) > 1 ){

                                                 teta.ind1 <- rep(teta.ind1, length(theta))      
                                                 teta.ind2 <- rep(teta.ind2, length(theta)) 
                                                 theta     <- matrix(theta, nrow = length(theta), ncol = 1)                                                                                                 
                                                 
                                                                         }


ms.dim <- sum(dim(as.matrix(theta)))


if(ms.dim == 2) {

theta1 <- theta[teta.ind1] 
theta2 <- theta[teta.ind2] 

                       }else{

theta1 <- theta[teta.ind1, ] 
theta2 <- theta[teta.ind2, ] 

                             } 


if(Cop1 %in% c("J0", "J180")){ theta1 <- ifelse(abs(theta1) > 30, 30, abs(theta1))
                               theta2 <- ifelse(abs(theta2) > 30, 30, abs(theta2))
                               theta2 <- -theta2
                             } 

if(Cop1 %in% c("C0", "C180")){ theta1 <- ifelse(abs(theta1) > 28, 28, abs(theta1))
                               theta2 <- ifelse(abs(theta2) > 28, 28, abs(theta2))
                               theta2 <- -theta2
                             }  
                     
if(Cop1 %in% c("G0", "G180")){ theta1 <- ifelse(abs(theta1) > 17, 17, abs(theta1))
                               theta2 <- ifelse(abs(theta2) > 17, 17, abs(theta2))
                               theta2 <- -theta2
                             }  
                     
if(Cop1 %in% c("GAL0", "GAL180")){ theta1 <- ifelse(abs(theta1) > 45, 45, abs(theta1))
                                   theta2 <- ifelse(abs(theta2) > 45, 45, abs(theta2))
                                   theta2 <- -theta2
                                 }                       
                     
if(ms.dim == 2) {

theta[teta.ind1] <- theta1
theta[teta.ind2] <- theta2

                        }else{

theta[teta.ind1, ] <- theta1
theta[teta.ind2, ] <- theta2 

                             }



if(tau.res == TRUE){

if(Cop1 %in% c("J0", "J180")) tau <- BiCopPar2Tau(family = 6, par = abs(theta))   ## this function is slow ##
if(Cop1 %in% c("G0", "G180")) tau <- 1-1/abs(theta)  
if(Cop1 %in% c("C0", "C180")) tau <- abs(theta)/(abs(theta)+2)  

if(ms.dim == 2)  tau <- as.numeric(tau)   

if(Cop1 %in% c("GAL0","GAL180")){

  if(ms.dim == 2)  tau <- as.numeric(copula::tau(galambosCopula(abs(theta))))   
  if(ms.dim > 2 ){ tau <- NA; for(i in 1:length(theta)) tau[i] <- as.numeric(copula::tau(galambosCopula(abs(theta)[i])))   }
  if(ms.dim > 2 )  tau <- matrix(tau, nrow = dim(theta)[1], ncol = dim(theta)[2] ) 
  

}

if(ms.dim == 2) {

tau[teta.ind1] <-  tau
tau[teta.ind2] <- -tau

                        }else{

tau[teta.ind1, ] <-  tau[teta.ind1, ]
tau[teta.ind2, ] <- -tau[teta.ind2, ] 

                             }


if(is.matrix(tau)) dimnames(tau)[[2]] <- rep("tau", dim(tau)[2])

tau.a <- mean(tau) 

} # tau condition





}else{



if(BivD %in% c("J0","J180","J90","J270")) theta <- ifelse(abs(theta) > 30, 30, abs(theta)) # changed from 50 to 30 - 05/08/2021, due to BiCopPar2Tau
if(BivD %in% c("J90","J270"))             theta <- -theta 

if(BivD %in% c("F")){ signs <- sign(theta) 
                      theta <- ifelse(abs(theta) > 35, 35, abs(theta)) # changed from 100 - 17/09/2021, due to BiCopPar2Tau
                      theta <- theta*signs }

# if(BivD %in% c("PL")) theta <- ifelse(abs(theta) > 100, 100, abs(theta))
# maybe there is no truncation to be done unless
# I encouter an issue


if(BivD %in% c("C0","C180","C90","C270"))                  theta <- ifelse(abs(theta) > 28, 28, abs(theta)) # based on BiCopPar2Tau 

if(BivD %in% c("G0","G180","G90","G270"))                  theta <- ifelse(abs(theta) > 17, 17, abs(theta)) # based on BiCopPar2Tau 


if(BivD %in% c("GAL0","GAL180","GAL90","GAL270"))                  theta <- ifelse(abs(theta) > 45, 45, abs(theta))

if(BivD %in% c("C90","C270","G90","G270","GAL90","GAL270"))        theta <- -theta


if(BivD %in% c("HO")){ theta <- ifelse(theta == 0, sqrt(.Machine$double.eps), theta)
                       theta <- ifelse(theta == 1, 0.999999,                  theta)
                      }


if(tau.res == TRUE){

tau <- 0
if(BivD %in% c("F", "J0","J180","J90","J270")) tau <- BiCopPar2Tau(family = nCa, par = theta)      ## this function is slow ##
if(BivD == "AMH")                              tau <- 1 - (2/3)/theta^2*(theta + (1-theta)^2*log(1-theta))
if(BivD == "FGM")                              tau <- 2/9*theta
if(BivD == "HO")                               tau <- 1 - theta  
if(BivD %in% c("T","N"))                       tau <- 2/pi*asin(theta) 
if(BivD %in% c("G0","G180"))                   tau <- 1-1/theta  
if(BivD %in% c("G90","G270"))                  tau <- -(1-1/abs(theta))  
if(BivD %in% c("C0","C180"))                   tau <- theta/(theta+2)  
if(BivD %in% c("C90","C270"))                  tau <- -(abs(theta)/(abs(theta)+2))  

if(length(theta)==1) tau <- as.numeric(tau)   



if(BivD == "PL"){
  if(length(theta)==1)   tau <- as.numeric(copula::tau(plackettCopula(theta)))   ## this function is slow ##
  if(length(theta) > 1){ tau <- NA; for(i in 1:length(theta)) tau[i] <- as.numeric(copula::tau(plackettCopula(theta[i])))   }
  if(dim(as.matrix(theta))[2] > 1)  tau <- matrix(tau, nrow = dim(theta)[1], ncol = dim(theta)[2] )  
}


GALs <- c("GAL0","GAL180","GAL90","GAL270")


if(BivD %in% GALs){

  thetaT <- abs(theta)

  if(length(theta)==1)   tau <- as.numeric(copula::tau(galambosCopula(thetaT)))   
  if(length(theta) > 1){ tau <- NA; for(i in 1:length(theta)) tau[i] <- as.numeric(copula::tau(galambosCopula(thetaT[i])))   }
  if(dim(as.matrix(theta))[2] > 1)  tau <- matrix(tau, nrow = dim(theta)[1], ncol = dim(theta)[2] ) 
  
  
  if(BivD %in% c("GAL90","GAL270")) tau <- -tau 
  
  
}


if(is.matrix(tau)) dimnames(tau)[[2]] <- rep("tau", dim(tau)[2])

tau.a   <- mean(tau) 

} # tau condition

} # else


theta.a <- mean(theta) 

list(theta.a = theta.a, theta = theta, tau = tau, tau.a = tau.a)

}
