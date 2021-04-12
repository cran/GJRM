ass.ms <- function(BivD, nCa, theta){


if(BivD %in% c("J0","J180","J90","J270")) theta <- ifelse(abs(theta) > 50, 50, abs(theta))
if(BivD %in% c("J90","J270"))             theta <- -theta 

if(BivD %in% c("F")){ signs <- sign(theta) 
                      theta <- ifelse(abs(theta) > 100, 100, abs(theta))
                      theta <- theta*signs }

# if(BivD %in% c("PL")) theta <- ifelse(abs(theta) > 100, 100, abs(theta))
# maybe there is no truncation to be done unless
# I encouter an issue

if(BivD %in% c("C0","C180","C90","C270","G0","G180","G90","G270",
               "GAL0","GAL180","GAL90","GAL270"))                  theta <- ifelse(abs(theta) > 100, 100, abs(theta))
if(BivD %in% c("C90","C270","G90","G270","GAL90","GAL270"))        theta <- -theta


if(BivD %in% c("HO")){ theta <- ifelse(theta == 0, sqrt(.Machine$double.eps), theta)
                       theta <- ifelse(theta == 1, 0.999999,                  theta)
                      }

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
  if(length(theta)==1)   tau <- as.numeric(tau(plackettCopula(theta)))   
  if(length(theta) > 1){ tau <- NA; for(i in 1:length(theta)) tau[i] <- as.numeric(tau(plackettCopula(theta[i])))   }
  if(dim(as.matrix(theta))[2] > 1)  tau <- matrix(tau, nrow = dim(theta)[1], ncol = dim(theta)[2] )  
}


GALs <- c("GAL0","GAL180","GAL90","GAL270")


if(BivD %in% GALs){

  thetaT <- abs(theta)

  if(length(theta)==1)   tau <- as.numeric(tau(galambosCopula(thetaT)))   
  if(length(theta) > 1){ tau <- NA; for(i in 1:length(theta)) tau[i] <- as.numeric(tau(galambosCopula(thetaT[i])))   }
  if(dim(as.matrix(theta))[2] > 1)  tau <- matrix(tau, nrow = dim(theta)[1], ncol = dim(theta)[2] ) 
  
  
  if(BivD %in% c("GAL90","GAL270")) tau <- -tau 
  
  
}



tau.a   <- mean(tau) 
theta.a <- mean(theta) 

list(theta.a = theta.a, theta = theta, tau = tau, tau.a = tau.a)

}
