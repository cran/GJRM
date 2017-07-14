ass.ms <- function(BivD, nCa, theta){

if(BivD %in% c("J0","J180","J90","J270"))  theta <- ifelse(abs(theta) > 50, 50, abs(theta))
if(BivD %in% c("J90","J270"))              theta <- -theta 

if(BivD %in% c("F")){ signs <- sign(theta) 
                      theta <- ifelse(abs(theta) > 100, 100, abs(theta))
                      theta <- theta*signs }

# f(BivD %in% c("PL")) theta <- ifelse(abs(theta) > 100, 100, abs(theta))
# maybe there is no truncation to be done unless
# I encouter an issue

if(BivD %in% c("C0","C180","G0","G180","C90","C270","G90","G270")) theta <- ifelse(abs(theta) > 100, 100, abs(theta))
if(BivD %in% c("C90","C270","G90","G270"))                         theta <- -theta


if(BivD %in% c("HO")){ theta <- ifelse(theta == 0, 0.0000001 , theta)
                       theta <- ifelse(theta == 1, 0.9999999 , theta)
                      }


if(!(BivD %in% c("AMH","FGM","PL","HO"))) tau <- BiCopPar2Tau(family = nCa, par = theta)
if(BivD == "AMH")                         tau <- t(as.numeric(1 - (2/3)/theta^2*(theta + (1-theta)^2*log(1-theta))))
if(BivD == "FGM")                         tau <- t(as.numeric(2/9*theta))
if(BivD == "HO")                          tau <- 1 - theta  

if(BivD == "PL"){
  if(length(theta)==1)   tau <- as.numeric(tau(plackettCopula(theta)))   
  if(length(theta) > 1){ tau <- NA; for(i in 1:length(theta)) tau[i] <- as.numeric(tau(plackettCopula(theta[i])))   }
  if(dim(as.matrix(theta))[2] > 1)  tau <- matrix(tau, nrow = dim(theta)[1] , ncol = dim(theta)[2] )
  #tau <- t(tau)
  
  #if(dim(as.matrix(theta))[2] == 1) tau <- t(t(tau)) 
  #if(length(theta) == 1) tau <- as.numeric(tau) 
  
  
}

tau.a   <- mean(tau) 
theta.a <- mean(theta) 

list(theta.a = theta.a, theta = theta, tau = tau, tau.a = tau.a)

}
