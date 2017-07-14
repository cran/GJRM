numch <- function(funcD, para1, para2){

#fd.prec  <- 10^(-7)
#fd.prec2 <- 10^-3
#
#
#para11 <- para1 + fd.prec
#para22 <- para2 + fd.prec2
#
#f00 <- funcD(para1, para2)
#
#
#
#F2.fd.dsigma.plus <- pPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus)
#      dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
#      dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
#      d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
#      d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
#

eps <- 1e-06 
eps2 <- eps*eps 

f00 <- funcD(para1, para2)  

t01 <- t10 <- t11 <- cbind(para1,para2)

t01[,1] <- t01[,1] + eps 
t10[,2] <- t10[,2] + eps 
t11[,1] <- t11[,1] + eps 
t11[,2] <- t11[,2] + eps 

f01 <- funcD(t01[,1],t01[,2]) 
f10 <- funcD(t10[,1],t10[,2]) 
f11 <- funcD(t11[,1],t11[,2]) 

(f11 - f01 - f10 + f00)/eps2  
 

}
    
  