copgHs2 <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD, par2 = NULL,
                    min.pr, max.pr, min.dn = 1e-40){



########################################################################################
# Rotations
########################################################################################

if(BivD %in% c("C90","J90","G90","GAL90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(BivD %in% c("C180","J180","G180","GAL180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(BivD %in% c("C270","J270","G270","GAL270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################



if(BivD %in% c("GAL0","GAL90","GAL180","GAL270")){

c.copula.be2 <- p1*(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))


c.copula2.be1be2 <- (1+teta*(-log(p2))^(1+teta)*(1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2)-((1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))/((-log(p1))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))+1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))))*
exp((1/(-log(p1))^teta+1/(-log(p2))^teta)^-(1/teta))




c.copula2.be1 <- p2*((teta*(1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)-(-log(p1))^teta*(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))/((-log(p1))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2-
(1-1/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))/p1

}







if(BivD=="HO"){

c.copula.be2 <- (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2 
 
 
c.copula2.be1be2 <- ((-log(p1))^(1/teta - 1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(2 * (teta - 1)) - (-log(p1))^(1/teta - 
    1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * (teta - 1)/teta) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/(p1 * 
    p2) 
    
 
c.copula2.be1 <- ((-log(p1))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(2 * 
    (teta - 1)) - ((-log(p1))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta + ((-log(p1))^(1/teta - 
    1) + (-log(p1))^(1/teta - 2) * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1))) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)/p1^2  
 
}


if(BivD=="PL"){

c.copula.be2 <- (teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * 
    teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * (teta - 1)) 

c.copula2.be1be2 <- -(((2 * (teta - 1) - 4 * teta)/sqrt(((p1 + p2) * (teta - 1) + 
    1)^2 - 4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * ((p1 + 
    p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p2 * teta)) * (teta - 1)/(((p1 + p2) * 
    (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5))/4) 
  
c.copula2.be1 <- -((2/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1))) - 0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 
    4 * (p2 * teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * 
    (p1 * p2 * teta * (teta - 1)))^1.5)) * (teta - 1)/4)
    

}





if(BivD=="FGM"){

c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))


c.copula2.be1be2 <- 1 + teta * (1 - 2 * p1) * (1 - 2 * p2)

c.copula2.be1 <- -(2 * (p2 * teta * (1 - p2)))



}


if(BivD=="AMH"){

c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))
 

c.copula2.be1be2 <- (1 - teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - p1) * (1 - 
    p2))


c.copula2.be1 <- -(p2 * teta * (1 - p2) * (2 - 2 * (p1 * teta * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2)


}





if(BivD=="N"){

c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   ) 


c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 

c.copula2.be1 <- dnorm((qnorm(p2)-teta*qnorm(p1))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p1)^2/2)     


}


if(BivD == "T"){
                    
c.copula.be2 <- BiCopHfunc2(p1, p2, family = 2, par = teta, par2 = par2) 

c.copula2.be1be2 <- BiCopPDF(p1, p2, family = 2, par = teta, par2 = par2) 

c.copula2.be1 <-  BiCopHfuncDeriv(p2, p1, 2, par=teta, par2 = par2, deriv="u2")     


}


if(BivD=="F"){
 
  c.copula.be2 <-   (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + 
                     exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))


c.copula2.be1be2 <- (exp((1 + p1 + p2)* teta)* (-1 + exp(teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2


c.copula2.be1 <-   (exp(teta + 
  p1* teta)* (-1 + exp(p2* teta))* (-exp(teta) + exp(
   p2* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2
  

}



if(BivD %in% c("C0","C90","C180","C270")){

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  
c.copula2.be1be2 <- p1^(-1 - teta)* p2^(-1 - teta)* (-1 + p1^-teta + p2^-teta)^(-2 - 1/
  teta) *(1 + teta)
     
c.copula2.be1 <- (p1^(-2 + teta)* p2^teta* (-1 + p1^-teta + p2^-teta)^(-1/
  teta)* (-1 + p2^teta)* (1 + teta))/(p2^teta - 
  p1^teta* (-1 + p2^teta))^2
 

  
}



if(BivD %in% c("G0","G90","G180","G270")){

  c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2

c.copula2.be1be2 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta) *(-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
   teta))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-log(p2))^(-1 + teta))/(p1 *p2)

c.copula2.be1 <-(1/(p1^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* (-log(p1))^(-2 + 
  teta) *((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* ((-log(p1))^
    teta* (log(p1) + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta)) + (1 - teta + log(p1))* (-log(p2))^teta)


}


if(BivD %in% c("J0","J90","J180","J270")){

  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))

c.copula2.be1be2 <- (1 - p1)^(-1 + 
  teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^(-1 + 
  teta)* (-(-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) + teta)


c.copula2.be1 <- (1 - p1)^(-2 + 
  teta)* (-1 + (1 - p2)^teta)* ((1 - p1)^
   teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^teta* (-1 + teta)

  

}



if(BivD %in% c("C90","J90","G90","GAL90") ) {


c.copula.be2  <- 1 - c.copula.be2

}  

if(BivD %in% c("C180","J180","G180","GAL180") ) {


c.copula.be2  <- 1 - c.copula.be2


}  


#if(BivD %in% c("C270","J270","G270") ) {
#
#
#c.copula.be2 <- c.copula.be2
#
#
#}   





  
c.copula.be2 <- mm(c.copula.be2, min.pr = min.pr, max.pr = max.pr) 

c.copula2.be1be2 <- ifelse(c.copula2.be1be2 < min.dn, min.dn, c.copula2.be1be2)

list( c.copula2.be1 = ifef(c.copula2.be1), c.copula.be2 = ifef(c.copula.be2), c.copula2.be1be2 = ifef(c.copula2.be1be2) )  





}
