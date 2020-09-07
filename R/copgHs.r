copgHs <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD, nu = NULL, nu.st = NULL, CLM = FALSE,
                   min.dn, min.pr, max.pr){

########################################################################################



c.copula.dof.st <- c.copula2.be1dof.st <- c.copula2.be2dof.st <- c.copula2.dof2.st <- c.copula2.thdof.st <- 1


cjg <- c("C0","J0","G0","C90","J90","G90","C180","J180","G180","C270","J270","G270","PL")

if(BivD %in% cjg) {

derteta.derteta.st <- der2teta.derteta.stteta.st <- exp(teta.st) 
   
}   


if(BivD %in% c("N", "T","FGM","AMH") ) {

derteta.derteta.st <- 1/cosh(teta.st)^2
der2teta.derteta.stteta.st <- -(2 * (sinh(teta.st) * cosh(teta.st))/(cosh(teta.st)^2)^2)

       
} 


if(BivD %in% c("F") ) {

derteta.derteta.st <- 1
der2teta.derteta.stteta.st <- 0
       
} 


if(BivD %in% c("HO") ) {

derteta.derteta.st         <- exp(-teta.st)/(1 + exp(-teta.st))^2 
der2teta.derteta.stteta.st <- -((1 - 2 * (exp(-teta.st)/(1 + exp(-teta.st)))) * exp(-teta.st)/(1 + exp(-teta.st))^2) 
       
} 




###################################################################################
     
#if(BivD %in% c("T") ) derdof.derdof.st <- der2dof.derdof.stdof.st <- exp(nu.st)

########################################################################################
# Rotations
########################################################################################

if(BivD %in% c("C90","J90","G90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(BivD %in% c("C180","J180","G180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(BivD %in% c("C270","J270","G270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################


if(BivD=="HO"){


c.copula.be1 <- (-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p1 

c.copula.be2 <- (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2 
  
  
  c.copula.thet <- -((((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1)/teta) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)) 
  


c.copula2.be1 <- ((-log(p1))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(2 * 
    (teta - 1)) - ((-log(p1))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta + ((-log(p1))^(1/teta - 
    1) + (-log(p1))^(1/teta - 2) * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1))) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)/p1^2 
 

 c.copula2.be2 <- ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(2 * 
    (teta - 1)) - ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta + ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((-log(p2))^(1/teta - 1) + 
    (-log(p2))^(1/teta - 2) * (1/teta - 1)))) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)/p2^2 


c.copula2.be1be2 <- ((-log(p1))^(1/teta - 1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(2 * (teta - 1)) - (-log(p1))^(1/teta - 
    1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * (teta - 1)/teta) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/(p1 * 
    p2) 
   
  
c.copula2.be1t <- ((-log(p1))^(1/teta - 1) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) + ((-log(p1))^(1/teta - 
    1) - ((-log(p1))^(1/teta - 1) + (-log(p1))^(1/teta - 1) * 
    log(-log(p1))/teta)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1)/teta - (-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
    ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
        log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1)/teta)) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p1
  
c.copula2.be2t <- ((-log(p2))^(1/teta - 1) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) + ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((-log(p2))^(1/teta - 1) - 
    ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta))/teta - 
    (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * 
        log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
        log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * 
        ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 1)/teta)) * 
    exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2 
  
  
bit1.th2ATE <- -(((((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
    log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1)/teta + (1 - ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta) * 
    log((-log(p1))^(1/teta) + (-log(p2))^(1/teta))) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
    ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
        log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1)/teta) - (((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
    log(-log(p2))) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) + (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta)/teta) + 
    ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 1) * 
        (((-log(p1))^(1/teta) - ((-log(p1))^(1/teta) * log(-log(p1))/teta + 
            2 * (-log(p1))^(1/teta))) * log(-log(p1)) + ((-log(p2))^(1/teta) - 
            ((-log(p2))^(1/teta) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta))) * 
            log(-log(p2)))/teta)/teta) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta))


}











if(BivD=="PL"){


c.copula.be1 <- (teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p2 * 
    teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * (teta - 1)) 

c.copula.be2 <- (teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * 
    teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * (teta - 1)) 
  
  
  c.copula.thet <- (p1 + p2 - 0.5 * ((2 * (((p1 + p2) * (teta - 1) + 1) * (p1 + 
    p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta))/sqrt(((p1 + 
    p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))))/(2 * 
    (teta - 1)) - 2 * (((p1 + p2) * (teta - 1) + 1 - sqrt(((p1 + 
    p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1))))/(2 * 
    (teta - 1))^2) 
  


c.copula2.be1 <- -((2/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1))) - 0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 
    4 * (p2 * teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * 
    (p1 * p2 * teta * (teta - 1)))^1.5)) * (teta - 1)/4)
 

 c.copula2.be2 <- -((2/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1))) - 0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 
    4 * (p1 * teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * 
    (p1 * p2 * teta * (teta - 1)))^1.5)) * (teta - 1)/4) 


c.copula2.be1be2 <- -(((2 * (teta - 1) - 4 * teta)/sqrt(((p1 + p2) * (teta - 1) + 
    1)^2 - 4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * ((p1 + 
    p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p2 * teta)) * (teta - 1)/(((p1 + p2) * 
    (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5))/4) 
   
  
c.copula2.be1t <- (1 - 0.5 * ((2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - p2 * (4 * 
    (teta - 1) + 4 * teta))/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * (((p1 + 
    p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 
    1) + 4 * teta)) * (2 * ((p1 + p2) * (teta - 1) + 1) - 4 * 
    (p2 * teta)) * (teta - 1)/(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))^1.5)))/(2 * (teta - 1)) - 
    2 * ((teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 
        4 * (p2 * teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 
        1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * 
        (teta - 1))^2) 
  
c.copula2.be2t <- (1 - 0.5 * ((2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - p1 * (4 * 
    (teta - 1) + 4 * teta))/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * (((p1 + 
    p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 
    1) + 4 * teta)) * (2 * ((p1 + p2) * (teta - 1) + 1) - 4 * 
    (p1 * teta)) * (teta - 1)/(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))^1.5)))/(2 * (teta - 1)) - 
    2 * ((teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 
        4 * (p1 * teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 
        1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * 
        (teta - 1))^2) 
  
  
bit1.th2ATE <- -(((2 * (p1 + p2)^2 - 8 * (p1 * p2))/sqrt(((p1 + p2) * (teta - 
    1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * 
    (((p1 + p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * 
    (teta - 1) + 4 * teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))^1.5))/(4 * (teta - 1)) + 
    (4 * (p1 + p2 - 0.5 * ((2 * (((p1 + p2) * (teta - 1) + 1) * 
        (p1 + p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta))/sqrt(((p1 + 
        p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 
        1))))) - 16 * (((p1 + p2) * (teta - 1) + 1 - sqrt(((p1 + 
        p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 
        1)))) * (teta - 1)/(2 * (teta - 1))^2))/(2 * (teta - 
        1))^2) 


}





if(BivD=="N"){


c.copula.be1 <- pnorm( (qnorm(p2) - teta*qnorm(p1))/sqrt(1 - teta^2)   )                            
c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   )  

c.copula.thet <- dbinorm(qnorm(p1),qnorm(p2), cov12=teta)

c.copula2.be1 <- dnorm((qnorm(p2)-teta*qnorm(p1))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p1)^2/2)     
c.copula2.be2 <- dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p2)^2/2) 

c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 

c.copula2.be1t <- (-(dnorm((qnorm(p2) - tanh(teta.st) * qnorm(p1))/sqrt(1 - tanh(teta.st)^2)) * 
     (1/cosh(teta.st)^2 * qnorm(p1)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p2) - 
         tanh(teta.st) * qnorm(p1)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
         tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
         tanh(teta.st)^2)^2)))/derteta.derteta.st                                                       
                                                                                                                                                                      
c.copula2.be2t <- (-(dnorm((qnorm(p1) - tanh(teta.st) * qnorm(p2))/sqrt(1 - tanh(teta.st)^2)) *  # this as well
    (1/cosh(teta.st)^2 * qnorm(p2)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p1) - 
        tanh(teta.st) * qnorm(p2)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
        tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
        tanh(teta.st)^2)^2)))/derteta.derteta.st
        
bit1.th2ATE <- (0.5 * (pi * teta/(pi * sqrt(1 - teta^2))^2) - 0.5 * 
    ((teta * (qnorm(p1) * (qnorm(p1) - 2 * (teta * qnorm(p2))) + 
        qnorm(p2)^2)/(1 - teta^2) - qnorm(p1) * qnorm(p2))/(pi * 
        (1 - teta^2)))) * exp(-(0.5 * ((qnorm(p1) * (qnorm(p1) - 
    2 * (teta * qnorm(p2))) + qnorm(p2)^2)/(1 - teta^2))))/sqrt(1 - 
    teta^2)


}


if(BivD=="T"){
 

c.copula.be1 <- BiCopHfunc1(p1, p2, family = 2, par = teta, par2 = nu)                         
c.copula.be2 <- BiCopHfunc2(p1, p2, family = 2, par = teta, par2 = nu) 

c.copula.thet <- ( 1 + (qt(p1,nu)^2+qt(p2,nu)^2-2*teta*qt(p1,nu)*qt(p2,nu))/(nu*(1-teta^2)))^(-nu/2)/(2*pi*sqrt(1-teta^2))

#funcD <- function(nu) BiCopCDF(p1, p2, family = 2, par = teta, par2 = nu)
 
# mayhave to be careful with nu values here as it must be > 2  
 
#nde <- numgh(funcD, nu) 
#
#c.copula.dof   <- nde$fi
#c.copula2.dof2 <- nde$se


c.copula2.be1 <-  BiCopHfuncDeriv(p2, p1, 2, par=teta, par2=nu, deriv="u2")     
c.copula2.be2 <-  BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=nu, deriv="u2") 

c.copula2.be1be2 <- BiCopPDF(p1, p2, 2, par=teta, par2=nu) 


c.copula2.be1t <- BiCopHfuncDeriv(p2, p1, 2, par=teta, par2=nu, deriv="par")                                                    
                                                                                                                                                                      
c.copula2.be2t <- BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=nu, deriv="par")


#c.copula2.be1dof <- BiCopHfuncDeriv(p2, p1, 2, par=teta, par2=nu, deriv="par2")                                                    
#                                                                                                                                                                      
#c.copula2.be2dof <- BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=nu, deriv="par2")
                                                                                                                                                                                                                 
bit1.th2ATE <- -((1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 2 * teta * qt(p1,nu) * qt(p2,nu))/(nu * (1 - teta^2)))^((-nu/2) - 
    1) * ((-nu/2) * (2 * qt(p1,nu) * qt(p2,nu)/(nu * (1 - teta^2)) - (qt(p1,nu)^2 + 
    qt(p2,nu)^2 - 2 * teta * qt(p1,nu) * qt(p2,nu)) * (nu * (2 * teta))/(nu * (1 - 
    teta^2))^2))/(2 * pi * sqrt(1 - teta^2)) - (1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 
    2 * teta * qt(p1,nu) * qt(p2,nu))/(nu * (1 - teta^2)))^(-nu/2) * (2 * pi * 
    (0.5 * (2 * teta * (1 - teta^2)^-0.5)))/(2 * pi * sqrt(1 - 
    teta^2))^2)       

#funcD1 <- function(teta, nu) BiCopCDF(p1, p2, family = 2, par = teta, par2 = nu)
#
#c.copula2.tetadof <- numch(funcD1, teta, nu)


}








if(BivD=="F"){

# 1 - exp(-teta) = -expm1(-teta)
# recall log1p as well if needed. 
# epcl <- -expm1(-teta) 

  c.copula.be1 <- (exp(teta)* (-1 + exp(p2* teta)))/(-exp((p1 + p2)* teta) + exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))

  c.copula.be2 <- (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))

  c.copula.thet <-  (exp(teta)* (1/(-1 + exp(teta)) + (-1 - exp(p2* teta)* (-1 + p1) + p1 - 
     exp(p1* teta)* (-1 + p2) + p2)/(
    exp((p1 + p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))* teta + 
 log((exp(-(p1 + p2)* teta)* (-exp((p1 + p2)* teta) + 
     exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))/(-1 + exp(teta))))/teta^2 


c.copula2.be1 <-   (exp(teta + 
  p1* teta)* (-1 + exp(p2* teta))* (-exp(teta) + exp(
   p2* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2
  
 
 c.copula2.be2 <-     (exp(teta + 
  p2* teta)* (-1 + exp(p1* teta))* (-exp(teta) + exp(
   p1* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2


c.copula2.be1be2 <- (exp((1 + p1 + p2)* teta)* (-1 + exp(teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2


c.copula2.be1t <- (exp(teta + 
  p1 *teta)* (exp(2* p2* teta)* (-1 + p1) + exp(teta)* p1 - 
   exp(p2* teta)* (-1 + p1 + exp(teta)* (p1 - p2) + p2)))/(exp((p1 + 
     p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2

c.copula2.be2t <-  (exp(teta + 
  p2 *teta)* (exp(2* p1* teta)* (-1 + p2) + exp(teta)* p2 - 
   exp(p1* teta)* (-1 + p2 + exp(teta)* (p2 - p1) + p1)))/(exp((p2 + 
     p1)* teta) - exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2


bit1.th2 <- bit1.th2ATE <- 1/teta^2 * ((1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))) - 2 * teta/(teta^2)^2 * 
    log(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * 
        p1)) * (1 - exp(-teta * p2)))) + (1/teta^2 * ((1/(1 - 
    exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - exp(-teta)) * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))) - -1/teta * ((1/(1 - exp(-teta)) * (exp(-teta) + 
    (exp(-teta * p1) * p1 * (exp(-teta * p2) * p2) - exp(-teta * 
        p1) * p1 * p1 * (1 - exp(-teta * p2)) + (exp(-teta * 
        p1) * p1 * (exp(-teta * p2) * p2) - (1 - exp(-teta * 
        p1)) * (exp(-teta * p2) * p2 * p2)))) + exp(-teta)/(1 - 
    exp(-teta))^2 * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) + (exp(-teta)/(1 - exp(-teta))^2 * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - (exp(-teta)/(1 - exp(-teta))^2 + 
    exp(-teta) * (2 * (exp(-teta) * (1 - exp(-teta))))/((1 - 
        exp(-teta))^2)^2) * ((1 - exp(-teta)) - (1 - exp(-teta * 
    p1)) * (1 - exp(-teta * p2)))))/(1/(1 - exp(-teta)) * ((1 - 
    exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) + 
    (1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * 
        (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) * (1/(1 - 
        exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
        exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2))))^2))
    
    
    
    
    
    
    


}





if(BivD %in% c("C0","C90","C180","C270")){



c.copula.be1 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta)))

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  
  
  c.copula.thet <- ((-1 + p1^-teta + p2^-teta)^(-1/
  teta) *((teta *(p2^teta *log(p1) + p1^teta* log(p2)))/(
   p2^teta - p1^teta* (-1 + p2^teta)) + 
   log(-1 + p1^-teta + p2^-teta)))/teta^2
  


c.copula2.be1 <- (p1^(-2 + teta)* p2^teta* (-1 + p1^-teta + p2^-teta)^(-1/
  teta)* (-1 + p2^teta)* (1 + teta))/(p2^teta - 
  p1^teta* (-1 + p2^teta))^2
 

 c.copula2.be2 <- (p2^(-2 + teta)* p1^teta* (-1 + p2^-teta + p1^-teta)^(-1/
  teta)* (-1 + p1^teta)* (1 + teta))/(p1^teta - 
  p2^teta* (-1 + p1^teta))^2


c.copula2.be1be2 <- p1^(-1 - teta)* p2^(-1 - teta)* (-1 + p1^-teta + p2^-teta)^(-2 - 1/
  teta) *(1 + teta)
   
  
c.copula2.be1t <- ((1 + 1/teta) * (log(p1)/p1^teta + log(p2)/p2^teta)/(1/p1^teta + 
    1/p2^teta - 1)^(1/teta + 2) + log(1/p1^teta + 1/p2^teta - 
    1)/(teta^2 * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)))/p1^(1 + 
    teta) + (1/p1^(1 + teta) - (1/p1^(1 + teta) + teta * log(p1)/p1^(1 + 
    teta)))/(teta * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)) 
  
 c.copula2.be2t <- ((1 + 1/teta) * (log(p1)/p1^teta + log(p2)/p2^teta)/(1/p1^teta + 
    1/p2^teta - 1)^(1/teta + 2) + log(1/p1^teta + 1/p2^teta - 
    1)/(teta^2 * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)))/p2^(1 + 
    teta) + (1/p2^(1 + teta) - (1/p2^(1 + teta) + teta * log(p2)/p2^(1 + 
    teta)))/(teta * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta))
  
bit1.th2ATE <- (teta * (2 * (p1^teta * log(p2) * (p2^teta * (p1^teta - 
    1) - p1^teta) * (teta - log(1/p1^teta + 1/p2^teta - 1))) + 
    2 * (p2^teta * ((p2^teta * (p1^teta - 1) - p1^teta) * (teta - 
        log(1/p1^teta + 1/p2^teta - 1)) + p1^teta * teta * (1 + 
        teta) * log(p2)) * log(p1)) + teta * (p1^teta * log(p2)^2 * 
    (p1^teta + p2^teta * teta * (p1^teta - 1)) + p2^teta * log(p1)^2 * 
    (p1^teta * teta * (p2^teta - 1) + p2^teta))) - (2 * teta - 
    log(1/p1^teta + 1/p2^teta - 1)) * log(1/p1^teta + 1/p2^teta - 
    1) * (p1^teta + p2^teta - p1^teta * p2^teta)^2)/(teta^4 * 
    (1/p1^teta + 1/p2^teta - 1)^(1/teta) * (p2^teta - p1^teta * 
    (p2^teta - 1))^2)
   
    
}








if(BivD %in% c("G0","G90","G180","G270")){


  c.copula.be1 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta)* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/teta))/p1


  c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2


  c.copula.thet <-   (1/(teta^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/
  teta)* (-teta* (-log(p1))^
    teta* log(-log(p1)) + ((-log(p1))^teta + (-log(p2))^
      teta)* log((-log(p1))^teta + (-log(p2))^teta) - 
   teta *(-log(p2))^teta* log(-log(p2)))

  
c.copula2.be1 <-(1/(p1^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* (-log(p1))^(-2 + 
  teta) *((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* ((-log(p1))^
    teta* (log(p1) + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta)) + (1 - teta + log(p1))* (-log(p2))^teta)

                
c.copula2.be2 <- (1/(p2^2))*exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/
  teta)))* (-log(p2))^(-2 + 
  teta) *((-log(p2))^teta + (-log(p1))^teta)^(-2 + 1/
  teta)* ((-log(p2))^
    teta* (log(p2) + ((-log(p2))^teta + (-log(p1))^teta)^(1/
      teta)) + (1 - teta + log(p2))* (-log(p1))^teta)


c.copula2.be1be2 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta) *(-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
   teta))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-log(p2))^(-1 + teta))/(p1 *p2)


 c.copula2.be1t <- ((-log(p1))^(-1 + teta) * (((-log(p1))^teta * log(-log(p1)) + 
    (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    2) * (1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    1) * log((-log(p1))^teta + (-log(p2))^teta)/teta^2) + ((-log(p1))^(-1 + 
    teta) * log(-log(p1)) - (-log(p1))^(-1 + teta) * (((-log(p1))^teta * 
    log(-log(p1)) + (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    log((-log(p1))^teta + (-log(p2))^teta)/teta)/teta) * ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1)) * exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta))/p1
  
 c.copula2.be2t <- ((-log(p2))^(-1 + teta) * (((-log(p1))^teta * log(-log(p1)) + 
    (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    2) * (1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    1) * log((-log(p1))^teta + (-log(p2))^teta)/teta^2) + ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1) * ((-log(p2))^(-1 + teta) * 
    log(-log(p2)) - (-log(p2))^(-1 + teta) * (((-log(p1))^teta * 
    log(-log(p1)) + (-log(p2))^teta * log(-log(p2))) * ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta - 1) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    log((-log(p1))^teta + (-log(p2))^teta)/teta)/teta)) * exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta))/p2


bit1.th2ATE <- ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 2) * 
    (((-log(p1))^teta + (-log(p2))^teta) * (((-log(p1))^teta + 
        (-log(p2))^teta) * (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 
        1) * log((-log(p1))^teta + (-log(p2))^teta) - 2 * (teta * 
        ((-log(p1))^teta + (-log(p2))^teta * ((((-log(p1))^teta + 
            (-log(p2))^teta)^(1/teta) - 1) * log(-log(p2)) + 
            1)))) * log((-log(p1))^teta + (-log(p2))^teta) + 
        teta * (2 * ((-log(p1))^teta * log(-log(p1)) * (teta * 
            ((-log(p1))^teta + (-log(p2))^teta * ((((-log(p1))^teta + 
                (-log(p2))^teta)^(1/teta) + teta - 1) * log(-log(p2)) + 
                1)) - ((-log(p1))^teta + (-log(p2))^teta) * (((-log(p1))^teta + 
            (-log(p2))^teta)^(1/teta) - 1) * log((-log(p1))^teta + 
            (-log(p2))^teta))) + teta * ((-log(p1))^teta * ((-log(p1))^teta * 
            (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 1) - 
            teta * (-log(p2))^teta) * log(-log(p1))^2 + (-log(p2))^teta * 
            ((-log(p1))^teta * (2 - teta * log(-log(p2))) + (-log(p2))^teta * 
                ((((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 
                  1) * log(-log(p2)) + 2)) * log(-log(p2))))) * 
    exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta))/teta^4


}






if(BivD %in% c("J0","J90","J180","J270")){

  c.copula.be1 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 
    1) * teta * (1 - p2)^teta))


  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))


  c.copula.thet <-    (((1 - p1)^
   teta - (-1 + (1 - p1)^teta) *(1 - p2)^
    teta)^(1/teta)* (log((1 - p1)^
     teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta) + (
   teta* ((1 - p1)^
       teta* (-1 + (1 - p2)^teta)* log(
        1 - p1) + (-1 + (1 - p1)^teta) *(1 - p2)^
       teta* log(1 - p2)))/((1 - p1)^
    teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)))/teta^2
        
  
c.copula2.be1 <- (1 - p1)^(-2 + 
  teta)* (-1 + (1 - p2)^teta)* ((1 - p1)^
   teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^teta* (-1 + teta)

       
c.copula2.be2 <- (1 - p2)^(-2 + 
  teta)* (-1 + (1 - p1)^teta)* ((1 - p2)^
   teta - (-1 + (1 - p2)^teta)* (1 - p1)^teta)^(-2 + 1/
  teta) *(1 - p1)^teta* (-1 + teta)


c.copula2.be1be2 <- (1 - p1)^(-1 + 
  teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^(-1 + 
  teta)* (-(-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) + teta)

      
 c.copula2.be1t <- ((((1 - p1)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - p1) + 
    ((1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - 
        p2)) * ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(1/teta - 2) * (1/teta - 1) - ((1 - p1)^teta + 
    (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1) * log((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)/teta^2) * ((1 - p1)^(teta - 1) - (1 - p1)^(teta - 
    1) * (1 - p2)^teta) + ((1 - p1)^(teta - 1) + (1 - p1)^(teta - 
    1) * (1 - p2)^teta + teta * ((1 - p1)^(teta - 1) * log(1 - 
    p1) - (1 - p1)^(teta - 1) * (1 - p2)^teta * log(1 - p2)) - 
    (((1 - p1)^(teta - 1) + teta * (1 - p1)^(teta - 1) * log(1 - 
        p1)) * (1 - p2)^teta + (1 - p1)^(teta - 1))) * ((1 - 
    p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1)/teta
  
  
 c.copula2.be2t <- ((((1 - p1)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - p1) + 
    ((1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta) * log(1 - 
        p2)) * ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(1/teta - 2) * (1/teta - 1) - ((1 - p1)^teta + 
    (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1) * log((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)/teta^2) * ((1 - p2)^(teta - 1) - (1 - p1)^teta * 
    (1 - p2)^(teta - 1)) + ((1 - p1)^teta * (1 - p2)^(teta - 
    1) + (1 - p2)^(teta - 1) + teta * ((1 - p2)^(teta - 1) * 
    log(1 - p2) - (1 - p1)^teta * (1 - p2)^(teta - 1) * log(1 - 
    p1)) - (((1 - p2)^(teta - 1) + teta * (1 - p2)^(teta - 1) * 
    log(1 - p2)) * (1 - p1)^teta + (1 - p2)^(teta - 1))) * ((1 - 
    p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta - 
    1)/teta

  
bit1.th2ATE <- ((1 - p1)^teta - ((1 - p1)^teta - 1) * (1 - p2)^teta)^(1/teta - 
    2) * ((2 * (teta * (((1 - p1)^teta - 1) * (1 - p2)^teta - 
    (1 - p1)^teta) * (((1 - p1)^teta - 1) * (1 - p2)^teta * (log(1 - 
    p2) - 1) + (1 - p1)^teta)) - ((1 - p1)^teta - ((1 - p1)^teta - 
    1) * (1 - p2)^teta)^2 * log((1 - p1)^teta - ((1 - p1)^teta - 
    1) * (1 - p2)^teta)) * log((1 - p1)^teta - ((1 - p1)^teta - 
    1) * (1 - p2)^teta) + teta * (2 * (((((1 - p1)^teta - 1) * 
    (1 - p2)^teta - (1 - p1)^teta) * ((1 - p2)^teta - 1) * log((1 - 
    p1)^teta - ((1 - p1)^teta - 1) * (1 - p2)^teta) + teta * 
    ((((1 - p1)^teta + teta - 1) * log(1 - p2) + 1 - 2 * (1 - 
        p1)^teta) * (1 - p2)^teta + (1 - p1)^teta - ((1 - p1)^teta - 
        1) * (1 - p2)^(2 * teta) * (log(1 - p2) - 1))) * (1 - 
    p1)^teta * log(1 - p1)) + teta * (((1 - p1)^teta - 1) * ((1 - 
    p1)^teta * (teta * log(1 - p2) - 2) - ((1 - p1)^teta - 1) * 
    (1 - p2)^teta * (log(1 - p2) - 2)) * (1 - p2)^teta * log(1 - 
    p2) - (((1 - p2)^teta - 1) * (1 - p1)^teta - teta * (1 - 
    p2)^teta) * ((1 - p2)^teta - 1) * (1 - p1)^teta * log(1 - 
    p1)^2)))/teta^4
 
}








if(BivD == "AMH"){

  c.copula.be1 <- p2 * (1 - p1 * teta * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))


  c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))


  c.copula.thet <- p1 * p2 * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))^2
        
  
c.copula2.be1 <- -(p2 * teta * (1 - p2) * (2 - 2 * (p1 * teta * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2)

       
c.copula2.be2 <- -(p1 * teta * (1 - p1) * (2 - 2 * (p2 * teta * (1 - p1)/(1 - 
    teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2)



c.copula2.be1be2 <- (1 - teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - p1) * (1 - 
    p2))

      
 c.copula2.be1t <- p2 * (1 - p1 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * 
    (1 - p1) * (1 - p2))))) * (1 - p2)/(1 - teta * (1 - p1) * 
    (1 - p2))^2

  
  
 c.copula2.be2t <- p1 * (1 - p1) * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2)))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2

  
bit1.th2ATE <- 2 * (p1 * p2 * (1 - p1)^2 * (1 - p2)^2/(1 - teta * (1 - p1) * 
    (1 - p2))^3)
 
}







if(BivD == "FGM"){

  c.copula.be1 <- p2 * (1 + teta * (1 - 2 * p1) * (1 - p2))


  c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))


  c.copula.thet <- p1 * p2 * (1 - p1) * (1 - p2)

        
  
c.copula2.be1 <- -(2 * (p2 * teta * (1 - p2)))


       
c.copula2.be2 <- -(2 * (p1 * teta * (1 - p1)))


c.copula2.be1be2 <- 1 + teta * (1 - 2 * p1) * (1 - 2 * p2)


      
 c.copula2.be1t <- p2 * (1 - 2 * p1) * (1 - p2)
  
  
 c.copula2.be2t <- p1 * (1 - 2 * p2) * (1 - p1)


  
bit1.th2ATE <- 0
 
}










#########################
# modular derivatives

c.copula.theta  <- c.copula.thet*derteta.derteta.st
c.copula2.be1th <- c.copula2.be1t*derteta.derteta.st
c.copula2.be2th <- c.copula2.be2t*derteta.derteta.st 


#if(BivD == "T"){
#
#c.copula.dof.st     <- c.copula.dof*derdof.derdof.st
#c.copula2.be1dof.st <- c.copula2.be1dof*derdof.derdof.st
#c.copula2.be2dof.st <- c.copula2.be2dof*derdof.derdof.st
#
#c.copula2.dof2.st <- c.copula2.dof2*derdof.derdof.st^2 + c.copula.dof*der2dof.derdof.stdof.st 
#
#c.copula2.thdof.st <- c.copula2.tetadof*derdof.derdof.st*derteta.derteta.st 
#
#}

                    
bit1.th2 <- bit1.th2ATE*derteta.derteta.st^2 + c.copula.thet*der2teta.derteta.stteta.st   

#########################





if(BivD %in% c("C90","J90","G90") ) {


c.copula.be2     <- 1 - c.copula.be2
c.copula.theta   <- - c.copula.theta
c.copula.thet    <- - c.copula.thet
c.copula2.be1    <- - c.copula2.be1
c.copula2.be2    <- - c.copula2.be2

c.copula2.be2th  <- - c.copula2.be2th
c.copula2.be2t   <- - c.copula2.be2t

bit1.th2ATE      <- - bit1.th2ATE  
bit1.th2         <- - bit1.th2 

} 



if(BivD %in% c("C180","J180","G180") ) {

c.copula.be1     <- 1 - c.copula.be1 
c.copula.be2     <- 1 - c.copula.be2

c.copula2.be1th  <- - c.copula2.be1th
c.copula2.be2th  <- - c.copula2.be2th

c.copula2.be1t  <- - c.copula2.be1t
c.copula2.be2t  <- - c.copula2.be2t
 

}  


if(BivD %in% c("C270","J270","G270") ) {

c.copula.be1     <- 1 - c.copula.be1

c.copula.theta   <- - c.copula.theta
c.copula.thet    <- - c.copula.thet
c.copula2.be1    <- - c.copula2.be1
c.copula2.be2    <- - c.copula2.be2

c.copula2.be1th  <- - c.copula2.be1th
c.copula2.be1t   <- - c.copula2.be1t

bit1.th2ATE      <- - bit1.th2ATE
bit1.th2         <- - bit1.th2

}   




  
  
if(CLM == FALSE){  
  
c.copula.be2 <- mm(c.copula.be2, min.pr = min.pr, max.pr = max.pr)
c.copula.be1 <- mm(c.copula.be1, min.pr = min.pr, max.pr = max.pr)
c.copula2.be1be2 <- ifelse(c.copula2.be1be2 < min.dn, min.dn, c.copula2.be1be2)


}



list(

c.copula.be1               = ifef(c.copula.be1    ) , 
c.copula.be2               = ifef(c.copula.be2    ) ,
c.copula.theta             = ifef(c.copula.theta  ) ,
c.copula.thet              = ifef(c.copula.thet   ) ,
c.copula2.be1              = ifef(c.copula2.be1   ) ,
c.copula2.be2              = ifef(c.copula2.be2   ) ,
c.copula2.be1be2           = ifef(c.copula2.be1be2) ,
c.copula2.be1th            = ifef(c.copula2.be1th ) ,
c.copula2.be2th            = ifef(c.copula2.be2th ) ,
c.copula2.be1t             = ifef(c.copula2.be1t  ) ,
c.copula2.be2t             = ifef(c.copula2.be2t  ) ,
bit1.th2ATE                = ifef(bit1.th2ATE     ) ,
bit1.th2                   = ifef(bit1.th2        ) ,
derteta.derteta.st         = ifef(derteta.derteta.st        ),
der2teta.derteta.stteta.st = ifef(der2teta.derteta.stteta.st),
c.copula.dof.st            = ifef(c.copula.dof.st    ) ,
c.copula2.be1dof.st        = ifef(c.copula2.be1dof.st) ,
c.copula2.be2dof.st        = ifef(c.copula2.be2dof.st) ,
c.copula2.dof2.st          = ifef(c.copula2.dof2.st),
c.copula2.thdof.st         = ifef(c.copula2.thdof.st)

)



}

