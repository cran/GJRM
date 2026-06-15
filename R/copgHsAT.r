copgHsAT <- function(p1, p2, teta, BivD, Ln = FALSE, par2 = NULL,
                     min.dn, min.pr, max.pr, teta.st = 1, T.BivD = "F"){



c.copula.be2 <- c.copula2.be1be2 <- c.copula.be1 <- c.copula2.be2th <- 1


BivD2mixed <- c("C0C90","C0C270","C180C90","C180C270", "J0J90","J0J270","J180J90","J180J270", "G0G90","G0G270","G180G90","G180G270", "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270") 


cjg <- c("C0","J0","G0","C90","J90","G90","C180","J180","G180","C270","J270","G270","PL",
         "GAL0","GAL90","GAL180","GAL270")

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



if(BivD %in% c("C0", "C180","GAL0","GAL180","J0", "J180","G0", "G180")  && T.BivD %in% BivD2mixed){ derteta.derteta.st         <- 1 # der teta.st wrt teta.st  
                                                                                                    der2teta.derteta.stteta.st <- 0} 

if(BivD %in% c("C90","C270","GAL90","GAL270","J90","J270","G90","G270") && T.BivD %in% BivD2mixed){ derteta.derteta.st         <- -1 #  
                                                                                                    der2teta.derteta.stteta.st <- 0}










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




if(Ln == FALSE){





if(BivD == "PL"){
                    
c.copula.be2 <- (teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * 
    teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * (teta - 1)) 

}





if(BivD == "HO"){
                    
c.copula.be2 <- (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2 

}








if(BivD == "AMH"){
                    
c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))

}

if(BivD == "FGM"){
                    
c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))
}




if(BivD == "N"){
                    
c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   )  

}

if(BivD == "T"){
                    
c.copula.be2 <- BiCopHfunc2(p1, p2, family = 2, par = teta, par2 = par2) 


}


if(BivD == "F"){

c.copula.be2 <- (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + 
                     exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))
    
}


if(BivD %in% c("C0","C90","C180","C270") ){

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))

}



if(BivD %in% c("GAL0","GAL90","GAL180","GAL270") ){

c.copula.be2 <- p1*(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))




}




if(BivD %in% c("G0","G90","G180","G270") ){

c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2

}


if(BivD %in% c("J0","J90","J180","J270") ){

  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))

}



if(BivD %in% c("C90","J90","G90","GAL90") )     c.copula.be2  <- 1 - c.copula.be2
if(BivD %in% c("C180","J180","G180","GAL180") ) c.copula.be2  <- 1 - c.copula.be2


c.copula.be2 <- mm(c.copula.be2, min.pr = min.pr, max.pr = max.pr)



}







if(Ln == TRUE){





if(BivD == "HO"){
                    
c.copula2.be1be2 <- ((-log(p1))^(1/teta - 1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(2 * (teta - 1)) - (-log(p1))^(1/teta - 
    1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * (teta - 1)/teta) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/(p1 * 
    p2) 
    

c.copula.be1 <- (-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p1 

c.copula.be2 <- (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2 
    
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
      
      

}





if(BivD == "PL"){
                    
c.copula2.be1be2 <- -(((2 * (teta - 1) - 4 * teta)/sqrt(((p1 + p2) * (teta - 1) + 
    1)^2 - 4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * ((p1 + 
    p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p2 * teta)) * (teta - 1)/(((p1 + p2) * 
    (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5))/4) 



c.copula.be1 <- (teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p2 * 
    teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * (teta - 1)) 

c.copula.be2 <- (teta - (0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * 
    teta)) * (teta - 1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 1))/(2 * (teta - 1)) 

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
    

}









if(BivD == "AMH"){
                    
c.copula2.be1be2 <- (1 - teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - p1) * (1 - 
    p2))
    

  c.copula.be1 <- p2 * (1 - p1 * teta * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))


  c.copula.be2 <- p1 * (1 - p2 * teta * (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
                  teta * (1 - p1) * (1 - p2))

    
 c.copula2.be2t <- p1 * (1 - p1) * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2)))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2


}

if(BivD == "FGM"){
                    
c.copula2.be1be2 <- 1 + teta * (1 - 2 * p1) * (1 - 2 * p2)


  c.copula.be1 <- p2 * (1 + teta * (1 - 2 * p1) * (1 - p2))


  c.copula.be2 <- p1 * (1 + teta * (1 - 2 * p2) * (1 - p1))

 c.copula2.be2t <- p1 * (1 - 2 * p2) * (1 - p1)


}






if(BivD == "N"){
                    
c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 



c.copula.be1 <- pnorm( (qnorm(p2) - teta*qnorm(p1))/sqrt(1 - teta^2)   )                            
c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   )  

c.copula2.be2t <- (-(dnorm((qnorm(p1) - tanh(teta.st) * qnorm(p2))/sqrt(1 - tanh(teta.st)^2)) *  # this as well
    (1/cosh(teta.st)^2 * qnorm(p2)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p1) - 
        tanh(teta.st) * qnorm(p2)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
        tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
        tanh(teta.st)^2)^2)))/derteta.derteta.st


}



if(BivD == "T"){
                    
c.copula2.be1be2 <- BiCopPDF(p1, p2, family = 2, teta, par2) 


c.copula.be1 <- BiCopHfunc1(p1, p2, family = 2, par = teta, par2 = par2)                         
c.copula.be2 <- BiCopHfunc2(p1, p2, family = 2, par = teta, par2 = par2) 

c.copula2.be2t <- BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=par2, deriv="par")

}




if(BivD == "F"){

c.copula2.be1be2 <- (exp((1 + p1 + p2)* teta)* (-1 + exp(teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2
 

  c.copula.be1 <- (exp(teta)* (-1 + exp(p2* teta)))/(-exp((p1 + p2)* teta) + exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))

  c.copula.be2 <- (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))
 
c.copula2.be2t <-  (exp(teta + 
  p2 *teta)* (exp(2* p1* teta)* (-1 + p2) + exp(teta)* p2 - 
   exp(p1* teta)* (-1 + p2 + exp(teta)* (p2 - p1) + p1)))/(exp((p2 + 
     p1)* teta) - exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2
 
 
}


if(BivD %in% c("C0","C90","C180","C270") ){

c.copula2.be1be2 <- p1^(-1 - teta)* p2^(-1 - teta)* (-1 + p1^-teta + p2^-teta)^(-2 - 1/
  teta) *(1 + teta)



c.copula.be1 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta)))

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  
 c.copula2.be2t <- ((1 + 1/teta) * (log(p1)/p1^teta + log(p2)/p2^teta)/(1/p1^teta + 
    1/p2^teta - 1)^(1/teta + 2) + log(1/p1^teta + 1/p2^teta - 
    1)/(teta^2 * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta)))/p2^(1 + 
    teta) + (1/p2^(1 + teta) - (1/p2^(1 + teta) + teta * log(p2)/p2^(1 + 
    teta)))/(teta * (1/p1^teta + 1/p2^teta - 1)^(1 + 1/teta))
    
  
}


if(BivD %in% c("GAL0","GAL90","GAL180","GAL270") ){

c.copula2.be1be2 <- (1+teta*(-log(p2))^(1+teta)*(1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2)-((1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))/((-log(p1))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))+1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))))*
exp((1/(-log(p1))^teta+1/(-log(p2))^teta)^-(1/teta))


c.copula.be1 <- p2*(1-1/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))


c.copula.be2 <- p1*(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))

c.copula2.be2t <- p1*(((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p2))-(-log(p2))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2))/((-log(p2))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2+
((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)))*
(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))/teta)*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))




  
}



if(BivD %in% c("G0","G90","G180","G270") ){

c.copula2.be1be2 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta) *(-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
   teta))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-log(p2))^(-1 + teta))/(p1 *p2)


  c.copula.be1 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta)* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/teta))/p1


  c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2

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



}


if(BivD %in% c("J0","J90","J180","J270") ){

c.copula2.be1be2 <- (1 - p1)^(-1 + 
  teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^(-1 + 
  teta)* (-(-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) + teta)


  c.copula.be1 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 
    1) * teta * (1 - p2)^teta))


  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))

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

 

}




####


c.copula2.be2th <- c.copula2.be2t*derteta.derteta.st 



if(BivD %in% c("C90","J90","G90","GAL90") ) {
   c.copula.be2     <- 1 - c.copula.be2
   c.copula2.be2th  <- - c.copula2.be2th

} 


if(BivD %in% c("C180","J180","G180","GAL180") ) {
  c.copula.be1     <- 1 - c.copula.be1 
  c.copula.be2     <- 1 - c.copula.be2
  c.copula2.be2th  <- - c.copula2.be2th

}  


if(BivD %in% c("C270","J270","G270","GAL270") ) {
  c.copula.be1     <- 1 - c.copula.be1
}   


c.copula.be1     <- ifef(mm(c.copula.be1, min.pr = min.pr, max.pr = max.pr))
c.copula.be2     <- ifef(mm(c.copula.be2, min.pr = min.pr, max.pr = max.pr))
c.copula2.be1be2 <- ifef(ifelse(c.copula2.be1be2 < min.dn, min.dn, c.copula2.be1be2))
c.copula2.be2th  <- ifef(c.copula2.be2th ) 



  
}



list(c.copula.be1 = c.copula.be1, c.copula.be2 = c.copula.be2, c.copula2.be1be2 = c.copula2.be1be2, c.copula2.be2th = c.copula2.be2th)     


}




     























