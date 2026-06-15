copgHsCountFM <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD, nu = NULL, nu.st = NULL, CLM = FALSE,
                   min.dn, min.pr, max.pr, T.BivD = NULL){

########################################################################################

BivD2mixed <- c("C0C90","C0C270","C180C90","C180C270", "J0J90","J0J270","J180J90","J180J270", "G0G90","G0G270","G180G90","G180G270", "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270") 


c.copula.dof.st <- c.copula2.be1dof.st <- c.copula2.be2dof.st <- c.copula2.dof2.st <- c.copula2.thdof.st <- 1


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







###################################################################################
     
#if(BivD %in% c("T") ) derdof.derdof.st <- der2dof.derdof.stdof.st <- exp(nu.st)

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


c.copula.thet <- p1*p2*((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)))*
exp((1/(-log(p1))^teta+1/(-log(p2))^teta)^-(1/teta))/teta



bit1.th2ATE <- p1*p2*(((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)))*
((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta))-
1)/teta-(((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)-
((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta))*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/(teta*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta))^2+(1/(teta*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))-((1+
1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2)/(1/(-log(p1))^teta+1/(-log(p2))^teta)^(2*
(1+1/teta)))*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(log(-log(p1))^2/(-log(p1))^teta+log(-log(p2))^2/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))/teta




}








if(BivD=="HO"){
  
  c.copula.thet <- -((((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1)/teta) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)) 
  
 
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


  c.copula.thet <- (p1 + p2 - 0.5 * ((2 * (((p1 + p2) * (teta - 1) + 1) * (p1 + 
    p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta))/sqrt(((p1 + 
    p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))))/(2 * 
    (teta - 1)) - 2 * (((p1 + p2) * (teta - 1) + 1 - sqrt(((p1 + 
    p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1))))/(2 * 
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


c.copula.thet <- dbinorm(qnorm(p1),qnorm(p2), cov12=teta)

        
bit1.th2ATE <- (0.5 * (pi * teta/(pi * sqrt(1 - teta^2))^2) - 0.5 * 
    ((teta * (qnorm(p1) * (qnorm(p1) - 2 * (teta * qnorm(p2))) + 
        qnorm(p2)^2)/(1 - teta^2) - qnorm(p1) * qnorm(p2))/(pi * 
        (1 - teta^2)))) * exp(-(0.5 * ((qnorm(p1) * (qnorm(p1) - 
    2 * (teta * qnorm(p2))) + qnorm(p2)^2)/(1 - teta^2))))/sqrt(1 - 
    teta^2)


}


if(BivD=="T"){
 

c.copula.thet <- ( 1 + (qt(p1,nu)^2+qt(p2,nu)^2-2*teta*qt(p1,nu)*qt(p2,nu))/(nu*(1-teta^2)))^(-nu/2)/(2*pi*sqrt(1-teta^2))

                                                                                                                                                                                                                 
bit1.th2ATE <- -((1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 2 * teta * qt(p1,nu) * qt(p2,nu))/(nu * (1 - teta^2)))^((-nu/2) - 
    1) * ((-nu/2) * (2 * qt(p1,nu) * qt(p2,nu)/(nu * (1 - teta^2)) - (qt(p1,nu)^2 + 
    qt(p2,nu)^2 - 2 * teta * qt(p1,nu) * qt(p2,nu)) * (nu * (2 * teta))/(nu * (1 - 
    teta^2))^2))/(2 * pi * sqrt(1 - teta^2)) - (1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 
    2 * teta * qt(p1,nu) * qt(p2,nu))/(nu * (1 - teta^2)))^(-nu/2) * (2 * pi * 
    (0.5 * (2 * teta * (1 - teta^2)^-0.5)))/(2 * pi * sqrt(1 - 
    teta^2))^2)       


}








if(BivD=="F"){

  c.copula.thet <-  (exp(teta)* (1/(-1 + exp(teta)) + (-1 - exp(p2* teta)* (-1 + p1) + p1 - 
     exp(p1* teta)* (-1 + p2) + p2)/(
    exp((p1 + p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))* teta + 
 log((exp(-(p1 + p2)* teta)* (-exp((p1 + p2)* teta) + 
     exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))/(-1 + exp(teta))))/teta^2 


bit1.th2ATE <- 1/teta^2 * ((1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * 
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


  c.copula.thet <- ((-1 + p1^-teta + p2^-teta)^(-1/
  teta) *((teta *(p2^teta *log(p1) + p1^teta* log(p2)))/(
   p2^teta - p1^teta* (-1 + p2^teta)) + 
   log(-1 + p1^-teta + p2^-teta)))/teta^2
  

  
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


  c.copula.thet <-   (1/(teta^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/
  teta)* (-teta* (-log(p1))^
    teta* log(-log(p1)) + ((-log(p1))^teta + (-log(p2))^
      teta)* log((-log(p1))^teta + (-log(p2))^teta) - 
   teta *(-log(p2))^teta* log(-log(p2)))


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


  c.copula.thet <-    (((1 - p1)^
   teta - (-1 + (1 - p1)^teta) *(1 - p2)^
    teta)^(1/teta)* (log((1 - p1)^
     teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta) + (
   teta* ((1 - p1)^
       teta* (-1 + (1 - p2)^teta)* log(
        1 - p1) + (-1 + (1 - p1)^teta) *(1 - p2)^
       teta* log(1 - p2)))/((1 - p1)^
    teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)))/teta^2
        
  

  
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

  c.copula.thet <- p1 * p2 * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))^2
        
  

bit1.th2ATE <- 2 * (p1 * p2 * (1 - p1)^2 * (1 - p2)^2/(1 - teta * (1 - p1) * 
    (1 - p2))^3)
 
}







if(BivD == "FGM"){

  c.copula.thet <- p1 * p2 * (1 - p1) * (1 - p2)

        
  
bit1.th2ATE <- 0
 
}











if(BivD %in% c("C90","J90","G90","GAL90") ) {


c.copula.thet    <- - c.copula.thet
bit1.th2ATE      <- - bit1.th2ATE  

} 




if(BivD %in% c("C270","J270","G270","GAL270") ) {

c.copula.thet    <- - c.copula.thet
bit1.th2ATE      <- - bit1.th2ATE

}   



list( c.copula.thet = ifef(c.copula.thet), bit1.th2ATE = ifef(bit1.th2ATE), derteta.derteta.st = ifef(derteta.derteta.st),
      der2teta.derteta.stteta.st = ifef(der2teta.derteta.stteta.st) )

}

