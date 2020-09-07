copgHs3 <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD, par2 = NULL){




########################################################################################

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







if(BivD == "HO"){




 c.copula2.be2 <- ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(2 * 
    (teta - 1)) - ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta + ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((-log(p2))^(1/teta - 1) + 
    (-log(p2))^(1/teta - 2) * (1/teta - 1)))) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)/p2^2 

der2h.derp2p2 <- ((-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(2 * (teta - 1)) - ((-log(p2))^(2 * 
    (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * (teta - 1)/teta + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 2) * 
    (1/teta - 1)))) + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * (((-log(p2))^(1/teta - 2) + (-log(p2))^(1/teta - 3) * 
    (1/teta - 2) + 2 * (-log(p2))^(1/teta - 2)) * (1/teta - 1) + 
    2 * (-log(p2))^(1/teta - 1)) + ((-log(p2))^(1/teta - 1) * 
    ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        3) * (teta - 2)/teta + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        2) * ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 
        2) * (1/teta - 1))) + 2 * ((-log(p2))^(1/teta - 1) * 
    ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 2) * 
    ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 2) * (1/teta - 
        1)))) * (teta - 1)/teta - 2 * ((-log(p2))^(1/teta - 1) * 
    ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 1) * 
    ((-log(p2))^(2 * (1/teta - 1)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        2) * (teta - 1)/teta + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1) * ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 
        2) * (1/teta - 1))))) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2^3 


}







if(BivD == "PL"){




 c.copula2.be2 <- -((2/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1))) - 0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 
    4 * (p1 * teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * 
    (p1 * p2 * teta * (teta - 1)))^1.5)) * (teta - 1)/4) 

der2h.derp2p2 <- (0.5 * (2 - 1.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * 
    teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * 
    teta * (teta - 1))))) + 2) * (2 * ((p1 + p2) * (teta - 1) + 
    1) - 4 * (p1 * teta)) * (teta - 1)^2/(4 * (((p1 + p2) * (teta - 
    1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5) 


}










if(BivD=="N"){




c.copula2.be2 <- dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p2)^2/2) 

der2h.derp2p2 <-  -(teta * dnorm((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - 
    teta^2)) * (qnorm(p2) + teta * (qnorm(p1) - teta * qnorm(p2))/(1 - 
    teta^2))/(dnorm(qnorm(p2))^2 * sqrt(1 - teta^2)))


}



if(BivD == "T"){
                    
c.copula2.be2 <-  BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=par2, deriv="u2") 

der2h.derp2p2 <- BiCopHfuncDeriv2(p1, p2, family = 2, teta, par2, deriv = "u2") 



}





if(BivD=="F"){

 c.copula2.be2 <-     (exp(teta + 
  p2* teta)* (-1 + exp(p1* teta))* (-exp(teta) + exp(
   p1* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2



t1 = exp(teta)
t2 = teta*p1
t3 = exp(t2)
t5 = t1*(t3-1)
t6 = teta*p2
t8 = exp(t6+t2)
t10 = exp(t6+teta)
t12 = exp(t2+teta)
t13 = t8-t10-t12+t1
t14 = t13*t13
t20 = (teta*t8-teta*t10)^2
t24 = teta^2

       der2h.derp2p2 <- -2.0*t5/t14/t13*t20+t5/t14*(t24*t8-t24*t10)

}





if(BivD %in% c("C0","C90","C180","C270")){


 

 c.copula2.be2 <- (p2^(-2 + teta)* p1^teta* (-1 + p2^-teta + p1^-teta)^(-1/
  teta)* (-1 + p1^teta)* (1 + teta))/(p1^teta - 
  p2^teta* (-1 + p1^teta))^2


t1 = -teta-1
t2 = p2^t1
t3 = t1^2
t5 = p2^2
t6 = 1/t5
t7 = p1^-teta
t8 = p2^-teta
t9 = t7+t8-1
t11 = -1-1/teta
t12 = t9^t11
t13 = t6*t12
t16 = t2*t1*t13
t18 = 1/t9
t23 = t2*t12
t24 = t11*t11
t26 = t8*t8
t27 = teta^2
t29 = t9*t9
t32 = t26*t27*t6/t29
t34 = t23*t11
t36 = t6*t18

    der2h.derp2p2 <- t2*t3*t13-t16-2*t16*t11*t8*teta*t18+t23*t24*t32+t34*t8*t27*t36+t34*t8*teta*t36-t34*t32
   
    
}








if(BivD %in% c("G0","G90","G180","G270")){



                
c.copula2.be2 <- (1/(p2^2))*exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/
  teta)))* (-log(p2))^(-2 + 
  teta) *((-log(p2))^teta + (-log(p1))^teta)^(-2 + 1/
  teta)* ((-log(p2))^
    teta* (log(p2) + ((-log(p2))^teta + (-log(p1))^teta)^(1/
      teta)) + (1 - teta + log(p2))* (-log(p1))^teta)


  
t1 = log(p2)
t2 = (-t1)^teta 
t3 = log(p1)
t4 = (-t3)^teta  
t5 = t2+t4
t6 = 1/teta
t7 = t5^t6
t8 = exp(-t7)
t9 = t6-1
t10 = t5^t9
t11 = t8*t10
t12 = p2^2
t14 = 1/t12/p2
t15 = t2*t14
t20 = t1*t1
t21 = 1/t20
t26 = 1/t20/t1
t30 = t7*t7
t31 = t2*t2
t32 = t31*t2
t34 = t5*t5
t36 = 1/t34
t37 = t26*t36
t38 = t37*t11
t40 = t11*t2
t41 = teta*t14
t48 = t7*t31
t49 = t48*t14
t50 = 1/t5
t51 = t21*t50
t55 = t26*t50
t59 = t7*t32
t62 = t14*t26
t65 = t10*teta
t69 = t59*t62
t70 = t36*t8
t79 = t11*t9*t31
t86 = t9*t9
t18 = teta^2
t16 = t18*t14
t13 = t16*t37
 

 
der2h.derp2p2 <- -2*t11*t15/t1-3*t11*t15*t21-2*t11*t15*t26-t30*t32*t14*t38+3*t40*t41*t21+3*t40*t41*t26-3*t49*t51*t11-3*t49*t55*t11+t59*t14*t38+3*t48*t62*t50*t8*t65-t69*t70*t65+2*t69*t70*t10*t9*teta+3*t79*t41*t51+3*t79*t41*t55-t11*t86*t32*t13-3*t79*t16*t55+t11*t9*t32*t13-t40*t16*t26 
 

}






if(BivD %in% c("J0","J90","J180","J270")){


       
c.copula2.be2 <- (1 - p2)^(-2 + 
  teta)* (-1 + (1 - p1)^teta)* ((1 - p2)^
   teta - (-1 + (1 - p2)^teta)* (1 - p1)^teta)^(-2 + 1/
  teta) *(1 - p1)^teta* (-1 + teta)

			t2 = (1 - p1)^teta
			t3 = 1 - p2
			t4 = t3^teta
			t5 = t2*t4
			t6 = t2+t4-t5
			t8 = 1/teta - 1
			t9 = t6^t8
			t10 = t8^2
			t12 = t4*teta
			t13 = 1/t3
			t14 = -t12*t13+t5*teta*t13
			t18 = t14^2
			t20 = t6^2
			t22 = teta - 1
			t23 = t3^t22
			t24 = 1 - t2
			t26 = 1/t20*t23*t24
			t27 = t9*t8
			t29 = teta^2
			t31 = t3*t3
			t32 = 1/t31
			t41 = 1/t6
			t51 = t9*t23
			t55 = t22*t22
			
der2h.derp2p2 <- t9*t10*t18*t26+t27*(t4*t29*t32-t12*t32-t5*t29*t32+t5*teta*t32)*t41*t23*t24-t27*t18*t26-2.0*t27*t14*t41*t23*t22*t13*t24+t51*t55*t32*t24-t51*t22*t32*t24 
 
 
}








if(BivD == "AMH"){


       
c.copula2.be2 <- -(p1 * teta * (1 - p1) * (2 - 2 * (p2 * teta * (1 - p1)/(1 - 
    teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - p1) * (1 - 
    p2))^2)

der2h.derp2p2 <- p1 * teta^2 * (1 - p1)^2 * (2 * (1 - p2 * teta * (1 - p1)/(1 - 
    teta * (1 - p1) * (1 - p2))) + 2 * (2 - 2 * (p2 * teta * 
    (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))))/(1 - teta * 
    (1 - p1) * (1 - p2))^3

 
}







if(BivD == "FGM"){


       
c.copula2.be2 <- -(2 * (p1 * teta * (1 - p1)))

der2h.derp2p2 <- 0
 
}










if(BivD %in% c("C90","J90","G90") ) {



c.copula2.be2  <- - c.copula2.be2
der2h.derp2p2  <- -der2h.derp2p2




} 



if(BivD %in% c("C180","J180","G180") ) {

der2h.derp2p2 <- -der2h.derp2p2

 

}  


if(BivD %in% c("C270","J270","G270") ) {


c.copula2.be2  <- - c.copula2.be2



}   






list( c.copula2.be2 = ifef(c.copula2.be2), der2h.derp2p2 = ifef(der2h.derp2p2) )     


}

