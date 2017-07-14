numgh <- function(funcD, para){

#fd.prec  <- 10^(-7)
#fd.prec2 <- 10^-3
#
#para1 <- para + fd.prec 
#
#f1 <- funcD(para)
#f2 <- funcD(para1)
#
#fi <- (f2 - f1)/fd.prec  
##
##
##para1 <- para + fd.prec2 
#
#
#f1 <- funcD(para1)
#f2 <- funcD(para1+fd.prec)
#    
#se <- ( (((f2 - f1) / fd.prec)*para1) - fi ) / fd.prec2 
#
#


eps <-  1e-06 # sqrt(.Machine$double.eps) # 1e-06  #  #  # 
eps2 <- eps*eps 

para1 <- para - eps/2
para2 <- para + eps/2 
f1 <- funcD(para1)
f2 <- funcD(para2)

fi <- (f2 - f1)/ eps  

t01 <- t10 <- t11 <- para + eps
t11 <- t11 + eps 

f00 <- funcD(para) 
f01 <- funcD(t01) 
f10 <- funcD(t10) 
f11 <- funcD(t11) 

se <- (f11 - f01 - f10 + f00)/eps2 

list(fi = fi, se = se)


}
    
  