numgh <- function(funcD, para){

para <- c(para)



#
# Test versions to check/compare
#
#
#
# the solution below is faster than the other option below
# but slightly less accurate but by a very very small amount
#
#eps <-  1e-06 # sqrt(.Machine$double.eps) # 1e-06  #  #  # 
#eps2 <- eps*eps 
#
#para1 <- para - eps/2
#para2 <- para + eps/2 
#f1 <- funcD(para1)
#f2 <- funcD(para2)
#
#fi <- (f2 - f1)/ eps  
#
#t01 <- t10 <- t11 <- para + eps
#t11 <- t11 + eps 
#
#f00 <- funcD(para) 
#f01 <- funcD(t01) 
#f10 <- funcD(t10) 
#f11 <- funcD(t11) 
#
#se <- (f11 - f01 - f10 + f00)/eps2 


if(length(para) == 1){

fi <- jacobian(function(para) funcD(para), para)
se <- jacobian(function(para) jacobian(function(para) funcD(para), para), para)

                     }

if(length(para) > 1){

fi <- grad(function(para) funcD(para), para)
se <- grad(function(para) grad(function(para) funcD(para), para), para)

                    }

list(fi = fi, se = se)


}
    
  