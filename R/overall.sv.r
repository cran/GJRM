overall.sv <- function(margins, M, vo, type = "copR", c.gam2 = NULL){
  


if(type == "copR"){


BivD <- M$BivD
if(M$surv == TRUE) BivD <- "N" # to make sure that dof is not estimated in surv context 


if(BivD == "T" && margins[1] %in% c(M$m2,M$m3) && margins[2] %in% c(M$m2,M$m3) ){ 

if( margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2) ) start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2,                           vo$dof.st, vo$i.rho)                         
if( margins[1] %in% c(M$m3) && margins[2] %in% c(M$m3) ) start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2, vo$log.nu.1, vo$log.nu.2, vo$dof.st, vo$i.rho)                           
if( margins[1] %in% c(M$m2) && margins[2] %in% c(M$m3) ) start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2,              vo$log.nu.2, vo$dof.st, vo$i.rho) 
if( margins[1] %in% c(M$m3) && margins[2] %in% c(M$m2) ) start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2, vo$log.nu.1,              vo$dof.st, vo$i.rho) 

}else{


if( margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2,                            vo$i.rho)                         
if( margins[1] %in% c(M$m3) && margins[2] %in% c(M$m3) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2, vo$log.nu.1, vo$log.nu.2,  vo$i.rho)                           
if( margins[1] %in% c(M$m2) && margins[2] %in% c(M$m3) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2,              vo$log.nu.2,  vo$i.rho) 
if( margins[1] %in% c(M$m3) && margins[2] %in% c(M$m2) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2, vo$log.nu.1,               vo$i.rho) 


if( margins[1] %in% c(M$m2) && margins[2] %in% c(M$bl) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1,               vo$i.rho)                         
if( margins[1] %in% c(M$m3) && margins[2] %in% c(M$bl) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1, vo$log.nu.1,  vo$i.rho)                           



if(margins[1] %in% c(M$m1d,M$bl) && margins[2] %in% c(M$m1d,M$bl))   start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients, vo$i.rho) 
if(margins[1] %in% c(M$m1d) && margins[2] %in% c(M$m2,M$m2d))        start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients,                vo$log.sig2.2,                                      vo$i.rho) 
if(margins[1] %in% c(M$m1d) && margins[2] %in% c(M$m3))              start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients,                vo$log.sig2.2,              vo$log.nu.2,            vo$i.rho)                        
if(margins[1] %in% c(M$m2d) && margins[2] %in% c(M$m2d) )            start.v <- c(vo$gam1$coefficients, vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2,                                      vo$i.rho)                         


}




 
}


if(type == "copSS"){


if(margins[2] %in% c(M$m1d))       start.v <- c(vo$gam1$coefficients, c.gam2, vo$i.rho) 
if(margins[2] %in% c(M$m2,M$m2d))  start.v <- c(vo$gam1$coefficients, c.gam2, vo$log.sig2.2, vo$i.rho)                            
if(margins[2] %in% M$m3)           start.v <- c(vo$gam1$coefficients, c.gam2, vo$log.sig2.2, vo$log.nu.2, vo$i.rho) 
      
 

 
}






        			
start.v

}

