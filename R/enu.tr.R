enu.tr <- function(vrb.st, margin){
 
 

if( !(margin %in% c("TW")) ){
 
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st ) # it was 20  # 485165195
   vrb.st <- ifelse( vrb.st < -13, -13, vrb.st ) # it was -17 # 2.260329e-06
   vrb    <- exp(vrb.st)
   
}   
   
   
   
if( margin %in% c("TW") ){
 
  aTW <- 1.001  
  bTW <- 1.999
 
   vrb.st <- ifelse( vrb.st >  15,  15, vrb.st )  
   vrb.st <- ifelse( vrb.st < -15, -15, vrb.st ) 
   vrb    <- (aTW + bTW*exp(vrb.st))/(1 + exp(vrb.st)) 
  
}   
   
   
    
    
 list(vrb = vrb, vrb.st = vrb.st )  
 
}    