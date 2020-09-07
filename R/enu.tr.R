enu.tr <- function(vrb.st, margin){
 
 

#if( !(margin %in% c("BE")) ){ 
# tw fine here since this is for phi (or nu)
 
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st ) # it was 20  # 485165195
   vrb.st <- ifelse( vrb.st < -13, -13, vrb.st ) # it was -17 # 2.260329e-06
   vrb    <- exp(vrb.st)
    



    
    
 list(vrb = vrb, vrb.st = vrb.st )  
 
}    