enu.tr <- function(vrb.st, margin){
 
if( margin %in% c("SICHEL") ){
 
   vrb.st <- ifelse( vrb.st > 125,   125, vrb.st )  
   vrb <- vrb.st <- ifelse( vrb.st < -125, -125, vrb.st ) 
   
}


if( margin %in% c("DEL") ){
 
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -16, -16, vrb.st ) 
   vrb    <- plogis(vrb.st)
    
}



if( !(margin %in% c("BE","DEL","SICHEL")) ){
 
   vrb.st <- ifelse( vrb.st > 28,   28, vrb.st )  
   vrb.st <- ifelse( vrb.st < -17, -17, vrb.st ) 
   vrb    <- exp(vrb.st)
    
}


if( margin %in% c("BE") ){
 
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -16, -16, vrb.st ) 
   vrb    <- plogis(vrb.st)
    
}


    
    
 list(vrb = vrb, vrb.st = vrb.st )  
 
}    