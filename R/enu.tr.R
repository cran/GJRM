enu.tr <- function(vrb.st, margin){
 
if( margin %in% c("SICHEL") ){ # not used
 
   vrb.st <- ifelse( vrb.st > 125,   125, vrb.st )  
   vrb <- vrb.st <- ifelse( vrb.st < -125, -125, vrb.st ) 
   
}


if( margin %in% c("DEL") ){ # not used
 
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -18, -18, vrb.st ) 
   vrb    <- plogis(vrb.st)
    
}



if( !(margin %in% c("BE","DEL","SICHEL")) ){ # tw fine here since this is for phi
 
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st ) # it was 20  
   vrb.st <- ifelse( vrb.st < -13, -13, vrb.st ) # it was -17 
   vrb    <- exp(vrb.st)
    
}


if( margin %in% c("BE") ){
 
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -18, -18, vrb.st ) 
   vrb    <- plogis(vrb.st)
    
}


    
    
 list(vrb = vrb, vrb.st = vrb.st )  
 
}    