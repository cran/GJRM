esp.tr <- function(vrb.st, margin){
 
 
mub <- c("BE", "TW") 
 

if(   !(margin %in% mub)    ){
 
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st )  # it was 28
   vrb.st <- ifelse( vrb.st < -13, -13, vrb.st )  # it was -17 
   
   vrb    <- exp(vrb.st)
    
}


if( margin %in% c("BE") ){
 
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -18, -18, vrb.st ) 
   vrb    <- plogis(vrb.st)
    
}



if( margin %in% c("TW") ){
 
 
 
aTW <- 1.001  
bTW <- 1.999

# sigma.stTW <- log( (sigma - aTW) / (bTW - sigma) ) # sigma = (aTW + bTW*exp(sigma.stTW))/(1 + exp(sigma.stTW))
 
   vrb.st <- ifelse( vrb.st >  15,  15, vrb.st )  
   vrb.st <- ifelse( vrb.st < -15, -15, vrb.st ) 
   vrb    <- (aTW + bTW*exp(vrb.st))/(1 + exp(vrb.st)) 
    
}


    
 list(vrb = vrb, vrb.st = vrb.st )  
 
}    