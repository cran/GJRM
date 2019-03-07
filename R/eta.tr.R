eta.tr <- function(vrb.st, margin){

mupos  <- c("LN","WEI","GO","iG","GA","GA2","GGA","DAGUM","SM","FISK","NBI","NBII","PO","ZTP","PIG")
mub    <- c("BE","logit")
muNpos <- c("GP","DGP")


if(margin %in% mupos){
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st )  # it was 28
   vrb.st <- ifelse( vrb.st < -14.5, -14.5, vrb.st ) # it was -17 then -9
}


if(margin %in% mub){
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -16, -16, vrb.st ) 
}


if(margin == "probit"){
   vrb.st <- ifelse( vrb.st >  5.3,  5.3, vrb.st )  
   vrb.st <- ifelse( vrb.st < -5.3, -5.3, vrb.st ) 
}

if(margin == "cloglog"){
   vrb.st <- ifelse( vrb.st >  2.9,  2.9, vrb.st )  
   vrb.st <- ifelse( vrb.st <  -16,  -16, vrb.st ) 
}

# stuff for probit, logit and cloglog not in use

vrb.st  
 
} 

