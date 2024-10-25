eta.tr <- function(vrb.st, margin, zero.tol = 1e-02){




mupos  <- c("LN","WEI","GO","IG","GA","GA2","GGA","DAGUM","SM","FISK","TW","NBI","NBII","tNBI","tNBII","P","tP","PIG","tPIG","DGP0")
mub    <- c("BE","logit")
# muNpos <- c("GP","DGP") # nothing done here for now


if(margin %in% mupos){
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st ) # it was 28  # 485165195 # to prevent eta1 to space in a too wide range
   vrb.st <- ifelse( vrb.st < -13, -13, vrb.st ) # it was -17 # 2.260329e-06
}


if(margin %in% mub){
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -18, -18, vrb.st ) 
}

if(margin %in% c("DGPII")){
   vrb.st <- ifelse( vrb.st < -745.2, -745.2, vrb.st ) # to include 0
}


#if(margin == "probit"){
#   vrb.st <- ifelse( vrb.st >  5.3,  5.3, vrb.st )  
#   vrb.st <- ifelse( vrb.st < -5.3, -5.3, vrb.st )  
#}
#
#if(margin == "cloglog"){
#   vrb.st <- ifelse( vrb.st >  2.9,  2.9, vrb.st )  
#   vrb.st <- ifelse( vrb.st <  -16,  -16, vrb.st )  
#}
#
#
# stuff for probit, logit and cloglog not in use

vrb.st 
 
} 

