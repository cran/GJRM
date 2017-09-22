susutsn <- function(object, bs, lf, cont1par, cont2par, cont3par, prob.lev, type = "copR", bin.link = NULL, n.sim = NULL ){


CIrs <- CIkt <- CIsig21 <- CIsig22 <- CInu1 <- CInu2 <- CIdof <- CInu <- CIsig2 <- est.RHOb <- NULL



if(type == "biv"){

  bs <- NULL


  if(object$VC$Model != "BPO0" && is.null(object$VC$theta.fx)) bs <- rMVN(n.sim, mean = object$coefficients, sigma=object$Vb)  
  if(object$VC$Model == "BPO0") epds <- rep(0, 10)
  if(object$VC$Model == "B" && !is.null(object$VC$theta.fx)) epds <- rep(object$VC$theta.fx, 10)


if( is.null(object$VC$theta.fx) ){############


  if(object$VC$margins[2] %in% c(bin.link, cont1par) && object$VC$Model != "BPO0"){
  
  if(object$VC$Model == "BSS")  X3x <- object$X3s else X3x <- object$X3 
  
  if( !is.null(object$X3) ) epds <- X3x%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
  if(  is.null(object$X3) ) epds <- bs[,lf]
  
                                                                                  }
  
  
  
  
  
  if(object$VC$margins[2] %in% cont2par ){
  
  if( !is.null(object$X3) ) sigma2.st <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
  if(  is.null(object$X3) ) sigma2.st <- bs[, lf-1]
  

   sigma2 <- esp.tr(sigma2.st, object$VC$margins[2])$vrb  
   
   if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
   
   CIsig2 <- rowQuantiles(sigma2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CIsig2 <- t(CIsig2) 

  if( !is.null(object$X4) ) epds <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  if(  is.null(object$X4) ) epds <- bs[, lf]  
   
  } 
  
  
  
  
    if(object$VC$margins[2] %in% cont3par ){
    
    if( !is.null(object$X3) ) sigma2.st <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
    if(  is.null(object$X3) ) sigma2.st <- bs[, lf - 2]
    
    if( !is.null(object$X4) ) nu.st     <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)]) 
    if(  is.null(object$X4) ) nu.st     <- bs[, lf - 1]    
    
    sigma2 <- esp.tr(sigma2.st, object$VC$margins[2])$vrb  
 
     
     if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
     
     CIsig2 <- rowQuantiles(sigma2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     if( is.null(object$X3) ) CIsig2 <- t(CIsig2) 

     
     if(object$VC$margins[2] %in% c("DAGUM","SM")){
     
     nu <- esp.tr(nu.st, object$VC$margins[2])$vrb  

                                                  }
     

     
     if(  is.null(object$X4) ) nu <- t(as.matrix(nu))
     
     CInu <- rowQuantiles(nu, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     if( is.null(object$X3) ) CInu <- t(CInu) 

    if( !is.null(object$X5) ) epds <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
    if(  is.null(object$X5) ) epds <- bs[, lf]  
     
  }
  
 

}#############

}# biv









if(type == "copR"){



BivD <- object$VC$BivD 
if(object$surv == TRUE) BivD <- "N"


if(BivD == "T") ec <- 1 else ec <- 0  


###############################################

  if(  is.null(object$X3) ) epds <- bs[, lf]
  
  if( !is.null(object$X3) ){ 
  
   
   if(object$VC$margins[1] %in% c(object$VC$m2,object$VC$m3) && object$VC$margins[2] %in% c(object$VC$m2,object$VC$m3) && BivD == "T"){
   
     if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont2par) epds <- object$X6%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2)])
     if((object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par) || (object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par) ) epds <- object$X7%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2)])
     if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par) epds <- object$X8%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2+object$X8.d2)])  
     
   }else{
  
  
  if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont1par) epds <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])  
  if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont2par) epds <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont3par) epds <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])  	
  
  if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont1par && object$surv.flex == TRUE) epds <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont1par && object$surv.flex == TRUE) epds <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
  
  
  if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont2par) epds <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
  if((object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par) || (object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par) ) epds <- object$X6%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2)])
  if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par) epds <- object$X7%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2)])
  
  }
  

  	                   }

 
###############################################


if(object$VC$margins[1] %in% c(object$VC$m2,object$VC$m3) && object$VC$margins[2] %in% c(object$VC$m2,object$VC$m3) && BivD == "T"){###

  if(  is.null(object$X3) ) epds1 <- bs[, lf - ec]
  
  if( !is.null(object$X3) ){ 
  	
  if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont2par) epds1 <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
  if((object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par) || (object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par) ) epds1 <- object$X6%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2)])
  if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par) epds1 <- object$X7%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2)])
  
  	                   }

 est.dof <- dof.tr(epds1)$vao
 if( is.null(object$X3) ) est.dof <- t(as.matrix(est.dof))
   
 CIdof <- rowQuantiles(est.dof, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
 if( is.null(object$X3) ) CIdof <- t(CIdof) 
  

}###



###############################################

      if( is.null(object$X3)  ) {
      
  
if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont1par ){ ps1 <- lf - 1 - ec; ps2 <- lf - 1 - ec}      
if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont2par ){ ps1 <- lf - 1 - ec; ps2 <- lf - 1 - ec}
if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont3par ){ ps1 <- lf - 2 - ec; ps2 <- lf - 2 - ec}
if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont2par ){ ps1 <- lf - 2 - ec; ps2 <- lf - 1 - ec}

if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont1par && object$surv.flex == TRUE  ){ ps1 <- ps2 <- lf - 1 - ec}
if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont1par && object$surv.flex == TRUE  ){ ps1 <- ps2 <- lf - 2 - ec}


if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par ){ ps1 <- lf - 4 - ec; ps2 <- lf - 3 - ec}
if((object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par) || (object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par) ){ ps1 <- lf - 3 - ec; ps2 <- lf - 2 - ec}
      
                                sigma2.1.star <- t(as.matrix(bs[, ps1]))
                                sigma2.2.star <- t(as.matrix(bs[, ps2]))
                                
                                }
  
  
      if( !is.null(object$X3) ) {

if(!(object$VC$margins[1] %in% cont1par)){
       sigma2.1.star <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
       sigma2.2.star <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)]) 
                                          }

if(object$VC$margins[1] %in% cont1par){
       sigma2.1.star <- c(0, 0)
       sigma2.2.star <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
                                      }
                                      
if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont1par){
       sigma2.1.star <- sigma2.2.star <- c(0, 0)
                                      } 
                                      
if(object$VC$margins[1] %in% c(cont2par,cont3par) && object$VC$margins[2] %in% cont1par && object$surv.flex == TRUE  ){
       sigma2.1.star <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
       sigma2.2.star <- c(0, 0) 
                                      }                                      
 
 

                                }  

    sigma21 <- esp.tr(sigma2.1.star, object$VC$margins[1])$vrb  
    sigma22 <- esp.tr(sigma2.2.star, object$VC$margins[2])$vrb  

   if(!(object$VC$margins[1] %in% cont1par)) CIsig21 <- rowQuantiles(sigma21, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) else CIsig21 <- c(0,0)
   if(!(object$VC$margins[2] %in% cont1par)) CIsig22 <- rowQuantiles(sigma22, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) else CIsig22 <- c(0,0)
   if( is.null(object$X3) ){ CIsig21 <- t(CIsig21); CIsig22 <- t(CIsig22) }
  
###############################################
  


if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont3par ){  
       
    
  if( is.null(object$X3) )  {    pn1 <- lf - 2 - ec
                                 pn2 <- lf - 1 - ec
                                 nu1.st <- t(as.matrix(bs[, pn1]))
                                 nu2.st <- t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu1.st <- object$X5%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2)]) 
       nu2.st <- object$X6%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2)])
  
                }   
   
nu1 <- esp.tr(nu1.st, object$VC$margins[1])$vrb   
nu2 <- esp.tr(nu2.st, object$VC$margins[2])$vrb   


   CInu1 <- rowQuantiles(nu1, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   CInu2 <- rowQuantiles(nu2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   
      if( is.null(object$X3) ){ CInu1 <- t(CInu1); CInu2 <- t(CInu2) }

   
} 


if(object$VC$margins[1] %in% cont2par && object$VC$margins[2] %in% cont3par ){  
  

  if( is.null(object$X3) )  {  pn2 <- lf - 1 - ec; nu2.st <- t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu2.st <- object$X5%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2)]) 
  
                }   
  
nu2 <- esp.tr(nu2.st, object$VC$margins[2])$vrb   
 
   
   CInu2 <- rowQuantiles(nu2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CInu2 <- t(CInu2) 

   
} 




if(object$VC$margins[1] %in% c(cont3par) && object$VC$margins[2] %in% cont1par && object$surv.flex == TRUE){  
  

  if( is.null(object$X3) )  {  pn1 <- lf - 1 - ec; nu1.st <- t(as.matrix(bs[, pn1]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu1.st <- object$X4%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2)]) 
  
                }   
  
nu1 <- esp.tr(nu1.st, object$VC$margins[1])$vrb   
 
   
   CInu1 <- rowQuantiles(nu1, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CInu1 <- t(CInu1) 

   
} 













if(object$VC$margins[1] %in% cont1par && object$VC$margins[2] %in% cont3par ){  
  

  if( is.null(object$X3) )  {  pn2 <- lf - 1 - ec; nu2.st <- t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu2.st <- object$X4%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2)]) 
  
                }   
  
nu2 <- esp.tr(nu2.st, object$VC$margins[2])$vrb   
 
   
   CInu2 <- rowQuantiles(nu2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CInu2 <- t(CInu2) 

   
} 








  
  
  
if(object$VC$margins[1] %in% cont3par && object$VC$margins[2] %in% cont2par ){  
    
  if( is.null(object$X3) )  {  pn1 <- lf - 1 - ec; nu1.st <- t(as.matrix(bs[, pn1]))  } 
  
  if( !is.null(object$X3) ) {  
 
       nu1.st <- object$X5%*%t(bs[,(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2)]) 
  
                }   
  
nu1 <- esp.tr(nu1.st, object$VC$margins[1])$vrb  

   CInu1 <- rowQuantiles(nu1, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CInu1 <- t(CInu1)

   
}  

###############################################




}




if(type == "gamls"){



if(object$VC$margins[1] %in% c(cont2par,cont3par)){

if(object$VC$margins[1] %in% cont3par) ps1 <- lf - 1 else ps1 <- lf

	if(  is.null(object$X2) ) sigma2.1.star <- t(as.matrix(bs[, ps1]))
	if( !is.null(object$X2) ) sigma2.1.star <- object$X2%*%t(bs[,(object$X1.d2+1):(object$X1.d2+object$X2.d2)]) 


        sigma21 <- esp.tr(sigma2.1.star, object$VC$margins[1])$vrb  

        CIsig21 <- rowQuantiles(sigma21, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
    
        if( is.null(object$X2) ) CIsig21 <- t(CIsig21)
}  
  

  
  
if(object$VC$margins[1] %in% cont3par){  
    
  if(  is.null(object$X3) ) nu1.st <- t(as.matrix(bs[, lf]))
  if( !is.null(object$X3) ) nu1.st <- object$X3%*%t(bs[,(object$X1.d2 + object$X2.d2 + 1):(object$X1.d2 + object$X2.d2 + object$X3.d2)]) 
  
  nu1 <- esp.tr(nu1.st, object$VC$margins[1])$vrb   
  
  CInu1 <- rowQuantiles(nu1, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

  if( is.null(object$X3) ) CInu1 <- t(CInu1)

}   








}






if(type == "copSS"){



  if(object$VC$margins[2] %in% cont1par){
    
  if( !is.null(object$X3) ) epds <- object$X3s%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
  if(  is.null(object$X3) ) epds <- bs[,lf]
  
  }  
  
  
  
  
  if(object$VC$margins[2] %in% cont2par ){
  
  if( !is.null(object$X3) ) sigma2.st <- object$X3s%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
  if(  is.null(object$X3) ) sigma2.st <- bs[, lf-1]
  

   sigma2 <- esp.tr(sigma2.st, object$VC$margins[2])$vrb  
   
   if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
   
   CIsig2 <- rowQuantiles(sigma2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
   if( is.null(object$X3) ) CIsig2 <- t(CIsig2) 

  if( !is.null(object$X4) ) epds <- object$X4s%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  if(  is.null(object$X4) ) epds <- bs[, lf]  
   
  } 
  
  
  
    if(object$VC$margins[2] %in% cont3par ){
    
    if( !is.null(object$X3) ) sigma2.st <- object$X3s%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
    if(  is.null(object$X3) ) sigma2.st <- bs[, lf - 2]
    
    if( !is.null(object$X4) ) nu.st     <- object$X4s%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)]) 
    if(  is.null(object$X4) ) nu.st     <- bs[, lf - 1]    
    
    sigma2 <- esp.tr(sigma2.st, object$VC$margins[2])$vrb  
 
     
     if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
     
     CIsig2 <- rowQuantiles(sigma2, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     if( is.null(object$X3) ) CIsig2 <- t(CIsig2) 

     
     if(object$VC$margins[2] %in% c("DAGUM","SM")){
     
     nu <- esp.tr(nu.st, object$VC$margins[2])$vrb  

     
     }
     

     
     if(  is.null(object$X4) ) nu <- t(as.matrix(nu))
     
     CInu <- rowQuantiles(nu, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     if( is.null(object$X3) ) CInu <- t(CInu) 

    if( !is.null(object$X5) ) epds <- object$X5s%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
    if(  is.null(object$X5) ) epds <- bs[, lf]  
     
  }
  
  
   


}



##################################################################
##################################################################


if(type != "gamls"){


if( !is.null(object$VC$theta.fx) && type == "biv"){
   est.RHOb <- epds 
   est.RHOb <- t(as.matrix(est.RHOb))
  }else{
         est.RHOb <- teta.tr(object$VC, epds)$teta
         
         est.RHOb <- ass.ms(object$VC$BivD, object$VC$nCa, est.RHOb)$theta

         if( is.null(object$X3) ) est.RHOb <- t(as.matrix(est.RHOb))
       }
       
 CIrs <- rowQuantiles(est.RHOb, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
 if( is.null(object$X3) ) CIrs <- t(CIrs) 


##################################################################

nCa   <- object$VC$nCa
BivDt <- object$VC$BivD

  if(object$BivD %in% object$BivD2){
  
  if(object$BivD %in% object$BivD2[1:4]){  BivDt <- "C0"; nCa <- 3} 
  if(object$BivD %in% object$BivD2[5:8]){  BivDt <- "J0"; nCa <- 6}
  if(object$BivD %in% object$BivD2[9:12]){ BivDt <- "G0"; nCa <- 4}
  
                                   }
  
ass.msR <- ass.ms(BivDt, nCa, est.RHOb)
tau     <- ass.msR$tau
if(!is.null(object$X3) && BivDt %in% c("AMH", "FGM")) tau <- matrix(tau, nrow(object$X3), nrow(bs))
CIkt <- rowQuantiles(tau, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
if( is.null(object$X3) ) CIkt <- t(CIkt) 



 if(object$BivD %in% object$BivD2){ 
 
   if(length(object$theta) > 1){
 
     if( length(object$teta2) != 0){ CIkt[object$teta.ind2, ] <- -CIkt[object$teta.ind2, ]; CIkt[object$teta.ind2, c(1,2)] <- CIkt[object$teta.ind2, c(2,1)] 
                                     CIrs[object$teta.ind2, ] <- -CIrs[object$teta.ind2, ]; CIrs[object$teta.ind2, c(1,2)] <- CIrs[object$teta.ind2, c(2,1)]} 
                               }else{
 
     if( length(object$teta2) != 0){ CIkt <- -CIkt; CIkt[, c(1,2)] <- CIkt[, c(2,1)]
                                     CIrs <- -CIrs; CIrs[, c(1,2)] <- CIrs[, c(2,1)]} 
                                    }
 }

##################################################################
##################################################################


}










list(CIrs = CIrs, CIkt = CIkt, CIsig21 = CIsig21, CIsig22 = CIsig22, CInu1 = CInu1, CInu2 = CInu2, CIdof = CIdof, 
     CInu = CInu, CIsig2 = CIsig2, bs = bs, est.RHOb = est.RHOb)




}

