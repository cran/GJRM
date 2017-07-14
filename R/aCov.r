aCov <- function(x){


# no correction for gamlss

if(x$VC$triv == TRUE){

if( !is.null(x$X4) ) {mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X4) )  mul2 <- mul3 <- mul4 <- 1 
            
            
if(x$Model == "TSS"){ X2m <- x$X2s; X3m <- x$X3s} else{X2m <- x$X2; X3m <- x$X3}  
                                               
scores <- cbind( c(x$fit$dl.de1)*x$X1, 
                 c(x$fit$dl.de2)*X2m, 
                 c(x$fit$dl.de3)*X3m,
                 c(x$fit$dl.dtheta12.st)*mul2,
                 c(x$fit$dl.dtheta13.st)*mul3,
                 c(x$fit$dl.dtheta23.st)*mul4  ) 

##############################                 
# needs to be extended to TESS 
##############################
                 
}






if(x$VC$triv == FALSE){


cont1par <- c(x$VC$m1d)   
cont2par <- c(x$VC$m2,x$VC$m2d)   
cont3par <- c(x$VC$m3,x$VC$m3d)  
bin.link <- x$VC$bl  
 



###################################################################################
###################################################################################

if(x$Cont == "NO"){ ###


if( x$margins[2] %in% bin.link && x$Model != "BPO0" && is.null(x$VC$theta.fx)){

if(is.null(x$X3) )  mul <- 1
if(!is.null(x$X3) ){ if(x$Model == "BSS") mul <- x$X3s else mul <- x$X3  }

if(x$Model == "BSS") X2m <- x$X2s else X2m <- x$X2

scores <- cbind( c(x$fit$dl.dbe1)*x$X1, c(x$fit$dl.dbe2)*X2m, c(x$fit$dl.drho)*mul )

                                                    }
                                                    
 

if( x$Model == "BPO0" || ( x$Model == "B" && !is.null(x$VC$theta.fx) )   ){

scores <- cbind( c(x$fit$dl.dbe1)*x$X1, c(x$fit$dl.dbe2)*x$X2 )

                       }
                       
                       
                       
                       


if( x$margins[2] %in% cont1par ){

if( !is.null(x$X3) ){ if(x$VC$ccss == "yes") mul <- x$X3s else mul <- x$X3} 
if(  is.null(x$X3) )  mul <- 1 
                      
                      
                      
if(x$VC$ccss == "yes") X2m <- x$X2s else X2m <- x$X2                       
                      
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*X2m, 
                 c(x$fit$dl.dteta.st)*mul       )                                           
                                 }




if( x$margins[2] %in% cont2par ){

if( !is.null(x$X3) && !is.null(x$X4) ) { if(x$VC$ccss == "yes"){ mul1 <- x$X3s; mul2 <- x$X4s} else {mul1 <- x$X3; mul2 <- x$X4}         } 
if(  is.null(x$X3) &&  is.null(x$X4) )  mul1 <- mul2 <- 1 
         
         
if(x$VC$ccss == "yes") X2m <- x$X2s else X2m <- x$X2                       
         
         
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*X2m, 
                 c(x$fit$dl.dsigma.st)*mul1,
                 c(x$fit$dl.dteta.st)*mul2       )                                           
                                 }

if( x$margins[2] %in% cont3par ){

if( !is.null(x$X3) && !is.null(x$X4) && !is.null(x$X5)) {if(x$VC$ccss == "yes"){mul1 <- x$X3s; mul2 <- x$X4s; mul3 <- x$X5s}else{mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5} } 
if(  is.null(x$X3) &&  is.null(x$X4) &&  is.null(x$X5))  mul1 <- mul2 <- mul3 <- 1 
      

if(x$VC$ccss == "yes") X2m <- x$X2s else X2m <- x$X2                       
      
      
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*X2m, 
                 c(x$fit$dl.dsigma.st)*mul1,
                 c(x$fit$dl.dnu.st)*mul2,
                 c(x$fit$dl.dteta.st)*mul3       )     
                                       
                                }




} ###

###################################################################################
###################################################################################



###################################################################################
###################################################################################

if(x$Cont == "YES"){ ###


if( x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) && x$BivD == "T"){



if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu.st)*mul3,                 
                 c(x$fit$dl.dteta.st)*mul4       ) 
                 
                 
}





if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6; mul5 <- x$X7; mul6 <- x$X8} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- mul5 <- mul6 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu1.st)*mul3,
                 c(x$fit$dl.dnu2.st)*mul4,
                 c(x$fit$dl.dnu.st)*mul5,
                 c(x$fit$dl.dteta.st)*mul6       ) 
                 
                 
}





if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6; mul5 <- x$X7} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- mul5 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu2.st)*mul3,
                 c(x$fit$dl.dnu.st)*mul4,
                 c(x$fit$dl.dteta.st)*mul5       ) 
                 
                 
}




if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6; mul5 <- x$X7} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- mul5 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu1.st)*mul3,
                 c(x$fit$dl.dnu.st)*mul4,
                 c(x$fit$dl.dteta.st)*mul5       ) 
                 
                 
                                                       }



}else{





if( x$margins[1] %in% cont1par && x$margins[2] %in% cont1par ){

if( !is.null(x$X3) ) mul1 <- x$X3 
if(  is.null(x$X3) ) mul1 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dteta.st)*mul1 ) 
                                 
}



if( x$margins[1] %in% cont1par && x$margins[2] %in% cont2par ){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4} 
if(  is.null(x$X3) )  mul1 <- mul2 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma22.st)*mul1,
                 c(x$fit$dl.dteta.st)*mul2       ) 
                                 
}


if( x$margins[1] %in% cont1par && x$margins[2] %in% cont3par ){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma22.st)*mul1,
                 c(x$fit$dl.dnu2.st)*mul2,
                 c(x$fit$dl.dteta.st)*mul3       ) 
                 
                 
}






if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dteta.st)*mul3       ) 
                 
                 
}




if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6; mul5 <- x$X7} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- mul5 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu1.st)*mul3,
                 c(x$fit$dl.dnu2.st)*mul4,
                 c(x$fit$dl.dteta.st)*mul5       ) 
                 
                 
}



if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu2.st)*mul3,
                 c(x$fit$dl.dteta.st)*mul4       ) 
                 
                 
}



if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu1.st)*mul3,
                 c(x$fit$dl.dteta.st)*mul4       ) 
                 
                 
}


}



}


###################################################################################
###################################################################################



} # end TRIV




scores 



}


