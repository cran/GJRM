pen <- function(qu.mag, sp, VC, univ, l.splist){

##############################################################
    ma1 <- matrix(0,VC$gp1,VC$gp1)   
    if(l.splist$l.sp1 == 0) EQ1P <- adiag(ma1)
    
    if(l.splist$l.sp1 != 0){
    ind <- 1:l.splist$l.sp1
    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
    S1 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    
    if( length(unique(offtemp)) != length(offtemp) ) S1 <- SS(offtemp, S1)
    
    
    S1 <- do.call(adiag, lapply(S1, unlist))	
    
    EQ1P <- adiag(ma1, S1)
     
                           } 
                   
##############################################################
    
##############################################################   

    if( is.null(VC$gp2.inf) ) ma2 <- matrix(0,VC$gp2,VC$gp2) else ma2 <- matrix(0,VC$gp2.inf,VC$gp2.inf) 
    
    if(l.splist$l.sp2 == 0) EQ2P <- adiag(ma2)    
    
    if(l.splist$l.sp2 != 0){
    ind <- (l.splist$l.sp1 + 1):(l.splist$l.sp1 + l.splist$l.sp2)
    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
    S2 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    
    if( length(unique(offtemp)) != length(offtemp) ) S2 <- SS(offtemp, S2)
   
    S2 <- do.call(adiag, lapply(S2, unlist))
    EQ2P <- adiag(ma2, S2)
                   }
                            
##############################################################
    
if(!is.null(VC$gp3)){ # this starts after first two equations

 EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- NULL 

    ma3 <- matrix(0,VC$gp3,VC$gp3) 
    if(l.splist$l.sp3 == 0) EQ3P <- adiag(ma3)    
    
    if(l.splist$l.sp3 != 0){
    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3)
    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
    S3 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    
    if( length(unique(offtemp)) != length(offtemp) ) S3 <- SS(offtemp, S3)

    S3 <- do.call(adiag, lapply(S3, unlist))
    EQ3P <- adiag(ma3, S3)
                   }
        
    if(!is.null(VC$gp4)){
    
    ma4 <- matrix(0,VC$gp4,VC$gp4) 
    if(l.splist$l.sp4 == 0) EQ4P <- adiag(ma4)    
    
    if(l.splist$l.sp4 != 0){
    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4)
    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
    S4 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    
    if( length(unique(offtemp)) != length(offtemp) ) S4 <- SS(offtemp, S4)    
    
    S4 <- do.call(adiag, lapply(S4, unlist))
    EQ4P <- adiag(ma4, S4)
                   }
                        }
                        
	    if(!is.null(VC$gp5)){
	    
	    ma5 <- matrix(0,VC$gp5,VC$gp5) 
	    if(l.splist$l.sp5 == 0) EQ5P <- adiag(ma5)    
	    
	    if(l.splist$l.sp5 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5)
	    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
	    S5 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    
            if( length(unique(offtemp)) != length(offtemp) ) S5 <- SS(offtemp, S5)
   	    
	    S5 <- do.call(adiag, lapply(S5, unlist))
	    EQ5P <- adiag(ma5, S5)
	                     } 
                                }
                                
            if(!is.null(VC$gp6)){
	    
	    ma6 <- matrix(0,VC$gp6,VC$gp6) 
	    if(l.splist$l.sp6 == 0) EQ6P <- adiag(ma6)    
	    
	    if(l.splist$l.sp6 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6)
	    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
	    S6 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    
            if( length(unique(offtemp)) != length(offtemp) ) S6 <- SS(offtemp, S6)   
	    
	    S6 <- do.call(adiag, lapply(S6, unlist))
	    EQ6P <- adiag(ma6, S6)
	                     } 
                                }  
                                
            if(!is.null(VC$gp7)){
	    
	    ma7 <- matrix(0,VC$gp7,VC$gp7) 
	    if(l.splist$l.sp7 == 0) EQ7P <- adiag(ma7)    
	    
	    if(l.splist$l.sp7 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + l.splist$l.sp7)
	    offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))
	    S7 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    
            if( length(unique(offtemp)) != length(offtemp) ) S7 <- SS(offtemp, S7)    
	    
	    S7 <- do.call(adiag, lapply(S7, unlist))
	    EQ7P <- adiag(ma7, S7)
	                     }               
                                }  
                                
            if(!is.null(VC$gp8)){
	    
	    ma8 <- matrix(0,VC$gp8,VC$gp8) 
	    if(l.splist$l.sp8 == 0) EQ8P <- adiag(ma8)    
	    
	    if(l.splist$l.sp8 != 0){
	    ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + l.splist$l.sp7 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + l.splist$l.sp7 + l.splist$l.sp8)
            offtemp <- as.numeric(as.factor( qu.mag$off[ind] ))    
	    S8 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)

            if( length(unique(offtemp)) != length(offtemp) ) S8 <- SS(offtemp, S8)    
	    
	    S8 <- do.call(adiag, lapply(S8, unlist))
	    EQ8P <- adiag(ma8, S8)
	                     }               
                                }                                  
                                
                                
}else{



if(VC$univ.gamls == FALSE){
            
   
   if(VC$margins[1] %in% c(VC$m2,VC$m3) && VC$margins[2] %in% c(VC$m2,VC$m3) && VC$BivD == "T"){
   
       if(VC$margins[1] %in% VC$m2 && VC$margins[2] %in% VC$m2 ) {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0; EQ7P <- EQ8P <- NULL     } 
       if(VC$margins[1] %in% VC$m3 && VC$margins[2] %in% VC$m3 ) {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0; EQ7P <- EQ8P <- 0        } 
       if(VC$margins[1] %in% VC$m2 && VC$margins[2] %in% VC$m3 ) {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0; EQ7P <- 0; EQ8P <- NULL  }  
       if(VC$margins[1] %in% VC$m3 && VC$margins[2] %in% VC$m2 ) {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0; EQ7P <- 0; EQ8P <- NULL  }  
       
   }else{
   
    if(VC$margins[1] %in% c(VC$bl)        && VC$Model != "BPO0")                 {EQ3P <- 0; EQ4P <- NULL; EQ5P <- NULL; EQ6P <- EQ7P <- EQ8P <- NULL  }
    
    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% c(VC$bl,VC$m1d)) {EQ3P <- 0; EQ4P <- NULL; EQ5P <- NULL; EQ6P <- EQ7P <- EQ8P <- NULL  }
    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% VC$m2d)          {EQ3P <- 0; EQ4P <- 0;    EQ5P <- NULL; EQ6P <- EQ7P <- EQ8P <- NULL  }
    

    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% VC$m2)           {EQ3P <- 0; EQ4P <- 0;    EQ5P <- NULL; EQ6P <- EQ7P <- EQ8P <- NULL  }
    if(VC$margins[1] %in% c(VC$bl,VC$m1d) && VC$margins[2] %in% VC$m3)           {EQ3P <- 0; EQ4P <- 0;    EQ5P <- 0;    EQ6P <- EQ7P <- EQ8P <- NULL  }     
    
    if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% c(VC$m2,VC$m2d)  )  {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- NULL; EQ7P <- EQ8P <- NULL  } 
    
    if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% c(VC$bl)  )  {EQ3P <- 0; EQ4P <- 0; EQ5P <- NULL; EQ6P <- NULL; EQ7P <- EQ8P <- NULL  }     
    if(VC$margins[1] %in% c(VC$m3) && VC$margins[2] %in% c(VC$bl)  )         {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- NULL; EQ7P <- EQ8P <- NULL  }     

    
    
    if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% VC$m3            )  {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0;    EQ7P <- 0; EQ8P <- NULL  }   
    if(VC$margins[1] %in% c(VC$m2,VC$m2d) && VC$margins[2] %in% VC$m3            )  {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0;    EQ7P <- EQ8P <- NULL  }  
    if(VC$margins[1] %in% VC$m3           && VC$margins[2] %in% c(VC$m2,VC$m2d)  )  {EQ3P <- 0; EQ4P <- 0; EQ5P <- 0; EQ6P <- 0;    EQ7P <- EQ8P <- NULL  }  
    
    if(VC$Model == "B" && !is.null(VC$theta.fx)) {EQ3P <- EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- NULL                 }
    if(VC$Model == "BPO0")                                                       {EQ3P <- EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- NULL                 }
    # BPO0 put here below to avoid conflict with bl bl second line above
    }


    # maybe not that efficient, it could be done above but done here for the moment


}


}  
    
    
    
    
    
    
if(VC$triv == FALSE){    ### TRIV
    
    
    if(univ == 0) S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P, EQ7P, EQ8P)
    
    
    if(univ == 2){
    
       if(VC$margins[1] %in% c(VC$m1d,VC$bl))        S.h <- adiag(EQ1P)   
       
       if(VC$margins[1] %in% c(VC$bl) && !is.null(VC$gp2.inf)) S.h <- adiag(EQ1P, EQ2P)   # this is for informative censoring
       
       if(VC$margins[1] %in% c(VC$m2,VC$m2d) ) S.h <- adiag(EQ1P, EQ2P)         
       if(VC$margins[1] %in% VC$m3           ) S.h <- adiag(EQ1P, EQ2P, EQ3P)   
   
    }
                
        
}    ### TRIV    
        
        
        
        
        
        
####################################################################################################        
        
if(VC$triv == TRUE){ # here I am setting up S.h but also changing/expanding qu.mag 

S.h <- adiag(EQ1P, EQ2P, EQ3P)   

if(VC$penCor %in% c("unpen") && VC$l.flist == 3) S.h <- adiag(S.h, matrix(0, 3, 3))        
 
if(VC$penCor %in% c("unpen") && VC$l.flist == 6) S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P)   
        
 
if(VC$penCor %in% c("ridge")){

  A <- diag(c(1, 1, 1))
  if( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0) qu.mag$Ss[[1]]                   <- A
  if( VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0) qu.mag$Ss[[length(qu.mag$Ss)+1]] <- A

  S.h <- adiag(S.h, sp[length(sp)]*qu.mag$Ss[[length(qu.mag$Ss)]]) 
  
                             }
  
                                          
        
} 

####################################################################################################        

        
list(S.h = S.h, qu.mag = qu.mag)
    
   
 
         }


