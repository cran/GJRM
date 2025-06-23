gamls.upsv <- function(gamlss1, gamlss2, margins, M, l.flist, nstv, VC, GAM, SP, type = "copR"){
  
sp1 <- SP$sp1
sp2 <- SP$sp2
sp3 <- SP$sp3
sp4 <- SP$sp4
sp5 <- SP$sp5
sp6 <- SP$sp6
sp7 <- SP$sp7
sp8 <- SP$sp8
sp9 <- SP$sp9

# this function should be cleaner as there are many useless reps but ok for now

if(type == "ROY"){

   
      if(margins[2] %in% c(M$m1d,M$bl) && margins[3] %in% c(M$m1d,M$bl)){
    
        b2 <- gamlss1$coefficients[1:VC$X2.d2] # this is actually gamlss2, done to avoid changing arguments of the function
        b3 <- gamlss2$coefficients[1:VC$X3.d2] # this is actually gamlss3
 
        start.v <- c(b2, b3)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp
                                                              }
                                                                  
                                                           
      if(margins[2] %in% c(M$m2d,M$m2) && margins[3] %in% c(M$m2d,M$m2)){
    
        b2.1 <- gamlss1$coefficients[1:VC$X2.d2] 
        b2.2 <- gamlss1$coefficients[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X4.d2)]
        
        b3.1 <- gamlss2$coefficients[1:VC$X3.d2]
        b3.2 <- gamlss2$coefficients[(VC$X3.d2 + 1):(VC$X3.d2 + VC$X5.d2)]
 
        start.v <- c(b2.1, b3.1, b2.2, b3.2)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp[1:VC$l.sp2]
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[1:VC$l.sp3]
        if( VC$l.sp4 != 0 ) sp4 <- gamlss1$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]
        if( VC$l.sp5 != 0 ) sp5 <- gamlss2$sp[(VC$l.sp3 + 1):(VC$l.sp3 + VC$l.sp5)]        
    
                                                                         }   
                                                                         
                                                                         
      if(margins[2] %in% c(M$m3) && margins[3] %in% c(M$m3)){
    
        b2.1 <- gamlss1$coefficients[1:VC$X2.d2] 
        b2.2 <- gamlss1$coefficients[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X4.d2)]
        b2.3 <- gamlss1$coefficients[(VC$X2.d2 + VC$X4.d2 + 1):(VC$X2.d2 + VC$X4.d2 + VC$X6.d2)]
        
        
        b3.1 <- gamlss2$coefficients[1:VC$X3.d2]
        b3.2 <- gamlss2$coefficients[(VC$X3.d2 + 1):(VC$X3.d2 + VC$X5.d2)]
        b3.3 <- gamlss2$coefficients[(VC$X3.d2 + VC$X5.d2 + 1):(VC$X3.d2 + VC$X5.d2 + VC$X7.d2)]
        
 
        start.v <- c(b2.1, b3.1, b2.2, b3.2, b2.3, b3.3)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp[1:VC$l.sp2]
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[1:VC$l.sp3]
        if( VC$l.sp4 != 0 ) sp4 <- gamlss1$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]
        if( VC$l.sp5 != 0 ) sp5 <- gamlss2$sp[(VC$l.sp3 + 1):(VC$l.sp3 + VC$l.sp5)]     
        if( VC$l.sp6 != 0 ) sp6 <- gamlss1$sp[(VC$l.sp2 + VC$l.sp4 + 1):(VC$l.sp2 + VC$l.sp4 + VC$l.sp6)]        
        if( VC$l.sp7 != 0 ) sp7 <- gamlss2$sp[(VC$l.sp3 + VC$l.sp5 + 1):(VC$l.sp3 + VC$l.sp5 + VC$l.sp7)]        

                                                             }  
                                                             
                                                             
      if(margins[2] %in% c(M$m2) && margins[3] %in% c(M$m3)){
    
        b2.1 <- gamlss1$coefficients[1:VC$X2.d2] 
        b2.2 <- gamlss1$coefficients[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X4.d2)]
        
        
        b3.1 <- gamlss2$coefficients[1:VC$X3.d2]
        b3.2 <- gamlss2$coefficients[(VC$X3.d2 + 1):(VC$X3.d2 + VC$X5.d2)]
        b3.3 <- gamlss2$coefficients[(VC$X3.d2 + VC$X5.d2 + 1):(VC$X3.d2 + VC$X5.d2 + VC$X6.d2)]
        
 
        start.v <- c(b2.1, b3.1, b2.2, b3.2, b3.3)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp[1:VC$l.sp2]
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[1:VC$l.sp3]
        if( VC$l.sp4 != 0 ) sp4 <- gamlss1$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]
        if( VC$l.sp5 != 0 ) sp5 <- gamlss2$sp[(VC$l.sp3 + 1):(VC$l.sp3 + VC$l.sp5)]     
        if( VC$l.sp6 != 0 ) sp6 <- gamlss2$sp[(VC$l.sp3 + VC$l.sp5 + 1):(VC$l.sp3 + VC$l.sp5 + VC$l.sp6)]        

                                                             }                                                               
                                                             

      if(margins[2] %in% c(M$m3) && margins[3] %in% c(M$m2)){
    
        b2.1 <- gamlss1$coefficients[1:VC$X2.d2] 
        b2.2 <- gamlss1$coefficients[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X4.d2)]
        b2.3 <- gamlss1$coefficients[(VC$X2.d2 + VC$X4.d2 + 1):(VC$X2.d2 + VC$X4.d2 + VC$X6.d2)]
        
        
        b3.1 <- gamlss2$coefficients[1:VC$X3.d2]
        b3.2 <- gamlss2$coefficients[(VC$X3.d2 + 1):(VC$X3.d2 + VC$X5.d2)]
        
 
        start.v <- c(b2.1, b3.1, b2.2, b3.2, b2.3)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp[1:VC$l.sp2]
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[1:VC$l.sp3]
        if( VC$l.sp4 != 0 ) sp4 <- gamlss1$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]
        if( VC$l.sp5 != 0 ) sp5 <- gamlss2$sp[(VC$l.sp3 + 1):(VC$l.sp3 + VC$l.sp5)]     
        if( VC$l.sp6 != 0 ) sp6 <- gamlss1$sp[(VC$l.sp2 + VC$l.sp4 + 1):(VC$l.sp2 + VC$l.sp4 + VC$l.sp6)]        

                                                             }                                                              
                                                                         

      if(margins[2] %in% c(M$m1d) && margins[3] %in% c(M$m2d)){
    
        b2.1 <- gamlss1$coefficients 
        
        b3.1 <- gamlss2$coefficients[1:VC$X3.d2]
        b3.2 <- gamlss2$coefficients[(VC$X3.d2 + 1):(VC$X3.d2 + VC$X4.d2)]
 
        start.v <- c(b2.1, b3.1, b3.2)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[1:VC$l.sp3]
        if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp3 + 1):(VC$l.sp3 + VC$l.sp4)]
    
                                                               } 


      if(margins[2] %in% c(M$m2d) && margins[3] %in% c(M$m1d)){
    
        b2.1 <- gamlss1$coefficients[1:VC$X2.d2] 
        b2.2 <- gamlss1$coefficients[(VC$X2.d2 + 1):(VC$X2.d2 + VC$X4.d2)]
        
        b3.1 <- gamlss2$coefficients
 
        start.v <- c(b2.1, b3.1, b2.2)
  
        if( VC$l.sp2 != 0 ) sp2 <- gamlss1$sp[1:VC$l.sp2]
        if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp
        if( VC$l.sp4 != 0 ) sp4 <- gamlss1$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]
    
                                                                         } 


}








if(type == "biv"){
  


  
  if( l.flist == 2 ) start.v <- c(GAM$gam1$coefficients, gamlss2$coefficients, VC$i.rho)              # PO, N, DAGUM etc
  if( l.flist == 3 ) start.v <- c(GAM$gam1$coefficients, gamlss2$coefficients, GAM$gam3$coefficients) # PO only
  if( l.flist == 4 ) start.v <- c(GAM$gam1$coefficients, gamlss2$coefficients, GAM$gam4$coefficients) # N
  if( l.flist == 5 ) start.v <- c(GAM$gam1$coefficients, gamlss2$coefficients, GAM$gam5$coefficients) # DAGUM
  if( l.flist == 6 ) start.v <- c(GAM$gam1$coefficients, gamlss2$coefficients, GAM$gam6$coefficients) # DAGUM
  

if( margins[2] %in% c(M$m1d) ) if(VC$l.sp2 != 0) sp2 <- gamlss2$sp[1:VC$l.sp2] 
                                

if( margins[2] %in% c(M$m2,M$m2d,M$m3) ){
  	if(VC$l.sp2 != 0) sp2 <- gamlss2$sp[1:VC$l.sp2] 
  	if(VC$l.sp3 != 0) sp3 <- gamlss2$sp[VC$l.sp2 + (1:VC$l.sp3)] 
                                    }
  
if( margins[2] %in% M$m3 ){ 
  	if(VC$l.sp4 != 0) sp4 <- gamlss2$sp[VC$l.sp2 + VC$l.sp3 + (1:VC$l.sp4)] 
  	                  } 


}







if(type == "copR"){

  if(margins[1] %in% c(M$m2) && margins[2] %in% c(M$bl) && l.flist == 2 && M$surv == TRUE){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[VC$X1.d2+1]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  start.v <- c(b1, b2, s1, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }
  
  
  if(margins[1] %in% c(M$bl) && margins[2] %in% c(M$bl) && M$surv == TRUE){#  && M$end.surv == TRUE){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  if(l.flist == 2) start.v <- c(b1, b2, VC$i.rho)
  if(l.flist  > 2) start.v <- c(b1, b2, GAM$gam3$coefficients)
  
  names(start.v) <- nstv
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }  
  
  
  
  
  
  
  
  if(margins[1] %in% c(M$m3) && margins[2] %in% c(M$bl) && l.flist == 2 && M$surv == TRUE){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[VC$X1.d2+1]
  n1 <- gamlss1$coefficients[VC$X1.d2+2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  start.v <- c(b1, b2, s1, n1, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }  
  
  
  
  
  
  if(margins[1] %in% c(M$m2) && margins[2] %in% c(M$bl) && l.flist > 2 && M$surv == TRUE){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  start.v <- c(b1, b2, s1, GAM$gam4$coefficients); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  
  }
  
  
  
  
  
  
  if(margins[1] %in% c(M$m3) && margins[2] %in% c(M$bl) && l.flist > 2 && M$surv == TRUE){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
  n1 <- gamlss1$coefficients[(VC$X1.d2+VC$X3.d2+1):(VC$X1.d2+VC$X3.d2+VC$X5.d2)]  
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  start.v <- c(b1, b2, s1, n1, GAM$gam5$coefficients); names(start.v) <- nstv  
  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  if( VC$l.sp5 != 0 ) sp5 <- gamlss1$sp[(VC$l.sp1 + VC$l.sp3 + 1):(VC$l.sp1 + VC$l.sp3 + VC$l.sp5)]
    
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
 
  }   
  
  
  
  




  if(margins[1] %in% c(M$m1d,M$bl) && margins[2] %in% c(M$m1d,M$bl) && l.flist == 2){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  start.v <- c(b1, b2, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  } 
  
  
  if(margins[1] %in% c(M$m1d,M$bl) && margins[2] %in% c(M$m1d,M$bl) && l.flist > 2){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
 
  start.v <- c(b1, b2, GAM$gam3$coefficients); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }   
  
  
  if(margins[1] %in% M$m1d && margins[2] %in% c(M$m2,M$m2d) && l.flist == 2){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[VC$X2.d2+1]  
 
  start.v  <- c(b1, b2, s2, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  } 
  
  if(margins[1] %in% M$m1d && margins[2] %in% c(M$m2,M$m2d) && l.flist > 2){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)]  
 
  start.v  <- c(b1, b2, s2, GAM$gam4$coefficients); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp3)]  

  }  
  
  #
  # 
  
  
  
  
  if(margins[1] %in% c(M$m1d) && margins[2] %in% M$m3 && l.flist == 2){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[VC$X2.d2+1] 
    n2 <- gamlss2$coefficients[VC$X2.d2+2]    
   
    start.v  <- c(b1, b2, s2, n2, VC$i.rho); names(start.v) <- nstv  
    
    if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
    if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
    
    }  
    
    if(margins[1] %in% c(M$m1d) && margins[2] %in% M$m3 && l.flist > 2){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)]
    n2 <- gamlss2$coefficients[(VC$X2.d2+VC$X3.d2+1):(VC$X2.d2+VC$X3.d2+VC$X4.d2)]  
   
    start.v  <- c(b1, b2, s2, n2, GAM$gam5$coefficients); names(start.v) <- nstv  
    
    
    if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
      
    if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
    if( VC$l.sp3 != 0 ) sp3 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp3)]  
    if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + VC$l.sp3 + 1):(VC$l.sp2 + VC$l.sp3 + VC$l.sp4)]  
    
    }     
    
  #
  #  
  
  
  
  

  
  if((margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d) && l.flist == 2 && M$BivD != "T") || (margins[1] %in% c(M$m2d) && margins[2] %in% c(M$m2d) && l.flist == 2 && M$BivD == "T")){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[VC$X1.d2+1]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[VC$X2.d2+1]  
 
  start.v  <- c(b1, b2, s1, s2, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }
  
  
  
  if((margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d) && l.flist > 2 && M$BivD != "T") || (margins[1] %in% c(M$m2d) && margins[2] %in% c(M$m2d) && l.flist > 2 && M$BivD == "T")){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, GAM$gam5$coefficients); names(start.v) <- nstv # dropped gam6 here   
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  

  } 
  
  
   
  

    
  if(margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2) && l.flist == 2 && M$BivD == "T"){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[VC$X1.d2+1]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[VC$X2.d2+1]  
   
    start.v  <- c(b1, b2, s1, s2, VC$dof.st, VC$i.rho); names(start.v) <- nstv  
    
    if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
    if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
    
    }    
    
    
    

  
  
     if(margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2) && l.flist > 2 && M$BivD == "T"){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]  
   
    start.v  <- c(b1, b2, s1, s2, GAM$gam5$coefficients, GAM$gam6$coefficients); names(start.v) <- nstv  
    
    if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
    if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
    
    if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
    if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
  
  }   
  
  
  
  
  
  #
  #
  
  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && l.flist == 2 && M$BivD != "T"){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[VC$X1.d2+1]
  n1 <- gamlss1$coefficients[VC$X1.d2+2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[VC$X2.d2+1] 
  n2 <- gamlss2$coefficients[VC$X2.d2+2]    
 
  start.v  <- c(b1, b2, s1, s2, n1, n2, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp  
  
  
  }  



  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && l.flist == 2 && M$BivD == "T"){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[VC$X1.d2+1]
  n1 <- gamlss1$coefficients[VC$X1.d2+2]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[VC$X2.d2+1] 
  n2 <- gamlss2$coefficients[VC$X2.d2+2]    
 
  start.v  <- c(b1, b2, s1, s2, n1, n2, VC$dof.st, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp  
  
  
  } 


  
  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && l.flist > 2 && M$BivD != "T"){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
  n1 <- gamlss1$coefficients[(VC$X1.d2+VC$X3.d2+1):(VC$X1.d2+VC$X3.d2+VC$X5.d2)]  
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]
  n2 <- gamlss2$coefficients[(VC$X2.d2+VC$X4.d2+1):(VC$X2.d2+VC$X4.d2+VC$X6.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, n1, n2, GAM$gam7$coefficients); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  if( VC$l.sp5 != 0 ) sp5 <- gamlss1$sp[(VC$l.sp1 + VC$l.sp3 + 1):(VC$l.sp1 + VC$l.sp3 + VC$l.sp5)]
  
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
  if( VC$l.sp6 != 0 ) sp6 <- gamlss2$sp[(VC$l.sp2 + VC$l.sp4 + 1):(VC$l.sp2 + VC$l.sp4 + VC$l.sp6)]  
  
  } 



  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && l.flist > 2 && M$BivD == "T"){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
  n1 <- gamlss1$coefficients[(VC$X1.d2+VC$X3.d2+1):(VC$X1.d2+VC$X3.d2+VC$X5.d2)]  
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]
  n2 <- gamlss2$coefficients[(VC$X2.d2+VC$X4.d2+1):(VC$X2.d2+VC$X4.d2+VC$X6.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, n1, n2, GAM$gam7$coefficients, GAM$gam8$coefficients); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  if( VC$l.sp5 != 0 ) sp5 <- gamlss1$sp[(VC$l.sp1 + VC$l.sp3 + 1):(VC$l.sp1 + VC$l.sp3 + VC$l.sp5)]
  
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
  if( VC$l.sp6 != 0 ) sp6 <- gamlss2$sp[(VC$l.sp2 + VC$l.sp4 + 1):(VC$l.sp2 + VC$l.sp4 + VC$l.sp6)]  
  
  } 

    
  
  #
  #
  
  if((margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% M$m3 && l.flist == 2 && M$BivD != "T") || (margins[1] %in% c(M$m2d) && margins[2] %in% M$m3 && l.flist == 2 && M$BivD == "T")){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[VC$X1.d2+1]
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[VC$X2.d2+1] 
  n2 <- gamlss2$coefficients[VC$X2.d2+2]    
 
  start.v  <- c(b1, b2, s1, s2, n2, VC$i.rho); names(start.v) <- nstv  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }  
  
  
  
  
  if((margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% M$m3 && l.flist > 2 && M$BivD != "T") || (margins[1] %in% c(M$m2d) && margins[2] %in% M$m3 && l.flist > 2 && M$BivD == "T")){
  
  b1 <- gamlss1$coefficients[1:VC$X1.d2]
  s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)] 
  b2 <- gamlss2$coefficients[1:VC$X2.d2]
  s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]
  n2 <- gamlss2$coefficients[(VC$X2.d2+VC$X4.d2+1):(VC$X2.d2+VC$X4.d2+VC$X5.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, n2, GAM$gam6$coefficients); names(start.v) <- nstv  
  
  
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
    
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
  if( VC$l.sp5 != 0 ) sp5 <- gamlss2$sp[(VC$l.sp2 + VC$l.sp4 + 1):(VC$l.sp2 + VC$l.sp4 + VC$l.sp5)]  
  
  } 
  
  
  
    
    
  if(margins[1] %in% c(M$m2) && margins[2] %in% M$m3 && l.flist == 2 && M$BivD == "T"){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[VC$X1.d2+1]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[VC$X2.d2+1] 
    n2 <- gamlss2$coefficients[VC$X2.d2+2]    
   
    start.v  <- c(b1, b2, s1, s2, n2, VC$dof.st, VC$i.rho); names(start.v) <- nstv  
    
    if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
    if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp
    
    }     
    
    
    

  

    if(margins[1] %in% c(M$m2) && margins[2] %in% M$m3 && l.flist > 2 && M$BivD == "T"){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)] 
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]
    n2 <- gamlss2$coefficients[(VC$X2.d2+VC$X4.d2+1):(VC$X2.d2+VC$X4.d2+VC$X5.d2)]  
   
    start.v  <- c(b1, b2, s1, s2, n2, GAM$gam6$coefficients, GAM$gam7$coefficients); names(start.v) <- nstv  
    
    
    if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
    if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
      
    if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
    if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
    if( VC$l.sp5 != 0 ) sp5 <- gamlss2$sp[(VC$l.sp2 + VC$l.sp4 + 1):(VC$l.sp2 + VC$l.sp4 + VC$l.sp5)]  
    
  }  
  
  
  
  
  #
  #  
  
  if((margins[1] %in% M$m3 && margins[2] %in% c(M$m2,M$m2d) && l.flist == 2 && M$BivD != "T") || (margins[1] %in% M$m3 && margins[2] %in% c(M$m2d) && l.flist == 2 && M$BivD == "T")){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[VC$X1.d2+1]
    n1 <- gamlss1$coefficients[VC$X1.d2+2]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[VC$X2.d2+1] 
   
    start.v  <- c(b1, b2, s1, s2, n1, VC$i.rho); names(start.v) <- nstv  
    
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp    
    
    }  
    
    
    
    
  if((margins[1] %in% M$m3 && margins[2] %in% c(M$m2,M$m2d) && l.flist > 2 && M$BivD != "T") || (margins[1] %in% M$m3 && margins[2] %in% c(M$m2d) && l.flist > 2 && M$BivD == "T")){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
    n1 <- gamlss1$coefficients[(VC$X1.d2+VC$X3.d2+1):(VC$X1.d2+VC$X3.d2+VC$X5.d2)]  
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]
   
    start.v  <- c(b1, b2, s1, s2, n1, GAM$gam6$coefficients); names(start.v) <- nstv  
   
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  if( VC$l.sp5 != 0 ) sp5 <- gamlss1$sp[(VC$l.sp1 + VC$l.sp3 + 1):(VC$l.sp1 + VC$l.sp3 + VC$l.sp5)]
    
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
   
    }  
    
    
   
    
  if(margins[1] %in% M$m3 && margins[2] %in% c(M$m2) && l.flist == 2 && M$BivD == "T"){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[VC$X1.d2+1]
    n1 <- gamlss1$coefficients[VC$X1.d2+2]
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[VC$X2.d2+1] 
   
    start.v  <- c(b1, b2, s1, s2, n1, VC$dof.st, VC$i.rho); names(start.v) <- nstv  
    
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp    
    
    }  
    




  if(margins[1] %in% M$m3 && margins[2] %in% c(M$m2) && l.flist > 2 && M$BivD == "T"){
    
    b1 <- gamlss1$coefficients[1:VC$X1.d2]
    s1 <- gamlss1$coefficients[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]
    n1 <- gamlss1$coefficients[(VC$X1.d2+VC$X3.d2+1):(VC$X1.d2+VC$X3.d2+VC$X5.d2)]  
    b2 <- gamlss2$coefficients[1:VC$X2.d2]
    s2 <- gamlss2$coefficients[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]
   
    start.v  <- c(b1, b2, s1, s2, n1, GAM$gam6$coefficients, GAM$gam7$coefficients); names(start.v) <- nstv  
   
  if( VC$l.sp1 != 0 ) sp1 <- gamlss1$sp[1:VC$l.sp1]
  if( VC$l.sp3 != 0 ) sp3 <- gamlss1$sp[(VC$l.sp1 + 1):(VC$l.sp1 + VC$l.sp3)]
  if( VC$l.sp5 != 0 ) sp5 <- gamlss1$sp[(VC$l.sp1 + VC$l.sp3 + 1):(VC$l.sp1 + VC$l.sp3 + VC$l.sp5)]
    
  if( VC$l.sp2 != 0 ) sp2 <- gamlss2$sp[1:VC$l.sp2]
  if( VC$l.sp4 != 0 ) sp4 <- gamlss2$sp[(VC$l.sp2 + 1):(VC$l.sp2 + VC$l.sp4)]  
   
    }     
    
    
    
    
    
}    
    
  #
  #
  
sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9)  

list(start.v = start.v, sp = sp)
  
  
}


