S.m <- function(GAM, L.SP, L.GAM, K1 = NULL){

  Ss <- list()
  off <- rank <- 0 
  i1 <- i2 <- i3 <- i4 <- i5 <- i6 <- i7 <- i8 <- 1
  k1 <- k2 <- k3 <- k4 <- k5 <- k6 <- k7 <- k8 <- 1

if (!is.null(K1)) {
  
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}




  
	for( j in 1:(L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7 + L.SP$l.sp8) ){
	
	
		if(j <= L.SP$l.sp1){     
		llg <- length(GAM$gam1$smooth[[i1]]$S) 
	                             
	                             if(llg == 1){ Ss[[j]] <- GAM$gam1$smooth[[i1]]$S[[1]]
	                                           rank[j] <- GAM$gam1$smooth[[i1]]$rank
	                                           off[j]  <- GAM$gam1$smooth[[i1]]$first.para + CLM.shift
	                                         }
	                             
	                             
	                             if(llg >  1){
	                                            Ss[[j]] <- GAM$gam1$smooth[[i1]]$S[[k1]]
	                                            rank[j] <- GAM$gam1$smooth[[i1]]$rank[k1]
	                                            off[j]  <- GAM$gam1$smooth[[i1]]$first.para + CLM.shift
		                                 } 
		                                
		             
                if(llg == 1)               i1 <- i1 + 1
                if(llg >  1 && llg == k1){ i1 <- i1 + 1; k1 <-      1; next} 
                if(llg >  1 && llg != k1)                k1 <- k1 + 1
                
                }
                
                
                
         
                if( (j >  L.SP$l.sp1 &&  j <= (L.SP$l.sp1 + L.SP$l.sp2) ) ){  
                
                llg <- length(GAM$gam2$smooth[[i2]]$S) 

                
                
                             if(llg == 1){
                     	     Ss[[j]] <- GAM$gam2$smooth[[i2]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + GAM$gam2$smooth[[i2]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam2$smooth[[i2]]$rank
                             }
                             
	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam2$smooth[[i2]]$S[[k2]]
                             off[j]  <- L.GAM$l.gam1 + GAM$gam2$smooth[[i2]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam2$smooth[[i2]]$rank[k2]

		                                 }                             
       
                if(llg == 1)               i2 <- i2 + 1
                if(llg >  1 && llg == k2){ i2 <- i2 + 1; k2 <-      1; next} 
                if(llg >  1 && llg != k2)                k2 <- k2 + 1

                }     
                     
                     

                if(j >  (L.SP$l.sp1 + L.SP$l.sp2) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3)  ){     
                
                llg <- length(GAM$gam3$smooth[[i3]]$S) 
                
                
                
                             if(llg == 1){
                     	     Ss[[j]] <- GAM$gam3$smooth[[i3]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + GAM$gam3$smooth[[i3]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam3$smooth[[i3]]$rank
                             }
 
	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam3$smooth[[i3]]$S[[k3]] 
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + GAM$gam3$smooth[[i3]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam3$smooth[[i3]]$rank[k3]

		                                 }   
 
  
                if(llg == 1)               i3 <- i3 + 1
                if(llg >  1 && llg == k3){ i3 <- i3 + 1; k3 <-      1; next} 
                if(llg >  1 && llg != k3)                k3 <- k3 + 1
                
                } 
                
                
 
                 
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4)  ){ 
                
                llg <- length(GAM$gam4$smooth[[i4]]$S) 
                
                
                             if(llg == 1){
                     	     Ss[[j]] <- GAM$gam4$smooth[[i4]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + GAM$gam4$smooth[[i4]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam4$smooth[[i4]]$rank
                             }
                             

	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam4$smooth[[i4]]$S[[k4]] 
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + GAM$gam4$smooth[[i4]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam4$smooth[[i4]]$rank[k4]

		                                 }                               

                if(llg == 1)               i4 <- i4 + 1
                if(llg >  1 && llg == k4){ i4 <- i4 + 1; k4 <-      1; next} 
                if(llg >  1 && llg != k4)                k4 <- k4 + 1
                
                }  
                
                
                
                 
      
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5)  ){  
                
                llg <- length(GAM$gam5$smooth[[i5]]$S) 
                
                
                             if(llg == 1){
                     	     Ss[[j]] <- GAM$gam5$smooth[[i5]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + GAM$gam5$smooth[[i5]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam5$smooth[[i5]]$rank
                             }

	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam5$smooth[[i5]]$S[[k5]] 
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + GAM$gam5$smooth[[i5]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam5$smooth[[i5]]$rank[k5]

		                                 }  
         
                if(llg == 1)               i5 <- i5 + 1
                if(llg >  1 && llg == k5){ i5 <- i5 + 1; k5 <-      1; next} 
                if(llg >  1 && llg != k5)                k5 <- k5 + 1
                
                }   
                
                
                
                

                              
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6)  ){        
                
                llg <- length(GAM$gam6$smooth[[i6]]$S) 
                     	    
                     	    
                     	     if(llg == 1){
                     	     Ss[[j]] <- GAM$gam6$smooth[[i6]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + GAM$gam6$smooth[[i6]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam6$smooth[[i6]]$rank
                             }
                             
	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam6$smooth[[i6]]$S[[k6]] 
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + GAM$gam6$smooth[[i6]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam6$smooth[[i6]]$rank[k6]

		                                 }                               
           
                if(llg == 1)               i6 <- i6 + 1
                if(llg >  1 && llg == k6){ i6 <- i6 + 1; k6 <-      1; next} 
                if(llg >  1 && llg != k6)                k6 <- k6 + 1
                
                
                }  
                
                
                
                

                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7)  ){        
                
                llg <- length(GAM$gam7$smooth[[i7]]$S)  
                     	     
                     	     
                     	     if(llg == 1){
                     	     Ss[[j]] <- GAM$gam7$smooth[[i7]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + L.GAM$l.gam6 + GAM$gam7$smooth[[i7]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam7$smooth[[i7]]$rank
                             }
                             
	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam7$smooth[[i7]]$S[[k7]] 
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + L.GAM$l.gam6 + GAM$gam7$smooth[[i7]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam7$smooth[[i7]]$rank[k7]

		                                 }                               
                             
                             
                if(llg == 1)               i7 <- i7 + 1
                if(llg >  1 && llg == k7){ i7 <- i7 + 1; k7 <-      1; next} 
                if(llg >  1 && llg != k7)                k7 <- k7 + 1
                
                
                }  
                
                
                
                   
                if(j >  (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7) && j <= (L.SP$l.sp1 + L.SP$l.sp2 + L.SP$l.sp3 + L.SP$l.sp4 + L.SP$l.sp5 + L.SP$l.sp6 + L.SP$l.sp7 + L.SP$l.sp8)  ){        
                
                llg <- length(GAM$gam8$smooth[[i8]]$S) 
                     	     
                     	     if(llg == 1){
                     	     Ss[[j]] <- GAM$gam8$smooth[[i8]]$S[[1]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + L.GAM$l.gam6 + L.GAM$l.gam7 + GAM$gam8$smooth[[i8]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam8$smooth[[i8]]$rank
                             }
                             
	                             if(llg >  1){

                     	     Ss[[j]] <- GAM$gam8$smooth[[i8]]$S[[k8]]  
                             off[j]  <- L.GAM$l.gam1 + L.GAM$l.gam2 + L.GAM$l.gam3 + L.GAM$l.gam4 + L.GAM$l.gam5 + L.GAM$l.gam6 + L.GAM$l.gam7 + GAM$gam8$smooth[[i8]]$first.para + CLM.shift2
                             rank[j] <- GAM$gam8$smooth[[i8]]$rank[k8]

		                                 }                               
                             
                             
                if(llg == 1)               i8 <- i8 + 1
                if(llg >  1 && llg == k8){ i8 <- i8 + 1; k8 <-      1; next} 
                if(llg >  1 && llg != k8)                k8 <- k8 + 1
                
                
                }    
                
                
             
                                                           
}


       
list(rank=rank, off=off, Ss=Ss)

}


