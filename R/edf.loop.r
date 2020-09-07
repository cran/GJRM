edf.loop <- function(VC, F, F1, GAM){

# CLM.shift accounts for the presence of the cut-points and the removed intercept in the ordinal model

if (!is.null(VC$K1)) {

        K1 <- VC$K1  
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}


edf <- edf1 <- NULL
na1 <- na2 <- na3 <- na4 <- na5 <- na6 <- na7 <- na8 <- NA

if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0 || VC$l.sp7!=0 || VC$l.sp8!=0) ){

  edf <- edf1 <- list(0, 0, 0, 0, 0, 0, 0, 0)
        
     for(i in 1:8){

       if(i==1) {mmm <- VC$lsgam1; if(mmm==0) next}
       if(i==2) {mmm <- VC$lsgam2; if(mmm==0) next} 
       if(i==3) {mmm <- VC$lsgam3; if(mmm==0) next} 
       if(i==4) {mmm <- VC$lsgam4; if(mmm==0) next} 
       if(i==5) {mmm <- VC$lsgam5; if(mmm==0) next}        
       if(i==6) {mmm <- VC$lsgam6; if(mmm==0) next} 
       if(i==7) {mmm <- VC$lsgam7; if(mmm==0) next}        
       if(i==8) {mmm <- VC$lsgam8; if(mmm==0) break}        

          for(k in 1:mmm){

              if(i==1){ gam <- GAM$gam1; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift } 
              if(i==2){ gam <- GAM$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 } 
              if(i==3){ gam <- GAM$gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 + VC$X2.d2 } 
              if(i==4){ gam <- GAM$gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 } 
              if(i==5){ gam <- GAM$gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 } 
              if(i==6){ gam <- GAM$gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 } 
              if(i==7){ gam <- GAM$gam7; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 } 
              if(i==8){ gam <- GAM$gam8; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2 } 
              
              
	      edf[[i]][k]  <-  sum(diag(F)[ind])
	      edf1[[i]][k] <- sum(diag(F1)[ind])
                        }
                  }
     
  if(VC$l.sp1!=0){ for(j in 1:VC$lsgam1) na1[j] <- GAM$gam1$smooth[[j]]$label; names(edf[[1]]) <- names(edf1[[1]]) <- na1 } 
  if(VC$l.sp2!=0){ for(j in 1:VC$lsgam2) na2[j] <- GAM$gam2$smooth[[j]]$label; names(edf[[2]]) <- names(edf1[[2]]) <- na2 }  
  if(VC$l.sp3!=0){ for(j in 1:VC$lsgam3) na3[j] <- GAM$gam3$smooth[[j]]$label; names(edf[[3]]) <- names(edf1[[3]]) <- na3 } 
  if(VC$l.sp4!=0){ for(j in 1:VC$lsgam4) na4[j] <- GAM$gam4$smooth[[j]]$label; names(edf[[4]]) <- names(edf1[[4]]) <- na4 }
  if(VC$l.sp5!=0){ for(j in 1:VC$lsgam5) na5[j] <- GAM$gam5$smooth[[j]]$label; names(edf[[5]]) <- names(edf1[[5]]) <- na5 }
  if(VC$l.sp6!=0){ for(j in 1:VC$lsgam6) na6[j] <- GAM$gam6$smooth[[j]]$label; names(edf[[6]]) <- names(edf1[[6]]) <- na6 }
  if(VC$l.sp7!=0){ for(j in 1:VC$lsgam7) na7[j] <- GAM$gam7$smooth[[j]]$label; names(edf[[7]]) <- names(edf1[[7]]) <- na7 } 
  if(VC$l.sp8!=0){ for(j in 1:VC$lsgam8) na8[j] <- GAM$gam8$smooth[[j]]$label; names(edf[[8]]) <- names(edf1[[8]]) <- na8 }
  
}


  
 
list(edf = edf, edf1 = edf1) 

}



