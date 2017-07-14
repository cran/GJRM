teta.tr <- function(VC, teta.st){
 
 epsilon <- 0.0000001 
 cjg <- c("C0","C180","C90","C270","J0","J180","J90","J270","G0","G180","G90","G270",VC$BivD2)
 
 if( VC$BivD %in% c("N","AMH","FGM","T") ) teta.st <- ifelse( abs(teta.st) > 8.75, sign(teta.st)*8.75, teta.st )  
      
      
      
 if( VC$BivD %in% c("HO") ){
    teta.st <- ifelse( teta.st >  16,  16, teta.st )  
    teta.st <- ifelse( teta.st < -16, -16, teta.st )    
 }    
      
 if( VC$BivD %in% cjg ) {
      teta.st <- ifelse( teta.st > 28, 28, teta.st )     # 709, maximum allowed
      teta.st <- ifelse( teta.st < -17, -17, teta.st )   # -20
                         }
     
     
     
 if(VC$BivD %in% c("HO") )                                     teta <- plogis(teta.st)                   
 if(VC$BivD %in% c("N","AMH","FGM","T") )                      teta <- tanh(teta.st)                        
 if(VC$BivD == "F")                                            teta <- ifelse( abs(teta.st) < epsilon, epsilon, teta.st )
 if(VC$BivD %in% c("C0", "C180","PL",VC$BivD2[1:4] ))          teta <-    exp(teta.st)      
 if(VC$BivD %in% c("C90","C270") )                             teta <- -( exp(teta.st) )                     
 if(VC$BivD %in% c("J0", "J180","G0", "G180",VC$BivD2[5:12]) ) teta <-    exp(teta.st) + 1     
 if(VC$BivD %in% c("J90","J270","G90","G270") )                teta <- -( exp(teta.st) + 1 ) 
    
 list(teta = teta, teta.st = teta.st)   
 
   
}   
