teta.tr <- function(VC, teta.st){
 
 
 epsilon <- c(sqrt(.Machine$double.eps)) # this works well especially for AMH
 
 cjg <- c("C0","C180","C90","C270","GAL0","GAL180","GAL90","GAL270","J0","J180","J90","J270","G0","G180","G90","G270",VC$BivD2)
 
 if( VC$BivD %in% c("N","AMH","FGM","T") ) teta.st <- ifelse( abs(teta.st) > 8.75, sign(teta.st)*8.75, teta.st )  # looks fine and gives almost 1
 if(VC$BivD %in% c("AMH"))                 teta.st <- ifelse( abs(teta.st) < epsilon, epsilon, teta.st )
      
      
      
 if( VC$BivD %in% c("HO") ){
    teta.st <- ifelse( teta.st >  16,  16, teta.st )    # 0.9999999     # they look sensible
    teta.st <- ifelse( teta.st < -16, -16, teta.st )    # 1.125352e-07
 }    
      
 if( VC$BivD %in% cjg ) {
      teta.st <- ifelse( teta.st > 20, 20, teta.st )     # 709, maximum allowed # these values look fine
      teta.st <- ifelse( teta.st < -17, -17, teta.st )   # -20                  # 
                         }
     
     
     
 if(VC$BivD %in% c("HO") )                                     teta <- plogis(teta.st)                   
 if(VC$BivD %in% c("N","AMH","FGM","T") )                      teta <- tanh(teta.st)                        
 if(VC$BivD %in% c("F"))                                       teta <- ifelse( abs(teta.st) < epsilon, epsilon, teta.st )
 if(VC$BivD %in% c("C0", "C180","GAL0","GAL180","PL",VC$BivD2[c(1:4,13:16)] ))          teta <-    exp(teta.st)      
 if(VC$BivD %in% c("C90","C270","GAL90","GAL270") )                             teta <- -( exp(teta.st) )                     
 if(VC$BivD %in% c("J0", "J180","G0", "G180",VC$BivD2[5:12]) ) teta <-    exp(teta.st) + 1     
 if(VC$BivD %in% c("J90","J270","G90","G270") )                teta <- -( exp(teta.st) + 1 ) 
    
 list(teta = teta, teta.st = teta.st)   
 
   
}   
