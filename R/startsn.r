startsn <- function(margins, y1, left.trunc){
       
      
       
log.nu.1 <- NULL      

#  if( margins %in% c("GO") )  log.sig2.1 <- log((1/mean(y1))^2)  
#  if( margins %in% c("GA2") ) log.sig2.1 <- log(mean(y1)/var(y1)) 


        # GO not used anyway   
        if( !(margins %in% c("GO")) ) par.est <- try( resp.check(y1, margin = margins, plots = FALSE, print.par = TRUE, i.f = TRUE, left.trunc = left.trunc), silent = TRUE)
        if( margins %in% c("GO") )    log.sig2.1 <- log((1/mean(y1))^2)              

        
        if( !(margins %in% c("GO")) ){
        
        if(   inherits(par.est, "try-error")    ) {
        
 		if( margins %in% c("NBI","PIG","tNBI","tPIG") )  log.sig2.1 <- log( max( sqrt((var(y1) - mean(y1))/mean(y1)^2), 0.1)   )
 		if( margins %in% c("NBII","tNBII") )      log.sig2.1 <- log( max( sqrt((var(y1)/mean(y1)) - 1), 0.1)            ) 		
 		if( margins %in% c("N","LN") )            log.sig2.1 <- log( sqrt( var(y1) )                                      )  
		if( margins %in% c("LO") )                log.sig2.1 <- log( sqrt( 3*var(y1)/pi^2 )                              )   
		if( margins %in% c("IG") )                log.sig2.1 <- log( sqrt( var(y1)/mean(y1)^3 )                           )      
		if( margins %in% c("GU","rGU") )          log.sig2.1 <- log( sqrt( 6*var(y1)/pi^2 )                               )    
		if( margins %in% c("WEI") )               log.sig2.1 <- log( (1.283/sqrt(var(log(y1))))                 )    	        
		if( margins %in% c("GA","GAi") )          log.sig2.1 <- log( sqrt( var(y1)/mean(y1)^2 )                          )
        	if( margins %in% c("DAGUM","SM","FISK") ) log.sig2.1 <- log( sqrt(2)                                      )  # log(0.01) #  log(sqrt(2))       # 0.1  
        	if( margins %in% c("BE"))                 log.sig2.1 <- qlogis( sqrt( var(y1)/( mean(y1)*(1-mean(y1)) ) )         ) 
        	if( margins %in% c("DGP","DGPII", "GP","GPII","GPo")) log.sig2.1 <- log(    mean((y1 + mean(y1))/2)^2                 )      
        	
        	
        	
                                        } else log.sig2.1 <- par.est[2]
        

              if( margins %in% c("DAGUM","SM") ){
        	if(  inherits(par.est, "try-error")           ) log.nu.1 <- log(1) else log.nu.1 <- par.est[3]
                                                }                                                                                    

        }
        
list(log.sig2.1 = log.sig2.1, log.nu.1 = log.nu.1)        


}

