startsn <- function(margins, y1, left.trunc){
       
      
       
log.nu.1 <- NULL      





#if (margins %in% c("BB")) {
#
#    y  <- y1[,1]
#    m  <- rowSums(y1)
#    p  <- (y + 0.5) / (m + 1)
#    mu <- p  
#
#    eps <- 1e-8
#    mu  <- pmin(pmax(mu, eps), 1 - eps)
#
#    num <- sum((y - m * mu)^2 / (m * mu * (1 - mu)))
#    den <- length(y) - 1
#    if (den <= 0) den <- length(y)
#
#    R <- num / den
#
#    m_bar <- mean(m)
#
#    rho <- (R - 1) / (m_bar - 1)
#    rho <- min(max(rho, 1e-6), 1 - 1e-6)
#    phi_start <- (1 - rho) / rho
#    if (!is.finite(phi_start) || phi_start <= 0) phi_start <- 20
#
#    log.sig2.1 <- log(phi_start)
#}
#


#if( !(margins %in% c("BB")) ) {


#  if( margins %in% c("GO") )  log.sig2.1 <- log((1/mean(y1))^2)  
#  if( margins %in% c("GA2") ) log.sig2.1 <- log(mean(y1)/var(y1)) 

        # GO not used anyway   
        if( !(margins %in% c("GO")) ) par.est <- try( resp.check(y1, margin = margins, plots = FALSE, print.par = TRUE, i.f = TRUE, left.trunc = left.trunc), silent = TRUE)
        if( margins %in% c("GO") )    log.sig2.1 <- log((1/mean(y1))^2)              

        
        if( !(margins %in% c("GO")) ){
        
        if(   inherits(par.est, "try-error")    ) {


if(margins %in% c("BB")) {

    y  <- y1[,1]
    m  <- rowSums(y1)
    p  <- (y + 0.5) / (m + 1)
    mu <- p  

    var_y <- var(y)
    bin_var <- mean(m * mu * (1 - mu))
    R <- var_y / bin_var
    m_bar <- mean(m)
    phi_start <- (m_bar - R) / (R - 1)
    if (!is.finite(phi_start) || phi_start <= 0) phi_start <- 20

    log.sig2.1 <- log(phi_start)
}

        
 		if( margins %in% c("NBI","PIG","tNBI","tPIG") )  log.sig2.1 <- log( max( sqrt((var(y1) - mean(y1))/mean(y1)^2), 0.1)   )
 		if( margins %in% c("NBII","tNBII") )      log.sig2.1 <- log( max( sqrt((var(y1)/mean(y1)) - 1), 0.1)            ) 		
 		if( margins %in% c("tN","N","LN") )            log.sig2.1 <- log( sqrt( var(y1) )                                      )  
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
        
        
        
#}        
        
        
        
        
list(log.sig2.1 = log.sig2.1, log.nu.1 = log.nu.1)        


}

