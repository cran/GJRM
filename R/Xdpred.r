Xdpred <- function(gamob, data, v.n){ 
  
    
    h <- .Machine$double.eps^(1/2) * ifelse(abs(data[, v.n ]) > 1, abs(data[, v.n ]), 1)   
    
    # choosing (1/2) instead of, e.g., (1/3) gave closer results to grad

    toleps <- 1e-05 # too close to zero means that derivative approx may be poor
  
    temp   <- data[, v.n ] + h
    h.hi   <- temp - data[, v.n ]  # h in regular situation
    temp   <- data[, v.n ] - h
    h.lo   <- data[, v.n ] - temp  # h
    twoeps <- h.hi + h.lo          # 2h

    data1 <- data2 <- data
    
    data1[, v.n ] <- data1[, v.n ] - h
    data2[, v.n ] <- data2[, v.n ] + h
    
    if( any(data1[, v.n ] < toleps) == TRUE ) stop("The time variable is either negative or too close to zero.")
        
    attr(data1, "terms") <- attr(data2, "terms") <- NULL 
    Xd <- (predict(gamob, data2, type = "lpmatrix") - predict(gamob, data1, type = "lpmatrix"))/twoeps 
  
    # this is for debugging
    #lvn <- dim(Xd)[1] # for testing/checking
    #lcl <- dim(Xd)[2]
    #v.temp <- data[, v.n] 
    #Xdi <- matrix(NA, lvn, lcl)
    #
    # for(i in 1:lvn){
    #
    #   for(j in 1:lcl){
    #
    # matpf <- function(v.temp){
    #      datat <- data 
    #      datat[, v.n] <- v.temp
    #      predict(gamob, eq = 1, datat, type = "lpmatrix")[i, j]  
    #                             }                               
    #      Xdi[i, j] <- grad(matpf, v.temp)[i] # 0.0001 smallest allowed 
    #                  }
    # 
    # }
  
  rm(data1, data2)
  Xd 

 }
