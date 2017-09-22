Xdpred <- function(gamob, data, v.n){ 
  
    h <- .Machine$double.eps^(1/3) * ifelse(abs(data[, v.n ]) > 1, abs(data[, v.n ]), 1)
    temp <- data[, v.n ] + h
    h.hi <- temp - data[, v.n ]
    temp <- data[, v.n ] - h
    h.lo <- data[, v.n ] - temp
    twoeps <- h.hi + h.lo
    data1 <- data2 <- data
    data1[, v.n ] <- data1[, v.n ] - h
    data2[, v.n ] <- data2[, v.n ] + h
    
    data1[, v.n ] <- ifelse(data1[, v.n ] < 1e-06, 1e-06, data1[, v.n ]) 
    data2[, v.n ] <- ifelse(data2[, v.n ] < 1e-06, 1e-06, data2[, v.n ]) 
    
    attr(data1, "terms") <- attr(data2, "terms") <- NULL 
  Xd <- (predict(gamob, data2, type = "lpmatrix") - predict(gamob, data1, type = "lpmatrix"))/twoeps 
  rm(data1, data2)
  Xd 

 }
