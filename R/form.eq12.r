form.eq12 <- function(formula.eq1, data, v1, margins, m1d, m2d, copSS = FALSE, inde = NULL){
  
    y1m <- f.eq1 <- NULL
  
    formula.eq1r <- formula.eq1   
    y1 <- y1.test <- data[, v1[1]]
    
    if(copSS == TRUE) y1 <- y1.test <- y1[inde]  
       
    if( v1[1] != as.character(formula.eq1r[2]) ) y1.test <- try(data[, as.character(formula.eq1r[2])], silent = TRUE)
    if(class(y1.test) == "try-error") stop("Please check the syntax of the equations' responses.") 

    if(margins %in% c(m1d,m2d) && min(y1.test, na.rm = TRUE) < 0) stop("The response of one or more margins must be positive.")
    if(margins %in% c(m1d,m2d)){
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if(sum(as.numeric(is.wholenumber(y1.test))) != length(y1.test)) stop("The response of one or more margins must be discrete.")     
    }
    if(margins %in% c("ZTP") && min(y1.test, na.rm = TRUE) < 1) stop("The response of one or more margins must be greater than 0.") 
    
    if(margins %in% c("probit","logit","cloglog","LN","WEI","GO","iG","GA","GAi","DAGUM","SM","FISK") && min(y1.test, na.rm = TRUE) <= 0) stop("The response of one or more margins must be positive.")
    
    if(margins %in% c("TW") && min(y1.test, na.rm = TRUE) < 0) stop("The response of one or more margins must be >= 0.")

    if(margins %in% c("BE") && (min(y1.test, na.rm = TRUE) <= 0 || max(y1.test, na.rm = TRUE) >= 1) ) stop("The response of one or more margins must be in the interval (0,1).")
     
    if( margins == "GEVlink" && length(table(y1.test))!=2 ) stop("The response must be binary.")
     
     
    # matrix useful for fitting
    
    if(margins %in% c("NBIa","NBIIa","NBI","PO","ZTP","DGP","DGPII","DGP0")){ # no PIG, NBII as these are numerical # a - all analytical
                                                                              # default for NBII all numerical (no need for y2m)
                                                                              # default for NBI half num half analyt (need y2m)
     
    ly1 <- length(y1)
    y1m <- list()
    my1 <- max(y1)
    
    if(margins != "ZTP") for(i in 1:ly1){ y1m[[i]] <- seq(0, y1[i]); length(y1m[[i]]) <- my1+1} 
    if(margins == "ZTP") for(i in 1:ly1){ y1m[[i]] <- seq(1, y1[i]); length(y1m[[i]]) <- my1} 

    
    y1m <- do.call(rbind, y1m)   
    
  
    if(max(y1) > 170 && margins %in% c("PO","ZTP","DGP0") ) y1m <- mpfr( y1m, pmax(53, getPrec(y1))) 
    
    }
     
    if( margins %in% c("N","LO","GU","rGU","GAi")      )                       formula.eq1 <- update(formula.eq1, (. + mean(.))/2 ~ . ) 
    if( margins %in% c(m1d, m2d) && margins != "GEVlink" && margins != "DGP")  formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ . )  
    if( margins %in% c("LN") )                                                 formula.eq1 <- update(formula.eq1, (log(.) + mean(log(.)))/2 ~ . )
    #if( margins %in% c("GO","GA2") )                            	       formula.eq1 <- update(formula.eq1, -(log(.) + mean(log(.)))/2 ~ . ) 
    if( margins %in% c("iG","GA","GGA","DAGUM","SM","FISK","TW") )             formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ . )    
    if( margins %in% c("WEI") )                                                formula.eq1 <- update(formula.eq1, log( exp(log(.) + 0.5772/(1.283/sqrt(var(log(.)))))  ) ~ . )     
    if( margins %in% c("BE") )                                                 formula.eq1 <- update(formula.eq1, qlogis((. + mean(.))/2) ~ . )    
    
    
    # changed 19/3/2019
    if( margins %in% c("GP","GPII","GPo","DGP","DGPII","DGP0") ) formula.eq1 <- update(formula.eq1, estobXiGP ~ . ) 
   
   
    f.eq1LI <- temp.respV ~ urcfcphmwicu # specific to surv model with L and I
       
    if( margins %in% c("probit") )  { f.eq1 <- update(formula.eq1r, . ~ urcfcphmwicu); formula.eq1 <- update(formula.eq1, -qnorm(Sh) ~ . )     }
    if( margins %in% c("logit") )   { f.eq1 <- update(formula.eq1r, . ~ urcfcphmwicu); formula.eq1 <- update(formula.eq1, -qlogis(Sh) ~ . )    } 
    if( margins %in% c("cloglog") ) { f.eq1 <- update(formula.eq1r, . ~ urcfcphmwicu); formula.eq1 <- update(formula.eq1, log(-log(Sh)) ~ . )  }   

  list(f.eq1LI = f.eq1LI, formula.eq1 = formula.eq1, formula.eq1r = formula.eq1r, y1 = y1, y1.test = y1.test, y1m = y1m, f.eq1 = f.eq1)
  
  }
  
  
  
