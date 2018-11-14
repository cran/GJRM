overall.svG <- function(formula, data, ngc, margins, M, vo, gam1, gam2, type = "copR", inde = NULL, c.gam2 = NULL, gam3 = NULL, knots = NULL){
  
X3 = X4 = X5 = X6 = X7 = X8 = X3.d2 = X4.d2 = X5.d2 = X6.d2 = X7.d2 = X8.d2 = NULL
gp3 = gp4 = gp5 = gp6 = gp7 = gp8 = NULL
gam4 = gam5 = gam6 = gam7 = gam8 = NULL
l.sp3 = l.sp4 = l.sp5 = l.sp6 = l.sp7 = l.sp8 = 0  
sp3 = sp4 = sp5 = sp6 = sp7 = sp8 = NULL  
X3s = X4s = X5s = NULL
  
if(type != "triv") gam3 <- NULL    
  
  
  
  
if(type == "triv"){

    
    formula.eq4 <- formula[[4]] 
    nad <- "theta12" 
    formula.eq4 <- as.formula( paste(nad,"~",formula.eq4[2],sep="") )
    
    set.seed(1)
    theta12 <- rnorm(vo$n, vo$theta12, 0.001)    
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    
    formula.eq5 <- formula[[5]] 
    nad <- "theta13" 
    formula.eq5 <- as.formula( paste(nad,"~",formula.eq5[2],sep="") )
    
    set.seed(1)
    theta13 <- rnorm(vo$n, vo$theta13, 0.001)    
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)       
      
    formula.eq6 <- formula[[6]] 
    nad <- "theta23" 
    formula.eq6 <- as.formula( paste(nad,"~",formula.eq6[2],sep="") )
    
    set.seed(1)
    theta23 <- rnorm(vo$n, vo$theta23, 0.001)    
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam6 <- gam(formula.eq6, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
 
    l.sp5 <- length(gam5$sp)    
    l.sp4 <- length(gam4$sp)    
    l.sp6 <- length(gam6$sp)    
    
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }
                   
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                     

                
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]    
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2] 

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf 
    
    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf   
    
    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf     
  
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients )
    
    
    
}  
  
  
  
  
if(type == "biv"){

    if(M$l.flist == 3){
    
    formula.eq3 <- formula[[3]] 
    nad <- "theta" 
    formula.eq3 <- as.formula( paste(nad,"~",formula.eq3[2],sep="") )
    
    set.seed(1)
    theta <- rnorm(vo$n, vo$i.rho, 0.001)    
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    
    
    
   if(M$Model=="BSS"){ 
    
    ######
    # TEST
    ######
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParBIV for more information.") 
    ######      
   
   }
    
    
    
    l.sp3 <- length(gam3$sp)
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    
    if(l.sp3 != 0) sp3 <- gam3$sp
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
  
  
   
if(M$Model == "BSS") start.v  <- c( gam1$coefficients, c.gam2, gam3$coefficients )           
if(M$Model != "BSS") start.v  <- c( gam1$coefficients, gam2$coefficients, gam3$coefficients )

  }
  
  
    if(M$l.flist == 4){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    nad1 <- "sigma2" 
    nad2 <- "theta" 
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") )        
  
    
    set.seed(1)
    sigma2 <- rnorm(vo$n, vo$log.sig2, 0.001) 
    theta  <- rnorm(vo$n, vo$i.rho, 0.001)           
    rm(list=".Random.seed", envir=globalenv()) 
    
    
   
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)  
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

               
    
    
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]    


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
  
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients )
    
    
    
  }  
  
  
    if(M$l.flist == 5){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]     
    nad1 <- "sigma2" 
    nad2 <- "nu" 
    nad3 <- "theta"     
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad3,"~",formula.eq5[2],sep="") )  
    
    
    set.seed(1)
    sigma2 <- rnorm(vo$n, vo$log.sig2, 0.001) 
    nu     <- rnorm(vo$n, vo$log.nu, 0.001)         
    theta  <- rnorm(vo$n, vo$i.rho, 0.001)       
    rm(list=".Random.seed", envir=globalenv()) 
    
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)  
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)    
  
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
  
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]     

    if(l.sp3 != 0) sp3 <- gam3$sp
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    if(l.sp4 != 0) sp4 <- gam4$sp
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    
    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf      
    
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients )
    
              }   
  
  
  





}












############################



if(type == "copR"){   
  

BivD <- M$BivD
if(M$surv == TRUE) BivD <- "N"  

    
    
    
    
if(margins[1] %in% c(M$m2,M$m3) && margins[2] %in% c(M$m2,M$m3) && BivD == "T"){




 if(margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2) ){

    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]     
    formula.eq6 <- formula[[6]]     
    
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.3  <- "dof"
    nad3   <- "theta"  
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.3,"~", formula.eq5[2],sep="") )     
    formula.eq6 <- as.formula( paste(  nad3,"~",formula.eq6[2],sep="") )   
    
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    dof      <- rnorm(vo$n, vo$dof.st, 0.001)  
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)      
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     

    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)
    l.sp6 <- length(gam6$sp)  
    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
    
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }       
    
    
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]     


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf      
    
    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf      
    
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients )
  
    
  }     
  
  
  
  
  
  
  if(margins[1] %in% c(M$m3) && margins[2] %in% c(M$m3)){
      
      formula.eq3 <- formula[[3]] 
      formula.eq4 <- formula[[4]] 
      formula.eq5 <- formula[[5]]   
      formula.eq6 <- formula[[6]] 
      formula.eq7 <- formula[[7]]  
      formula.eq8 <- formula[[8]]     

      nad2.1 <- "sigma2.1" 
      nad2.2 <- "sigma2.2" 
      nad.1 <- "nu.1" 
      nad.2 <- "nu.2"
      nad.3 <- "dof"
      nad3 <- "theta"     
      
      formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
      formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
      formula.eq5 <- as.formula( paste(nad.1,"~", formula.eq5[2],sep="") ) 
      formula.eq6 <- as.formula( paste(nad.2,"~", formula.eq6[2],sep="") ) 
      formula.eq7 <- as.formula( paste(nad.3,"~", formula.eq7[2],sep="") )     
      formula.eq8 <- as.formula( paste( nad3,"~",formula.eq8[2],sep="") )      
          
      set.seed(1)
      sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
      sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
      nu.1     <- rnorm(vo$n, vo$log.nu.1, 0.001)   
      nu.2     <- rnorm(vo$n, vo$log.nu.2, 0.001)     
      dof      <- rnorm(vo$n, vo$dof.st, 0.001)          
      theta    <- rnorm(vo$n, vo$i.rho, 0.001)        
      rm(list=".Random.seed", envir=globalenv()) 
      
      
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
      gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
      gam7 <- gam(formula.eq7, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam8 <- gam(formula.eq8, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    

      l.sp3 <- length(gam3$sp)   
      l.sp4 <- length(gam4$sp)    
      l.sp5 <- length(gam5$sp)   
      l.sp6 <- length(gam6$sp)    
      l.sp7 <- length(gam7$sp) 
      l.sp8 <- length(gam8$sp)  
      
      
      if(l.sp3 != 0){
      ngc <- 2
      while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     } 
                     
      if(l.sp4 != 0){
      ngc <- 2
      while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                    
  
      if(l.sp5 != 0){
      ngc <- 2
      while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                    
   
      if(l.sp6 != 0){
      ngc <- 2
      while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }   
  
      if(l.sp7 != 0){
      ngc <- 2
      while( any(round(summary(gam7)$edf, 1) > 1) ) {gam7 <- gam(formula.eq7, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }   
                     
      if(l.sp8 != 0){
      ngc <- 2
      while( any(round(summary(gam8)$edf, 1) > 1) ) {gam8 <- gam(formula.eq8, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                          
      
      
      
      
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]  
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]  
      X6 <- model.matrix(gam6)
      X6.d2 <- dim(X6)[2]
      X7 <- model.matrix(gam7)
      X7.d2 <- dim(X7)[2]   
      X8 <- model.matrix(gam8)
      X8.d2 <- dim(X8)[2]        
  
  
      if(l.sp3 != 0) sp3 <- gam3$sp 
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf 
      
  
      if(l.sp4 != 0) sp4 <- gam4$sp 
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf     
      
  
      if(l.sp5 != 0) sp5 <- gam5$sp 
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf    
      
  
      if(l.sp6 != 0) sp6 <- gam6$sp 
      environment(gam6$formula) <- environment(gam2$formula)
      gp6 <- gam6$nsdf   
      
  
      if(l.sp7 != 0) sp7 <- gam7$sp 
      environment(gam7$formula) <- environment(gam2$formula)
      gp7 <- gam7$nsdf 
      
      if(l.sp8 != 0) sp8 <- gam8$sp 
      environment(gam8$formula) <- environment(gam2$formula)
      gp8 <- gam8$nsdf         
      
      start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients, gam7$coefficients, gam8$coefficients )
    
      
  }    
  
  
  
  
  
    if(margins[1] %in% c(M$m2) && margins[2] %in% c(M$m3) ){
        
        formula.eq3 <- formula[[3]] 
        formula.eq4 <- formula[[4]] 
        formula.eq5 <- formula[[5]]   
        formula.eq6 <- formula[[6]] 
        formula.eq7 <- formula[[7]]  

  
        nad2.1 <- "sigma2.1" 
        nad2.2 <- "sigma2.2" 
        nad.2 <- "nu.2"
        nad.3 <- "dof"
        nad3 <- "theta"     
        
        formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
        formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
        formula.eq5 <- as.formula( paste(nad.2,"~", formula.eq5[2],sep="") ) 
        formula.eq6 <- as.formula( paste(nad.3,"~", formula.eq6[2],sep="") )     
        formula.eq7 <- as.formula( paste( nad3,"~", formula.eq7[2],sep="") )      
            
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)   
        nu.2     <- rnorm(vo$n, vo$log.nu.2, 0.001)     
        dof      <- rnorm(vo$n, vo$dof.st, 0.001)          
        theta    <- rnorm(vo$n, vo$i.rho, 0.001)        
        rm(list=".Random.seed", envir=globalenv()) 
        
        
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
        gam7 <- gam(formula.eq7, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)
   
  
        l.sp3 <- length(gam3$sp)   
        l.sp4 <- length(gam4$sp)    
        l.sp5 <- length(gam5$sp)   
        l.sp6 <- length(gam6$sp)    
        l.sp7 <- length(gam7$sp) 

        
        
        if(l.sp3 != 0){
        ngc <- 2
        while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                       } 
                       
        if(l.sp4 != 0){
        ngc <- 2
        while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                       }                    
    
        if(l.sp5 != 0){
        ngc <- 2
        while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                       }                    
     
        if(l.sp6 != 0){
        ngc <- 2
        while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                       }   
    
        if(l.sp7 != 0){
        ngc <- 2
        while( any(round(summary(gam7)$edf, 1) > 1) ) {gam7 <- gam(formula.eq7, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                       }   
                         
        
        
        
        
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]  
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]  
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        X7 <- model.matrix(gam7)
        X7.d2 <- dim(X7)[2]   
       
    
    
        if(l.sp3 != 0) sp3 <- gam3$sp 
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf 
        
    
        if(l.sp4 != 0) sp4 <- gam4$sp 
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf     
        
    
        if(l.sp5 != 0) sp5 <- gam5$sp 
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf    
        
    
        if(l.sp6 != 0) sp6 <- gam6$sp 
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf   
        
    
        if(l.sp7 != 0) sp7 <- gam7$sp 
        environment(gam7$formula) <- environment(gam2$formula)
        gp7 <- gam7$nsdf 
        
        
        
        start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients, gam7$coefficients)
      
        
  }    
  
  
  
  
  
  
  
 


  
  
  
  if(margins[1] %in% c(M$m3) && margins[2] %in% c(M$m2) ){
      
      formula.eq3 <- formula[[3]] 
      formula.eq4 <- formula[[4]] 
      formula.eq5 <- formula[[5]]   
      formula.eq6 <- formula[[6]] 
      formula.eq7 <- formula[[7]]  
 

      nad2.1 <- "sigma2.1" 
      nad2.2 <- "sigma2.2" 
      nad.1 <- "nu.1" 
      nad.3 <- "dof"
      nad3 <- "theta"     
      
      formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
      formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
      formula.eq5 <- as.formula( paste(nad.1,"~", formula.eq5[2],sep="") ) 
      formula.eq6 <- as.formula( paste(nad.3,"~", formula.eq6[2],sep="") )     
      formula.eq7 <- as.formula( paste( nad3,"~", formula.eq7[2],sep="") )      
          
      set.seed(1)
      sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
      sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
      nu.1     <- rnorm(vo$n, vo$log.nu.1, 0.001)     
      dof      <- rnorm(vo$n, vo$dof.st, 0.001)          
      theta    <- rnorm(vo$n, vo$i.rho, 0.001)        
      rm(list=".Random.seed", envir=globalenv()) 
      
      
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
      gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
      gam7 <- gam(formula.eq7, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)
  

      l.sp3 <- length(gam3$sp)   
      l.sp4 <- length(gam4$sp)    
      l.sp5 <- length(gam5$sp)   
      l.sp6 <- length(gam6$sp)    
      l.sp7 <- length(gam7$sp) 
 
      
      
      if(l.sp3 != 0){
      ngc <- 2
      while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     } 
                     
      if(l.sp4 != 0){
      ngc <- 2
      while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                    
  
      if(l.sp5 != 0){
      ngc <- 2
      while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                    
   
      if(l.sp6 != 0){
      ngc <- 2
      while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }   
  
      if(l.sp7 != 0){
      ngc <- 2
      while( any(round(summary(gam7)$edf, 1) > 1) ) {gam7 <- gam(formula.eq7, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }   
                     
                        
      
      
      
      
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]  
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]  
      X6 <- model.matrix(gam6)
      X6.d2 <- dim(X6)[2]
      X7 <- model.matrix(gam7)
      X7.d2 <- dim(X7)[2]   
      
  
  
      if(l.sp3 != 0) sp3 <- gam3$sp 
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf 
      
  
      if(l.sp4 != 0) sp4 <- gam4$sp 
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf     
      
  
      if(l.sp5 != 0) sp5 <- gam5$sp 
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf    
      
  
      if(l.sp6 != 0) sp6 <- gam6$sp 
      environment(gam6$formula) <- environment(gam2$formula)
      gp6 <- gam6$nsdf   
      
  
      if(l.sp7 != 0) sp7 <- gam7$sp 
      environment(gam7$formula) <- environment(gam2$formula)
      gp7 <- gam7$nsdf 
      

      start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients, gam7$coefficients)
     
  }    
  





}else{
    


    if(margins[1] %in% c(M$m1d,M$bl) && margins[2] %in% c(M$m1d,M$bl)){
    
    formula.eq3 <- formula[[3]] 
    
    nad3   <- "theta"     
    
    formula.eq3 <- as.formula( paste(  nad3,"~",formula.eq3[2],sep="") )   
    
    set.seed(1)
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)      
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    l.sp3 <- length(gam3$sp)   
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                             
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
 
    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
       
          
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients )   
    
  }   
  
  
  

  
  
    if(margins[1] %in% c(M$m1d) && margins[2] %in% c(M$m2,M$m2d)){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]]  
    
    nad2.2 <- "sigma2.2" 
    nad3   <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.2,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(  nad3,"~",formula.eq4[2],sep="") )   
    
    set.seed(1)
    sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)      
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

                  
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    
    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
          
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients )
  
    
  }  
  
  
  
  
  
  
  
   if(margins[1] %in% c(M$m1d) && margins[2] %in% c(M$m3)){
      
      formula.eq3 <- formula[[3]] 
      formula.eq4 <- formula[[4]] 
      formula.eq5 <- formula[[5]]   
      nad2.2 <- "sigma2.2" 
      nad.2 <- "nu.2"     
      nad3 <- "theta"     
      
      formula.eq3 <- as.formula( paste(nad2.2,"~",formula.eq3[2],sep="") ) 
      formula.eq4 <- as.formula( paste(nad.2,"~", formula.eq4[2],sep="") )     
      formula.eq5 <- as.formula( paste(  nad3,"~",formula.eq5[2],sep="") ) 
      
      set.seed(1)
      sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
      nu.2     <- rnorm(vo$n, vo$log.nu.2, 0.001)       
      theta    <- rnorm(vo$n, vo$i.rho, 0.001)         
      rm(list=".Random.seed", envir=globalenv()) 
      
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
      l.sp3 <- length(gam3$sp)    
      l.sp4 <- length(gam4$sp)    
      l.sp5 <- length(gam5$sp)    
      
      if(l.sp3 != 0){
      ngc <- 2
      while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     } 
                     
      if(l.sp4 != 0){
      ngc <- 2
      while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                    
  
      if(l.sp5 != 0){
      ngc <- 2
      while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                     }                    
   
 
  
       
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]  
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]  
    
  
  
      if(l.sp3 != 0) sp3 <- gam3$sp 
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf 
      
  
      if(l.sp4 != 0) sp4 <- gam4$sp 
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf     
      
  
      if(l.sp5 != 0) sp5 <- gam5$sp 
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf    
      

      start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients )
    
      
    }    


    
    
    
    
    
    


    if(margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d)){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]     
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad3   <- "theta"     
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(  nad3,"~",formula.eq5[2],sep="") )   
    
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)      
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)  
    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
    
    
    
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]     


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf      
    
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients )
  
    
  }   
  
  
  
    if(margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$bl)){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    nad2.1 <- "sigma2.1" 
    nad3   <- "theta"     
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(  nad3,"~",formula.eq4[2],sep="") )   
    
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)      
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
 
    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
     
    
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients)
  
    
  }  
  
  

    if(margins[1] %in% c(M$m3,M$m3d) && margins[2] %in% c(M$bl) ){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]      
    nad2.1 <- "sigma2.1" 
    nad.1  <- "nu.1" 
    nad3   <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad.1,"~", formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(  nad3,"~",formula.eq5[2],sep="") )      
        
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    nu.1     <- rnorm(vo$n, vo$log.nu.1, 0.001)   
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)         
    rm(list=".Random.seed", envir=globalenv()) 
    
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)   
 
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
                     
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
      


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients)
   
    
  }    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    if(margins[1] %in% c(M$m3,M$m3d) && margins[2] %in% c(M$m3,M$m3d) ){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]   
    formula.eq6 <- formula[[6]] 
    formula.eq7 <- formula[[7]]     
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.1 <- "nu.1" 
    nad.2 <- "nu.2"     
    nad3 <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.1,"~", formula.eq5[2],sep="") ) 
    formula.eq6 <- as.formula( paste(nad.2,"~", formula.eq6[2],sep="") )     
    formula.eq7 <- as.formula( paste(  nad3,"~",formula.eq7[2],sep="") )      
        
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    nu.1     <- rnorm(vo$n, vo$log.nu.1, 0.001)   
    nu.2     <- rnorm(vo$n, vo$log.nu.2, 0.001)         
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)         
    rm(list=".Random.seed", envir=globalenv()) 
    
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
    gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
    gam7 <- gam(formula.eq7, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)   
    l.sp6 <- length(gam6$sp)    
    l.sp7 <- length(gam7$sp)  
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }   

    if(l.sp7 != 0){
    ngc <- 2
    while( any(round(summary(gam7)$edf, 1) > 1) ) {gam7 <- gam(formula.eq7, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                      
    
    
    
    
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]
    X7 <- model.matrix(gam7)
    X7.d2 <- dim(X7)[2]       


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf   
    

    if(l.sp7 != 0) sp7 <- gam7$sp 
    environment(gam7$formula) <- environment(gam2$formula)
    gp7 <- gam7$nsdf       
    
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients, gam7$coefficients )
   
    
  }    
  
  

  
  
  
  
  
    if(margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m3,M$m3d) ){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]   
    formula.eq6 <- formula[[6]]    
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.2 <- "nu.2"     
    nad3 <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.2,"~", formula.eq5[2],sep="") )     
    formula.eq6 <- as.formula( paste(  nad3,"~",formula.eq6[2],sep="") ) 
    
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    nu.2     <- rnorm(vo$n, vo$log.nu.2, 0.001)        
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)       
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
    gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)    
    l.sp6 <- length(gam6$sp)  
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }   

     
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]      


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf   
         
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients )
   
    
  }    
    
  

    if(margins[1] %in% c(M$m3,M$m3d) && margins[2] %in% c(M$m2,M$m2d) ){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]   
    formula.eq6 <- formula[[6]]    
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.1 <- "nu.1"     
    nad3 <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.1,"~", formula.eq5[2],sep="") )     
    formula.eq6 <- as.formula( paste(  nad3,"~",formula.eq6[2],sep="") )     
    
    set.seed(1)
    sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    nu.1     <- rnorm(vo$n, vo$log.nu.1, 0.001)   
    theta    <- rnorm(vo$n, vo$i.rho, 0.001)       
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
    gam6 <- gam(formula.eq6, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)   
    l.sp5 <- length(gam5$sp)   
    l.sp6 <- length(gam6$sp)  
    

    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }   
  
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]      


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf   
         
    start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients, gam6$coefficients )
    
    
  }   
  
  
  
  
}  
  
  


}







if(type == "gaml"){ 

sp2 <- NULL


    if(margins %in% c(M$m2,M$m2d)){
    
    formula.eq2 <- formula[[2]] 
    nad2.1      <- "sigma2" 
    formula.eq2 <- as.formula( paste(nad2.1,"~",formula.eq2[2], sep="") ) 
  
    
    set.seed(1)
    sigma2 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam2  <- gam(formula.eq2, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)    
    l.sp2 <- length(gam2$sp)   
 
    if(l.sp2 != 0){
    ngc <- 2
    while( any(round(summary(gam2)$edf, 1) > 1) ) {gam2 <- gam(formula.eq2, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    X2    <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]

    if(l.sp2 != 0) sp2 <- gam2$sp 
    environment(gam2$formula) <- environment(gam1$formula)
    gp2 <- gam2$nsdf 
    
    start.v <- c(gam1$coefficients, gam2$coefficients )  
    
  }   
  
  
    if(margins %in% c(M$m3,M$m3d)){
    
    formula.eq2 <- formula[[2]] 
    formula.eq3 <- formula[[3]] 
   
    nad2.1 <- "sigma2" 
    nad.1  <- "nu" 
   
    formula.eq2 <- as.formula( paste(nad2.1,"~",formula.eq2[2],sep="") ) 
    formula.eq3 <- as.formula( paste(nad.1,"~", formula.eq3[2],sep="") )    
        
    set.seed(1)
    sigma2 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
    nu     <- rnorm(vo$n, vo$log.nu.1, 0.001)   
    rm(list=".Random.seed", envir=globalenv()) 
    
    
    gam2 <- gam(formula.eq2, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
   
    l.sp2 <- length(gam2$sp)   
    l.sp3 <- length(gam3$sp)    

    
    if(l.sp2 != 0){
    ngc <- 2
    while( any(round(summary(gam2)$edf, 1) > 1) ) {gam2 <- gam(formula.eq2, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    
    X2 <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]  
   

    if(l.sp2 != 0) sp2 <- gam2$sp 
    environment(gam2$formula) <- environment(gam1$formula)
    gp2 <- gam2$nsdf 
    

    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam1$formula)
    gp3 <- gam3$nsdf     
    
    start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients )  
    
  }    
     

}






if(type == "copSS"){ 


if(M$l.flist == 3){
    
    formula.eq3 <- formula[[3]] 
    nad1 <- "theta"
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    
    set.seed(1)
    theta  <- rnorm(vo$n, vo$i.rho, 0.001)           
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    
    ######
    # TEST
    ######
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######   

    l.sp3 <- length(gam3$sp)    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                                
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
  
    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
         
    start.v  <- c(gam1$coefficients, c.gam2, gam3$coefficients )
 
  }      
    
    
    if(M$l.flist == 4){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    nad1 <- "sigma2" 
    nad2 <- "theta" 
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") )        
    
    set.seed(1)
    sigma2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    theta  <- rnorm(vo$n, vo$i.rho, 0.001)           
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)  
    
    ######
    # TEST
    ######
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    X4s <- try(predict.gam(gam4, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X4s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######    
    
    
    
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]    

    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    
    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
  
    start.v  <- c(gam1$coefficients, c.gam2, gam3$coefficients, gam4$coefficients )


  }  
  
  
    if(M$l.flist == 5){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]     
    nad1 <- "sigma2" 
    nad2 <- "nu" 
    nad3 <- "theta"     
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad3,"~",formula.eq5[2],sep="") )  
    
    set.seed(1)
    sigma2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
    nu     <- rnorm(vo$n, vo$log.nu.2, 0.001)         
    theta  <- rnorm(vo$n, vo$i.rho, 0.001)       
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)  
    
    ######
    # TEST
    ######
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    X4s <- try(predict.gam(gam4, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X4s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    X5s <- try(predict.gam(gam5, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X5s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######       
    
    
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)    
  
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, subset=inde, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
  
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]     

    if(l.sp3 != 0) sp3 <- gam3$sp
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    if(l.sp4 != 0) sp4 <- gam4$sp
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    
    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf      
    
    start.v  <- c(gam1$coefficients, c.gam2, gam3$coefficients, gam4$coefficients, gam5$coefficients )
     
  }   
  
  

}







L <- list(start.v = start.v,
X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2,
X7.d2 =	X7.d2, X8.d2 =	X8.d2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8,
gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8,
l.sp3 = l.sp3, l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8,
sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8, X3s = X3s, X4s = X4s, X5s = X5s)


if(type == "gaml") {L$X2 <- X2; L$X2.d2 <- X2.d2; L$gp2 <- gp2; L$gam2 <- gam2; L$l.sp2 <- l.sp2; L$sp2 <- sp2}     


L

}

