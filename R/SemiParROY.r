SemiParROY <- function(formula, data = list(), weights = NULL, subset = NULL, surv = FALSE, cens = NULL,
                       BivD1 = "N", BivD2 = "N", margins = c("probit","P","P"), 
                       dof1 = 3, dof2 = 3, left.trunc2 = 0, left.trunc3 = 0, gamlssfit = FALSE,
                       fp = FALSE, infl.fac = 1, 
                       rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                       gc.l = FALSE, parscale, extra.regI = "t", knots = NULL,
                       drop.unused.levels = TRUE,
                       min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- y2m <- y3m <- theta.fx <- c.gam2 <- Sl.sf <- y23 <- NULL  
  X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- Xd2 <- Xd3 <- NULL
  y10.y20 <- y10.y21 <- y11.y30 <- y11.y31 <- NULL
  y10.y20R <- y10.y21R <- y11.y30R <- y11.y31R <- X2s <- X3s <- NULL
  
  left.trunc1 <- left.trunc2
  left.trunc2 <- left.trunc3
  
  Model <- "ROY"
  
  hess <- TRUE
  intf <- surv.flex <- FALSE
  sp.method <- "perf"
  
  mono.sm.pos2 <- mono.sm.pos3 <- NULL

  
  # **** clean up the NULL stuff above?
  # **** make sure intf does not do anything
  # TW leave out for now as it involves even more complexity in the log.lik
  
  y00 <- y10 <- y0p <- y1p <- gam2TW <- NULL # for binary - tweedie margins 
  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- X9.d2 <- l.sp1 <- l.sp2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- i.rho1 <- i.rho2 <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- gam9 <- dof.st1 <- dof.st2 <- NULL
  
  gamlss2 <- gamlss3 <- sp1 <- sp2 <- NULL
    
  sp3 <- gp3 <- X3 <- NULL  
  sp4 <- gp4 <- X4 <- NULL  
  sp5 <- gp5 <- X5 <- NULL   
  sp6 <- gp6 <- X6 <- NULL  
  sp7 <- gp7 <- X7 <- NULL  
  sp8 <- gp8 <- X8 <- NULL  
  sp9 <- gp9 <- X9 <- NULL  
  

  log.sig1 <- log.sig2 <- log.nu1 <- log.nu2 <- NULL
  
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM","T","PL","HO","GAL0", "GAL90", "GAL180", "GAL270")
  scc  <- c("C0", "C180", "GAL0" , "GAL180","J0", "J180", "G0", "G180")
  sccn <- c("C90", "C270", "GAL90", "GAL270","J90", "J270", "G90", "G270")
    
  m2   <- c("tN","N","GU","rGU","LO","LN","WEI","IG","GA","BE","FISK","GP","GPII","GPo")
  m3   <- c("DAGUM","SM","TW")
  m1d  <- c("P", "tP","DGP0") 
  m2d  <- c("tNBI", "tNBII", "tPIG","NBI", "NBII", "PIG","DGP","DGPII") 
  bl   <- c("probit", "logit", "cloglog")   
  
  
  if(margins[2] %in% c(bl,m1d) && margins[3] %in% c(bl,m1d))  gp4 <- gp5 <- gp6 <- gp7 <- gp8 <- gp9 <- 0  
  if(margins[2] %in% c(m2,m2d) && margins[3] %in% c(m2,m2d)){ gp4 <- gp5 <- gp6 <- gp7 <- 1; gp8 <- gp9 <- 0}  
  if(margins[2] %in% c(m3)     && margins[3] %in% c(m3)    ){ gp4 <- gp5 <- gp6 <- gp7 <- gp8 <- gp9 <- 1}  
  
  
  
  M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, BivD1 = BivD1, BivD2 = BivD2, 
               opc = opc, extra.regI = extra.regI, margins = margins, bl = bl, intf = intf,
               theta.fx = theta.fx, Model = "ROY", mb = NULL, dof1 = dof1, dof2 = dof2, K1 = NULL, surv = surv)  
  

  ct  <- data.frame( c(opc),
                     c(1:14,55,56,57,60,61,62:65) 
                     )
                     
  cta <- data.frame( c(opc),
                     c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2,60,61,62:65) 
                     ) 
                     
                     
  nC1  <-  ct[which( ct[,1] == BivD1),2]
  nCa1 <- cta[which(cta[,1] == BivD1),2]     
 
  nC2  <-  ct[which( ct[,1] == BivD2),2]
  nCa2 <- cta[which(cta[,1] == BivD2),2]   
   
  
 #######################################################################################  
  
  if(!is.list(formula)) stop("You must specify a list of equations.")
  M$l.flist <- l.flist <- length(formula)
  
  pream.wm(formula, margins, M, l.flist, type = "ROY")
  form.check(formula, l.flist, ROY = TRUE)   
  

 #######################################################################################
  if(surv == TRUE) data$cens <- cens
  data$weights <- weights
  data$subset  <- subset # ok if null so no change here
  
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
            
  pred.varR <- pred.var(formula, l.flist, ROY = TRUE) 
   
  v1 <- pred.varR$v1  
  v2 <- pred.varR$v2
  v3 <- pred.varR$v3
  pred.n <- pred.varR$pred.n    
  
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
 
  mf$surv <- mf$left.trunc2 <- mf$left.trunc3 <- mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$ordinal <- mf$knots <- mf$dof1 <- mf$dof2 <- mf$intf <- mf$theta.fx <- mf$BivD1 <- mf$BivD2 <- mf$margins <- mf$fp <- mf$hess <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  if(surv == FALSE) mf$cens <- NULL
  
  mf$drop.unused.levels <- drop.unused.levels  
  mf[[1]] <- as.name("model.frame")
  
  #data <- eval(mf, parent.frame())
  data <- eval(mf)
  
  
    if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                          data$weights <- weights
                          names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"] 
 
 if(surv == TRUE){
 
    if(!("(cens)" %in% names(data))) {cens <- rep(0,dim(data)[1]) 
                          data$cens <- cens
                          names(data)[length(names(data))] <- "(cens)"}  else cens <- data[,"(cens)"]   
                          
    cens.t <- ifelse(cens == 0, 1e-07, cens)                      
 }
 

  if(gc.l == TRUE) gc()  
                                                                          
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  formula.eq3 <- formula[[3]] 

 ##############################################################  
 # Equation 1
 ##############################################################  
   
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = margins[1]), gamma=infl.fac, weights=weights, 
                              data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights))) 

  X1    <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp); if(l.sp1 != 0) sp1 <- gam1$sp
  y1    <- gam1$y
  n     <- length(y1) 
  gp1   <- gam1$nsdf 
  cy    <- 1 - y1
  inde0 <- !as.logical(y1)
  inde1 <-  as.logical(y1)
  
  n.se0 <- sum(as.numeric(inde0))
  n.se1 <- sum(as.numeric(inde1))


##############################################################
# Equations 2 and 3 for Survival outcomes 
##############################################################  

if(surv == TRUE){

 data[, v2[1]] <- ifelse(data[, v2[1]] < 0.0001, 0.0001, data[, v2[1]]) # only need once, so not also for v3[1] as variable is the same

 ###############################################################################
 ################
 # second margin
 ################


 form.eq12R   <- form.eq12(formula.eq2, data, v2, margins[2], m1d, m2d, eq1.binsurv = FALSE)   
 formula.eq2  <- form.eq12R$formula.eq1
 formula.eq2r <- form.eq12R$formula.eq1r 
 f.eq2        <- form.eq12R$f.eq1
 
 data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
 tempb <- eval(substitute(gam(f.eq2, family = cox.ph(), data = data, weights = cens, subset = inde0, drop.unused.levels = drop.unused.levels),list(cens = cens, inde0 = inde0)))
 data$Sh[inde0] <- as.vector(mm(predict(tempb, type = "response"), min.pr = min.pr, max.pr = max.pr))
  
  
 gam2  <- eval(substitute(scam(formula.eq2, gamma = infl.fac, weights = weights[inde0]*cens.t[inde0], data = data[inde0,]), list(weights = weights, cens.t = cens.t, inde0 = inde0)))
 l.sp2 <- length(gam2$sp) 
  
  
 lsgam2 <- length(gam2$smooth)
 if(lsgam2 == 0) stop("You must at least use a monotonic smooth function of time in the second equation.")
  
 clsm <- ggr <- NA 
 for(i in 1:lsgam2) clsm[i] <- class(gam2$smooth[[i]])[1] 
                    
 if( sum(as.numeric(clsm %in% c("mpi.smooth")))==0 ) stop("You must have a monotonic smooth of time, mpi, in the second equation.")
 pos.mpi <- which(clsm == "mpi.smooth")
  
     
 ###########################################################  
 
 if(l.sp2 != 0) sp2 <- gam2$sp 
 sp2[pos.mpi] <- 1 

 gam.call <- gam2$call
 gam.call$sp <- sp2
 gam2 <- eval(gam.call)
  
 ###########################################################

  for(i in 1:lsgam2){ 
    if( max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) != 0 && clsm[i] == "mpi.smooth" ) mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para) )    
                    }
       
  X2  <- predict(gam2, type = "lpmatrix")
  Xd2 <- Xdpred(gam2, data[inde0,], v2[1])

  gam2$y <- data[inde0, v2[1]]
  gam2$formula <- formula.eq2r 
 
  attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
  gp2 <- gam2$nsdf 
  X2.d2 <- dim(X2)[2]

  X2s <- try(predict(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
  if(any(class(X2s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
  
  y2 <- gam2$y 
  
 ###############################################################################
 ###############################################################################
 # third margin
 ################

 form.eq12R   <- form.eq12(formula.eq3, data, v3, margins[3], m1d, m2d, eq1.binsurv = FALSE)   
 formula.eq3  <- form.eq12R$formula.eq1
 formula.eq3r <- form.eq12R$formula.eq1r 
 f.eq3        <- form.eq12R$f.eq1
 
 data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
 tempb <- eval(substitute(gam(f.eq3, family = cox.ph(), data = data, weights = cens, subset = inde1, drop.unused.levels = drop.unused.levels),list(cens = cens, inde1 = inde1)))
 data$Sh[inde1] <- as.vector(mm(predict(tempb, type = "response"), min.pr = min.pr, max.pr = max.pr))
  
  
 gam3  <- eval(substitute(scam(formula.eq3, gamma = infl.fac, weights = weights[inde1]*cens.t[inde1], data = data[inde1,]), list(weights = weights, cens.t = cens.t, inde1 = inde1)))
 l.sp3 <- length(gam3$sp) 
  
  
 lsgam3 <- length(gam3$smooth)
 if(lsgam3 == 0) stop("You must at least use a monotonic smooth function of time in the third equation.")
  
 clsm <- ggr <- NA 
 for(i in 1:lsgam3) clsm[i] <- class(gam3$smooth[[i]])[1] 
                    
 if( sum(as.numeric(clsm %in% c("mpi.smooth")))==0 ) stop("You must have a monotonic smooth of time, mpi, in the third equation.")
 pos.mpi <- which(clsm == "mpi.smooth")
  
     
 ###########################################################  
 
 if(l.sp3 != 0) sp3 <- gam3$sp 
 sp3[pos.mpi] <- 1 

 gam.call <- gam3$call
 gam.call$sp <- sp3
 gam3 <- eval(gam.call)
  
 ###########################################################

  for(i in 1:lsgam3){ 
    if( max(as.numeric(grepl(v3[1], gam3$smooth[[i]]$term))) != 0 && clsm[i] == "mpi.smooth" ) mono.sm.pos3 <- c(mono.sm.pos3, c(gam3$smooth[[i]]$first.para:gam3$smooth[[i]]$last.para) )    
                    }
       
  X3  <- predict(gam3, type = "lpmatrix")
  Xd3 <- Xdpred(gam3, data[inde1,], v3[1])

  gam3$y <- data[inde1, v3[1]]
  gam3$formula <- formula.eq3r 
 
  attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
  gp3 <- gam3$nsdf 
  X3.d2 <- dim(X3)[2]

  X3s <- try(predict(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
  if(any(class(X3s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
  
  y3 <- gam3$y 
  
############################
# SET UP indicator variables 
############################

# cens = 1 event occured, 0 event did not occur
# cbind(y10.y20, y10.y21,y11.y30,y11.y31)[1:10,] for testing
# cbind(VC$y10.y20, VC$y10.y21,VC$y11.y30,VC$y11.y31)[1:10,]

    y10.y20 <- y10.y21 <- y11.y30 <- y11.y31 <- rep(FALSE, n)  
    
    y10.y20[inde0] <- y10.y20R <- as.logical( (1 - y1)[inde0]*(1 - cens[inde0]) )
    y10.y21[inde0] <- y10.y21R <- as.logical( (1 - y1)[inde0]*cens[inde0]       )
    y11.y30[inde1] <- y11.y30R <- as.logical( y1[inde1]*(1 - cens[inde1])       )
    y11.y31[inde1] <- y11.y31R <- as.logical( y1[inde1]*cens[inde1]             )


}




 ##############################################################
 # Equations 2 and 3 for continuous/discrete response 
 ##############################################################  

 if( margins[2] %in% bl && surv == FALSE){

 gam2 <- eval(substitute(gam(formula.eq2, binomial(link = margins[2]), gamma=infl.fac, weights=weights, 
                             data=data, subset=inde0, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights,inde0=inde0)))  
   
    X2    <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp); if(l.sp2 != 0) sp2 <- gam2$sp
    gp2   <- gam2$nsdf
    y2 <- gam2$y

    #* this will be useful for prediction, to calculate effects, not required for fitting *#
    
    X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X2s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
    
  }
 
 
  if( margins[3] %in% bl && surv == FALSE){
 
  gam3 <- eval(substitute(gam(formula.eq3, binomial(link = margins[3]), gamma=infl.fac, weights=weights, 
                              data=data, subset=inde1, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights,inde1=inde1)))  
    
    X3    <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    l.sp3 <- length(gam3$sp); if(l.sp3 != 0) sp3 <- gam3$sp
    gp3   <- gam3$nsdf
    y3 <- gam3$y
    
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X3s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
  
    y10.y20 <- y10.y21 <- y11.y30 <- y11.y31 <- rep(FALSE, n)  
    
    y10.y20[inde0] <- y10.y20R <- as.logical( (1 - y1)[inde0]*(1 - y2) )
    y10.y21[inde0] <- y10.y21R <- as.logical( (1 - y1)[inde0]*y2       )
    y11.y30[inde1] <- y11.y30R <- as.logical( y1[inde1]*(1 - y3)       )
    y11.y31[inde1] <- y11.y31R <- as.logical( y1[inde1]*y3             )
    
    
  }
  
  


  if( !(margins[2] %in% bl) ){
  
    form.eq12R   <- form.eq12(formula.eq2, data[inde0, ], v2, margins[2], m1d, m2d)   
    formula.eq2  <- form.eq12R$formula.eq1
    formula.eq2r <- form.eq12R$formula.eq1r
    y2           <- form.eq12R$y1
    y2.test      <- form.eq12R$y1.test 
    y2m          <- form.eq12R$y1m  
    
    gam2         <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data, subset=inde0, knots = knots, drop.unused.levels = drop.unused.levels),list(weights = weights, inde0 = inde0)))
    gam2$formula <- formula.eq2r  
    names(gam2$model)[1] <- as.character(formula.eq2r[2])

    y2 <- y2.test  
    if( margins[2] %in% c("LN") ) y2 <- log(y2) 

    X2    <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp); if(l.sp2 != 0) sp2 <- gam2$sp
    gp2   <- gam2$nsdf
    
    #* this will be useful for prediction, to calculate effects, not required for fitting *#
    
    X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X2s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
    
  } 
  
  
  if( !(margins[3] %in% bl) ){
  
    form.eq12R   <- form.eq12(formula.eq3, data[inde1, ], v3, margins[3], m1d, m2d) 
    formula.eq3  <- form.eq12R$formula.eq1
    formula.eq3r <- form.eq12R$formula.eq1r
    y3           <- form.eq12R$y1
    y3.test      <- form.eq12R$y1.test 
    y3m          <- form.eq12R$y1m  
    

    gam3         <- eval(substitute(gam(formula.eq3, gamma=infl.fac, weights=weights, data=data, subset=inde1, knots = knots, drop.unused.levels = drop.unused.levels),list(weights = weights, inde1 = inde1)))
    gam3$formula <- formula.eq3r  
    names(gam3$model)[1] <- as.character(formula.eq3r[2])

    y3 <- y3.test  
    if( margins[3] %in% c("LN") ) y3 <- log(y3) 

    X3    <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    l.sp3 <- length(gam3$sp); if(l.sp3 != 0) sp3 <- gam3$sp
    gp3   <- gam3$nsdf
    
    #* this will be useful for prediction, to calculate effects, not required for fitting *#
    
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X3s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
    
  }   
  
  
  
X4 <- X5 <- X6 <- X7 <- X8 <- X9 <- X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- matrix(1, n, 1)  
X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- X9.d2 <- 1 


if( margins[2] %in% c(m1d, bl) && margins[3] %in% c(m1d, bl) ){
  X4 <- as.matrix(X4[inde0,])
  X5 <- as.matrix(X5[inde1,])
}
   
if( margins[2] %in% c(m2, m2d) && margins[3] %in% c(m2, m2d) ){
  X4 <- as.matrix(X4[inde0,])
  X5 <- as.matrix(X5[inde1,])
  X6 <- as.matrix(X6[inde0,])
  X7 <- as.matrix(X7[inde1,])    
}   

if( margins[2] %in% c(m3) && margins[3] %in% c(m3) ){
  X4 <- as.matrix(X4[inde0,])
  X5 <- as.matrix(X5[inde1,])
  X6 <- as.matrix(X6[inde0,])
  X7 <- as.matrix(X7[inde1,]) 
  X8 <- as.matrix(X8[inde0,])
  X9 <- as.matrix(X9[inde1,])       
}  

if( margins[2] %in% c(m1d) && margins[3] %in% c(m2d) ){
  X4 <- as.matrix(X4[inde1,])
  X5 <- as.matrix(X5[inde0,])
  X6 <- as.matrix(X6[inde1,])     
} 


if( margins[2] %in% c(m2d) && margins[3] %in% c(m1d) ){
  X4 <- as.matrix(X4[inde0,])
  X5 <- as.matrix(X5[inde0,])
  X6 <- as.matrix(X6[inde1,])
}  

if( margins[2] %in% c(m2) && margins[3] %in% c(m3) ){
  X4 <- as.matrix(X4[inde0,])
  X5 <- as.matrix(X5[inde1,])
  X6 <- as.matrix(X6[inde1,])
  X7 <- as.matrix(X7[inde0,])  
  X8 <- as.matrix(X8[inde1,])   
} 

if( margins[2] %in% c(m3) && margins[3] %in% c(m2) ){
  X4 <- as.matrix(X4[inde0,])
  X5 <- as.matrix(X5[inde1,])
  X6 <- as.matrix(X6[inde0,])
  X7 <- as.matrix(X7[inde0,])  
  X8 <- as.matrix(X8[inde1,]) 
} 






##############################################################
# Starting values for dependence parameter
##############################################################
res1 <- residuals(gam1); res1 <- res1 + rnorm(length(res1), sd = 0.01)
res2 <- residuals(gam2); res2 <- res2 + rnorm(length(res2), sd = 0.01)
res3 <- residuals(gam3); res3 <- res3 + rnorm(length(res3), sd = 0.01)

ass.s  <- cor(res1[inde0], res2)#, method = "kendall")
ass.s  <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))                
i.rho1 <- ass.dp(ass.s, BivD1, scc, sccn, nCa1)

ass.s  <- cor(res1[inde1], res3)#, method = "kendall")
ass.s  <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))                
i.rho2 <- ass.dp(ass.s, BivD2, scc, sccn, nCa2)

names(i.rho1) <- "theta12.star"
names(i.rho2) <- "theta13.star"

               
##############################################################
# Starting values for whole parameter vector
##############################################################
          

if( margins[2] %in% c(m1d, bl) && margins[3] %in% c(m1d, bl) ) start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, i.rho1, i.rho2)
      

if( margins[2] %in% c(m2, m2d) && margins[3] %in% c(m2, m2d) ){

   start.snR1 <- startsn(margins[2], y2, left.trunc = left.trunc1)  
   start.snR2 <- startsn(margins[3], y3, left.trunc = left.trunc2)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma2.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma3.star"

start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, i.rho1, i.rho2)
    
}   


if( margins[2] %in% c(m3) && margins[3] %in% c(m3) ){

   start.snR1 <- startsn(margins[2], y2, left.trunc = left.trunc1)  
   start.snR2 <- startsn(margins[3], y3, left.trunc = left.trunc2)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma2.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma3.star"
   
   log.nu1    <- start.snR1$log.nu.1; names(log.nu1) <- "nu2.star"
   log.nu2    <- start.snR2$log.nu.1; names(log.nu2) <- "nu3.star"
   
start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, log.nu1, log.nu2, i.rho1, i.rho2)
    
}  


if( margins[2] %in% c(m1d) && margins[3] %in% c(m2d) ){

   start.snR2 <- startsn(margins[3], y3, left.trunc = left.trunc2)
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma3.star"

start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig2, i.rho1, i.rho2)
    
} 


if( margins[2] %in% c(m2d) && margins[3] %in% c(m1d) ){

   start.snR1 <- startsn(margins[2], y2, left.trunc = left.trunc1)  
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma2.star"

start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, i.rho1, i.rho2)
    
}  


if( margins[2] %in% c(m2) && margins[3] %in% c(m3) ){

   start.snR1 <- startsn(margins[2], y2, left.trunc = left.trunc1)  
   start.snR2 <- startsn(margins[3], y3, left.trunc = left.trunc2)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma2.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma3.star"
   
   log.nu2    <- start.snR2$log.nu.1; names(log.nu2) <- "nu3.star"
   
start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, log.nu2, i.rho1, i.rho2)
    
} 


if( margins[2] %in% c(m3) && margins[3] %in% c(m2) ){

   start.snR1 <- startsn(margins[2], y2, left.trunc = left.trunc1)  
   start.snR2 <- startsn(margins[3], y3, left.trunc = left.trunc2)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma2.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma3.star"
   
   log.nu1    <- start.snR1$log.nu.1; names(log.nu1) <- "nu2.star"
   
start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, log.nu1, i.rho1, i.rho2)
    
}  



##############################################################  
  
if(l.flist > 3){  


### *** ONLY DONE MAIN COMBINATIONS FOR THE TIME BEING 16/7/2021 ***###
 
vo <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, i.rho1 = i.rho1, i.rho2 = i.rho2, 
           log.sig1 = log.sig1, log.sig2 = log.sig2, log.nu1 = log.nu1, log.nu2 = log.nu2, n = n, 
           inde0 = inde0, inde1 = inde1, drop.unused.levels = drop.unused.levels)  
  
    overall.svGR <- overall.svG(formula, data, ngc = 2, margins, M, vo, gam1, gam2, type = "ROY", inde = inde, c.gam2 = c.gam2, gam3 = gam3, knots = knots)

    start.v = overall.svGR$start.v 
    
    X4 = overall.svGR$X4
    X5 = overall.svGR$X5
    X6 = overall.svGR$X6
    X7 = overall.svGR$X7 
    X8 = overall.svGR$X8
    X9 = overall.svGR$X9
    
    X4.d2 = overall.svGR$X4.d2
    X5.d2 = overall.svGR$X5.d2
    X6.d2 = overall.svGR$X6.d2
    X7.d2 = overall.svGR$X7.d2 
    X8.d2 = overall.svGR$X8.d2
    X9.d2 = overall.svGR$X9.d2
    
    gp4 = overall.svGR$gp4
    gp5 = overall.svGR$gp5
    gp6 = overall.svGR$gp6
    gp7 = overall.svGR$gp7 
    gp8 = overall.svGR$gp8
    gp9 = overall.svGR$gp9
    
    gam4 = overall.svGR$gam4
    gam5 = overall.svGR$gam5
    gam6 = overall.svGR$gam6
    gam7 = overall.svGR$gam7 
    gam8 = overall.svGR$gam8
    gam9 = overall.svGR$gam9
    
    l.sp4 = overall.svGR$l.sp4
    l.sp5 = overall.svGR$l.sp5
    l.sp6 = overall.svGR$l.sp6
    l.sp7 = overall.svGR$l.sp7
    l.sp8 = overall.svGR$l.sp8
    l.sp9 = overall.svGR$l.sp9
    
    sp4 = overall.svGR$sp4
    sp5 = overall.svGR$sp5
    sp6 = overall.svGR$sp6
    sp7 = overall.svGR$sp7
    sp8 = overall.svGR$sp8
    sp9 = overall.svGR$sp9
    
    X4s = overall.svGR$X4s
    X5s = overall.svGR$X5s 
    X6s = overall.svGR$X6s
    X7s = overall.svGR$X7s
    X8s = overall.svGR$X8s
    X9s = overall.svGR$X9s
  
}  
  
  
  
  
##########################################################
  
GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, gam9 = gam9)   


if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0 || l.sp9!=0) && fp==FALSE ){ 

L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), l.gam4 = length(gam4$coefficients),
              l.gam5 = length(gam5$coefficients), l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), l.gam8 = length(gam8$coefficients),
              l.gam9 = length(gam9$coefficients))

L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9)

                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9)
                 qu.mag <- S.m(GAM, L.SP, L.GAM)                   
                 
                                                        }  
  
##########################################################


if(missing(parscale)) parscale <- 1   

  respvec <- list(y1 = y1, y2 = y2, y3 = y3, cy = cy, univ = 0)

  my.env <- new.env()
  my.env$signind <- 1

  lsgam1 <- length(gam1$smooth)
  lsgam2 <- length(gam2$smooth)
  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)
  lsgam9 <- length(gam9$smooth)


  VC <- list(lsgam1 = lsgam1, robust = FALSE, sp.fixed = NULL, K1 = NULL, 
             lsgam2 = lsgam2, Sl.sf = Sl.sf, sp.method = sp.method,
             lsgam3 = lsgam3, left.trunc1 = left.trunc1, left.trunc2 = left.trunc2,
             lsgam4 = lsgam4, mono.sm.pos2 = mono.sm.pos2, mono.sm.pos3 = mono.sm.pos3, Xd3 = Xd3, Xd2 = Xd2,
             lsgam5 = lsgam5,
             lsgam6 = lsgam6,
             lsgam7 = lsgam7,
             lsgam8 = lsgam8, 
             lsgam9 = lsgam9,
             X1 = X1, inde0 = inde0, n.se0 = n.se0, inde1 = inde1, n.se1 = n.se1, my.env=my.env,
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,
             X7 = X7, 
             X8 = X8, 
             X9 = X9,
             X2s = X2s, 
             X3s = X3s,
             X4s = X4s, 
             X5s = X5s, 
             X6s = X6s,
             X7s = X7s, 
             X8s = X8s, 
             X9s = X9s,             
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,       
             X7.d2 = X7.d2,
             X8.d2 = X8.d2, 
             X9.d2 = X9.d2,
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6,  
             gp7 = gp7, 
             gp8 = gp8, 
             gp9 = gp9,  
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6,    
             l.sp7 = l.sp7,
             l.sp8 = l.sp8, 
             l.sp9 = l.sp9, 
             infl.fac = infl.fac,
             weights = weights,
             fp = fp, univ.gamls = FALSE,
             hess = hess, nCa1 = nCa1, nCa2 = nCa2,
             Model = Model, gamlssfit = gamlssfit,
             end = end,
             BivD = BivD1, BivD1 = BivD1, BivD2 = BivD2, dof.st1 = log(dof1 - 2), dof.st2 = log(dof2 - 2), dof1 = dof1, dof2 = dof2, 
             nC1 = nC1, nC2 = nC2, gc.l = gc.l, n = n, n.se0 = n.se0, n.se1 = n.se1, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "NO", ccss = "no", m2 = m2, m3 = m3, m2d = m2d, m1d = m1d, m3d = NULL, bl = bl,
             X2s = X2s, X3s = X3s, 
             triv = FALSE, y2m = y2m, y3m = y3m,
             theta.fx = theta.fx, i.rho1 = i.rho1, i.rho2 = i.rho2, 
             cta = cta, ct = ct, zerov = -10, surv.flex = surv.flex, gp2.inf = NULL,
             informative = "no", 
             zero.tol = 1e-02,
             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, 
             y00 = y00, y10 = y10, y0p = y0p, y1p = y1p,
             l.flist = l.flist,
             y10.y20 = y10.y20,
             y10.y21 = y10.y21,
             y11.y30 = y11.y30,
             y11.y31 = y11.y31,       
             y10.y20R = y10.y20R,             
             y10.y21R = y10.y21R,            
             y11.y30R = y11.y30R,            
             y11.y31R = y11.y31R, end.surv = FALSE, surv = surv ) # original n only needed in SemiParBIV.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################

# never thought of it, but subset may cause trouble here? don't think so

if(gamlssfit == TRUE && margins[2] %in% c(m1d, m2d, m2, m3) && margins[3] %in% c(m1d, m2d, m2, m3)){ 

  nstv <- names(start.v) 
  
  form.gamlR <- form.gaml(formula, l.flist, M, type = "ROY")

  gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = inde0,  
                   family = margins[2], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, drop.unused.levels = drop.unused.levels), list(inde0 = inde0, weights = weights)))   
                   
  gamlss3 <- eval(substitute(gamlss(form.gamlR$formula.gamlss3, data = data, weights = weights, subset = inde1,  
                   family = margins[3], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, drop.unused.levels = drop.unused.levels), list(inde1 = inde1, weights = weights)))                    
                      
  # updated starting values                   
  
  SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8, sp9 = sp9)
  
  gamls.upsvR <- gamls.upsv(gamlss1 = gamlss2, gamlss2 = gamlss3, margins, M, l.flist, nstv = NULL, VC, GAM, SP, type = "ROY")
  
  
  main.sv                    <- c(gam1$coefficients, gamls.upsvR$start.v) 
  start.v[1:length(main.sv)] <- main.sv
  names(start.v)             <- nstv                                           # this really needed? does not hurt
  
  sp <- gamls.upsvR$sp
 
  
  }  
  
  
  if(surv == TRUE){ # should also work with varying theta although never tried, also in the case above, I think it was tested at some point
  
    nstv <- names(start.v) 
    
    form.gamlR <- form.gaml(formula, l.flist, M, type = "ROY")
  
    marginss <- margins
    marginss[2] <- paste("-", marginss[2], sep = "")
    marginss[3] <- paste("-", marginss[3], sep = "")
  
    gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = inde0,  
                     family = marginss[2], infl.fac = infl.fac, cens = cens,
                     rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                     gc.l = gc.l, parscale = 1, drop.unused.levels = drop.unused.levels), list(inde0 = inde0, weights = weights, cens = cens)))   
                     
    gamlss3 <- eval(substitute(gamlss(form.gamlR$formula.gamlss3, data = data, weights = weights, subset = inde1,  
                     family = marginss[3], infl.fac = infl.fac, cens = cens,
                     rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                     gc.l = gc.l, parscale = 1, drop.unused.levels = drop.unused.levels), list(inde1 = inde1, weights = weights, cens = cens)))                    
                        
    # updated starting values                   
    
    SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8, sp9 = sp9)
    
    gamls.upsvR <- gamls.upsv(gamlss1 = gamlss2, gamlss2 = gamlss3, margins, M, l.flist, nstv = NULL, VC, GAM, SP, type = "ROY")
    
    
    main.sv                    <- c(gam1$coefficients, gamls.upsvR$start.v) 
    start.v[1:length(main.sv)] <- main.sv
    names(start.v)             <- nstv                                           # this really needed? does not hurt
    
    sp <- gamls.upsvR$sp
   
    
  }
  
  
  
  
  
  ##########################################################################################################################
  ##########################################################################################################################

  # joint.probs, not relevant for now as we are just interested in calculating effects for the time being


  func.opt <- func.OPT(margins, M, type = "ROY")

  SemiParFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
    
  
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParROY.fit.post(SemiParFit = SemiParFit, Model = Model, VC = VC, GAM)
                                            
  SemiParFit <- SemiParFit.p$SemiParFit # this may be useful but not sure yet at this stage
 
  y2.m <- y2  
  y3.m <- y3  

  if(margins[2] == "LN")  y2.m <- exp(y2)
  if(margins[3] == "LN")  y3.m <- exp(y3)

  ##########################################################################################################################


if(gc.l == TRUE) gc()

  ##########################################################################################################################

cov.c(SemiParFit)


  ##########################################################################################################################

  gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- gam9$call$data <- cl$data 
  # for all.terms
  ##########################################################################################################################


#** not sure yet I need the stuff below, will figure out when calculating effects?

if( !(Model=="B" && !(margins[2] %in% bl) && end == 2) ) {dataset <- NULL; rm(data) } else { attr(data, "terms") <- NULL; dataset <- data; rm(data) } 


L <- list(fit = SemiParFit$fit, dataset = dataset, formula = formula, SemiParFit = SemiParFit, 
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, robust = FALSE,
          gam7 = gam7, gam8 = gam8, gam9 = gam9, gam2TW = gam2TW,
          coefficients = SemiParFit$fit$argument, coef.t = SemiParFit.p$coef.t, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, left.trunc1 = left.trunc1, left.trunc2 = left.trunc2,
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta12 = SemiParFit.p$theta12, theta13 = SemiParFit.p$theta13, 
          theta12.a = SemiParFit.p$theta12.a, theta13.a = SemiParFit.p$theta13.a,   
          n = n, n.se0 = n.se0, n.se1 = n.se1, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X9 = X9, X4.d2 = X4.d2, X5.d2 = X5.d2, 
          X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, X9.d2 = X9.d2,           
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf11 = SemiParFit.p$edf11,
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8, edf9 = SemiParFit.p$edf9,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6,
          edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, edf1.9 = SemiParFit.p$edf1.9,
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p1 = SemiParFit$fit$p1, p0 = SemiParFit$fit$p0,
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, eta3 = SemiParFit$fit$eta3,           
          etad12 = SemiParFit$fit$etad1, etad13 = SemiParFit$fit$etad2, etas2 = SemiParFit$fit$etas1, etas3 = SemiParFit$fit$etas2,
          etan2 = SemiParFit$fit$etan1, etan3 = SemiParFit$fit$etan2,
          y1 = y1, y2 = y2.m, y3 = y3.m, 
          margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC1 = nC1, nC2 = nC2, hess = hess, 
          respvec = respvec, inde0 = inde0, inde1 = inde1, 
          qu.mag = qu.mag,
          sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a,
          sigma3 = SemiParFit.p$sigma3, sigma3.a = SemiParFit.p$sigma3.a,
          nu2 = SemiParFit.p$nu2, nu2.a = SemiParFit.p$nu2.a, 
          nu3 = SemiParFit.p$nu3, nu3.a = SemiParFit.p$nu3.a, 
          tau12 = SemiParFit.p$tau12, tau12.a = SemiParFit.p$tau12.a, 
          tau13 = SemiParFit.p$tau13, tau13.a = SemiParFit.p$tau13.a,
          Vb.t = SemiParFit.p$Vb.t,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, gp9 = gp9, 
          X2s = X2s, X3s = X3s, X4s = X4s, X5s = X5s, X6s = X6s, X7s = X7s, X8s = X8s, X9s = X9s, 
          VC = VC, Model = Model, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "NO",  
          l.flist = l.flist, v1 = v1, v2 = v2, v3 = v3, triv = FALSE, univar.gamlss = FALSE,
          gamlss2 = gamlss2, gamlss3 = gamlss3, BivD1 = BivD1, BivD2 = BivD2, dof12 = dof1, dof13 = dof2, 
          dof12.a = dof1, dof13.a = dof2, call = cl,
          surv = surv, surv.flex = surv.flex, end.surv = FALSE)

class(L) <- c("SemiParROY", "SemiParBIV", "gjrm") 

L

}

