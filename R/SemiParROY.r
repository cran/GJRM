SemiParROY <- function(formula, data = list(), weights = NULL, subset = NULL,
                       BivD1 = "N", BivD2 = "N", margins = c("probit","PO","PO"), 
                       dof1 = 3, dof2 = 3, gamlssfit = FALSE,
                       fp = FALSE, infl.fac = 1, 
                       rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                       gc.l = FALSE, parscale, extra.regI = "t", knots = NULL,
                       drop.unused.levels = TRUE,
                       min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- y2m <- y3m <- theta.fx <- c.gam2 <- Sl.sf <- y23 <- NULL  
  
  X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- NULL
  
  y10.y20 <- y10.y21 <- y11.y30 <- y11.y31 <- NULL
  
  y10.y20R <- y10.y21R <- y11.y30R <- y11.y31R <- NULL
  
  Model <- "ROY"
  
  hess <- TRUE
  intf <- surv.flex <- FALSE
  sp.method <- "perf"

  
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
    
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK","GP","GPII","GPo")
  m3   <- c("DAGUM","SM","TW")
  m1d  <- c("PO", "ZTP","DGP0") 
  m2d  <- c("NBI", "NBII", "PIG","DGP","DGPII") 
  bl   <- c("probit", "logit", "cloglog")   
  
  M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, BivD1 = BivD1, BivD2 = BivD2, 
               opc = opc, extra.regI = extra.regI, margins = margins, bl = bl, intf = intf,
               theta.fx = theta.fx, Model = "ROY", mb = NULL, dof1 = dof1, dof2 = dof2, K1 = NULL)  
                 
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
  mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$ordinal <- mf$knots <- mf$dof1 <- mf$dof2 <- mf$intf <- mf$theta.fx <- mf$BivD1 <- mf$BivD2 <- mf$margins <- mf$fp <- mf$hess <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  mf$drop.unused.levels <- drop.unused.levels  
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  

  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
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

 ##############################################################
 # Equations 2 and 3 for continuous/discrete response 
 ##############################################################  

 if( margins[2] %in% bl ){

 gam2 <- eval(substitute(gam(formula.eq2, binomial(link = margins[2]), gamma=infl.fac, weights=weights, 
                             data=data, subset=inde0, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights,inde0=inde0)))  
   
    X2    <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp); if(l.sp2 != 0) sp2 <- gam2$sp
    gp2   <- gam2$nsdf
    y2 <- gam2$y

    
    n.se0 <- sum(as.numeric(inde0))

    #* this will be useful for prediction, to calculate effects, not required for fitting *#
    
    X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X2s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
    
  }
 
 
  if( margins[3] %in% bl ){
 
  gam3 <- eval(substitute(gam(formula.eq3, binomial(link = margins[3]), gamma=infl.fac, weights=weights, 
                              data=data, subset=inde1, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights,inde1=inde1)))  
    
    X3    <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    l.sp3 <- length(gam3$sp); if(l.sp3 != 0) sp3 <- gam3$sp
    gp3   <- gam3$nsdf
    y3 <- gam3$y

    n.se1 <- sum(as.numeric(inde1))
    
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
    
    n.se0 <- sum(as.numeric(inde0))

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
    
    n.se1 <- sum(as.numeric(inde1))

    #* this will be useful for prediction, to calculate effects, not required for fitting *#
    
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X3s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
    
  }   
  
  
  
X4 <- X5 <- X6 <- X7 <- X8 <- X9 <- matrix(1, n, 1)  
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
res1 <- residuals(gam1)
res2 <- residuals(gam2)
res3 <- residuals(gam3)

ass.s  <- cor(res1[inde0], res2, method = "kendall")
ass.s  <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))                
i.rho1 <- ass.dp(ass.s, BivD1, scc, sccn, nCa1)

ass.s  <- cor(res1[inde1], res3, method = "kendall")
ass.s  <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))                
i.rho2 <- ass.dp(ass.s, BivD2, scc, sccn, nCa2)

names(i.rho1) <- "theta1.star"
names(i.rho2) <- "theta2.star"

               
##############################################################
# Starting values for whole parameter vector
##############################################################
          

if( margins[2] %in% c(m1d, bl) && margins[3] %in% c(m1d, bl) ) start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, i.rho1, i.rho2)
      

if( margins[2] %in% c(m2, m2d) && margins[3] %in% c(m2, m2d) ){

   start.snR1 <- startsn(margins[2], y2)  
   start.snR2 <- startsn(margins[3], y3)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma1.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma2.star"

start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, i.rho1, i.rho2)
    
}   


if( margins[2] %in% c(m3) && margins[3] %in% c(m3) ){

   start.snR1 <- startsn(margins[2], y2)  
   start.snR2 <- startsn(margins[3], y3)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma1.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma2.star"
   
   log.nu1    <- start.snR1$log.nu.1; names(log.nu1) <- "nu1.star"
   log.nu2    <- start.snR2$log.nu.1; names(log.nu2) <- "nu2.star"
   
start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, log.nu1, log.nu2, i.rho1, i.rho2)
    
}  


if( margins[2] %in% c(m1d) && margins[3] %in% c(m2d) ){

   start.snR2 <- startsn(margins[3], y3)
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma2.star"

start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig2, i.rho1, i.rho2)
    
} 


if( margins[2] %in% c(m2d) && margins[3] %in% c(m1d) ){

   start.snR1 <- startsn(margins[2], y2)  
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma1.star"

start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, i.rho1, i.rho2)
    
}  


if( margins[2] %in% c(m2) && margins[3] %in% c(m3) ){

   start.snR1 <- startsn(margins[2], y2)  
   start.snR2 <- startsn(margins[3], y3)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma1.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma2.star"
   
   log.nu2    <- start.snR2$log.nu.1; names(log.nu2) <- "nu2.star"
   
start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, log.sig1, log.sig2, log.nu2, i.rho1, i.rho2)
    
} 


if( margins[2] %in% c(m3) && margins[3] %in% c(m2) ){

   start.snR1 <- startsn(margins[2], y2)  
   start.snR2 <- startsn(margins[3], y3)
   
   log.sig1   <- start.snR1$log.sig2.1; names(log.sig1) <- "sigma1.star"
   log.sig2   <- start.snR2$log.sig2.1; names(log.sig2) <- "sigma2.star"
   
   log.nu1    <- start.snR1$log.nu.1; names(log.nu1) <- "nu1.star"
   
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
             lsgam3 = lsgam3,
             lsgam4 = lsgam4,
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
             y11.y31R = y11.y31R ) # original n only needed in SemiParBIV.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################

# never thought of it, but subset may cause trouble here? 

if(gamlssfit == TRUE && margins[2] %in% c(m1d, m2d, m2, m3) && margins[3] %in% c(m1d, m2d, m2, m3)){ 

  nstv <- names(start.v) 
  
  form.gamlR <- form.gaml(formula, l.flist, M, type = "ROY")

  gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = inde0,  
                   margin = margins[2], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, drop.unused.levels = drop.unused.levels), list(inde0 = inde0, weights = weights)))   
                   
  gamlss3 <- eval(substitute(gamlss(form.gamlR$formula.gamlss3, data = data, weights = weights, subset = inde1,  
                   margin = margins[2], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, drop.unused.levels = drop.unused.levels), list(inde1 = inde1, weights = weights)))                    
                      
  # updated starting values                   
  
  SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8, sp9 = sp9)
  
  gamls.upsvR <- gamls.upsv(gamlss1 = gamlss2, gamlss2 = gamlss3, margins, M, l.flist, nstv = NULL, VC, GAM, SP, type = "ROY")
  
  
  main.sv                    <- c(gam1$coefficients, gamls.upsvR$start.v) 
  start.v[1:length(main.sv)] <- main.sv
  names(start.v)             <- nstv                                           # this really needed? double check
  
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
          coefficients = SemiParFit$fit$argument, coef.t = NULL, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta1 = SemiParFit.p$theta1, theta2 = SemiParFit.p$theta2, 
          theta1.a = SemiParFit.p$theta1.a, theta2.a = SemiParFit.p$theta2.a,   
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
          etad1 = SemiParFit$fit$etad1, etad2 = SemiParFit$fit$etad2, etas1 = SemiParFit$fit$etas1, etas2 = SemiParFit$fit$etas2,
          etan1 = SemiParFit$fit$etan1, etan2 = SemiParFit$fit$etan2,
          y1 = y1, y2 = y2.m, y3 = y3.m, 
          margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC1 = nC1, nC2 = nC2, hess = hess, 
          respvec = respvec, inde0 = inde0, inde1 = inde1, 
          qu.mag = qu.mag,
          sigma1 = SemiParFit.p$sigma1, sigma1.a = SemiParFit.p$sigma1.a,
          sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a,
          nu1 = SemiParFit.p$nu1, nu1.a = SemiParFit.p$nu1.a, 
          nu2 = SemiParFit.p$nu2, nu2.a = SemiParFit.p$nu2.a, 
          tau1 = SemiParFit.p$tau1, tau1.a = SemiParFit.p$tau1.a, 
          tau2 = SemiParFit.p$tau2, tau2.a = SemiParFit.p$tau2.a,
          Vb.t = SemiParFit.p$Vb.t,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, gp9 = gp9, 
          X2s = X2s, X3s = X3s, X4s = X4s, X5s = X5s, X6s = X6s, X7s = X7s, X8s = X8s, X9s = X9s, 
          VC = VC, Model = Model, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "NO",  
          l.flist = l.flist, v1 = v1, v2 = v2, v3 = v3, triv = FALSE, univar.gamlss = FALSE,
          gamlss2 = gamlss2, gamlss3 = gamlss3, BivD1 = BivD1, BivD2 = BivD2, dof1 = dof1, dof2 = dof2, 
          dof1.a = dof1, dof2.a = dof2, call = cl,
          surv = FALSE, surv.flex = surv.flex)

class(L) <- c("SemiParBIV", "gjrm", "SemiParROY") # check that this is ok

L

}

