SemiParTRIV <- function(formula, data = list(), weights = NULL, subset = NULL,
                             Model = "T", margins = c("probit","probit","probit"), 
                             penCor = "unpen", sp.penCor = 3, approx = FALSE, Chol = FALSE, 
                             infl.fac = 1, gamma = 1, w.alasso = NULL, 
                             rinit = 1, rmax = 100, 
                             iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t", knots = NULL){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- qu.mag2 <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss <- inde <- spgamlss <- n.sel1 <- n.sel2 <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  ngc <- 2; hess <- TRUE
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- NULL
  opc  <- scc <- sccn <- m2 <- m2d <- m3 <- m3d <- bl <- NULL  
   
  fp <- FALSE
  surv.flex <- FALSE

  X2s <- X3s <- ct <- cta <- nC  <- nCa <- NULL
    
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL  
  sp7 <- gp7 <- gam7 <- X7 <- NULL 
  sp8 <- gp8 <- gam8 <- X8 <- NULL 

  spCor <- NULL

  y1.y2.y3    <- NULL
  y1.y2.cy3   <- NULL
  cy1.y2.y3   <- NULL 
  cy1.y2.cy3  <- NULL
  cy1.cy2.cy3 <- NULL
  cy1.cy2.y3  <- NULL
  y1.cy2.cy3  <- NULL
  y1.cy2.y3   <- NULL

  cy1         <- NULL
  y1.cy2      <- NULL
  y1.y2.cy3   <- NULL
  y1.y2.y3    <- NULL
  
  
  ###################################
  
  if(penCor != "unpen"){ spCor <- sp.penCor; names(spCor) <- "spCor"} 
 
  M <- list(mb = c("T", "TSS", "TESS"), margins = margins, penCor = penCor, w.alasso = w.alasso, 
            extra.regI = extra.regI, Model = Model, Chol = Chol)
  
  if(!is.list(formula)) stop("You must specify a list of equations.")
  l.flist <- length(formula)
  if(l.flist > 3 && l.flist != 6) stop("You have to specify six equations.")

  pream.wm(formula, margins = margins, M, l.flist, type = "triv")
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)        
  pred.varR <- pred.var(formula, l.flist, triv = TRUE) 
   
  v1 <- pred.varR$v1  
  v2 <- pred.varR$v2
  v3 <- pred.varR$v3
  pred.n <- pred.varR$pred.n  
  
              
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$knots <- mf$Chol <- mf$margins <- mf$infl.fac <- mf$rinit <- mf$approx <- mf$gamma <- mf$w.alasso <- mf$rmax <- mf$Model <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$penCor <- mf$sp.penCor <- NULL                           
  mf$drop.unused.levels <- TRUE 
  if(Model=="TSS") mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
  
  
if(Model=="TESS"){
 
     data[is.na(data[, v1[1]]), v1[1]] <- 0
     indS1 <- data[, v1[1]]  
     indS1[is.na(indS1)] <- 0      
     indS1 <- as.logical(indS1)
     
     data[indS1 == FALSE, v2[1]] <- 0
     data[indS1 == FALSE, v3[1]] <- 0  
     
     data <- na.omit(data)
    
                   }

  
  if(Model=="TSS"){
  
     data[is.na(data[, v1[1]]), v1[1]] <- 0
     indS1 <- data[, v1[1]] 
     indS2 <- data[, v2[1]]     
     indS1[is.na(indS1)] <- 0   
     indS2[is.na(indS2)] <- 0       
     indS1 <- as.logical(indS1)
     indS2 <- as.logical(indS2) 
     
     data[indS1 == FALSE, v2[1]] <- 0
     data[indS1 == FALSE, v3[1]] <- 0
     data[indS2 == FALSE, v3[1]] <- 0     
     
     data <- na.omit(data)
     
                   }
  

  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  formula.eq3 <- formula[[3]]  
  
 ##############################################################  
 # Equations 1, 2 and 3
 ##############################################################  
 
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = margins[1]), gamma=infl.fac, weights=weights, data=data, knots = knots),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp
  
  inde1 <- inde2 <- inde2.1 <- rep(TRUE, n) # useful for double ss
  
  ########################
  
  if(Model == "TSS")  inde1 <- inde2 <- as.logical(y1)
  if(Model == "TESS") inde1 <- inde2 <- inde2.1 <- as.logical(y1)
  
  ###########

  gam2 <- eval(substitute(gam(formula.eq2, binomial(link = margins[2]), gamma=infl.fac, weights=weights, data=data, subset=inde1, knots = knots),list(weights=weights,inde1=inde1))) 

  if(Model %in% c("TSS","TESS")){
  
  X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
  if(class(X2s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParTRIV for more information.")    
  
  }
  
  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 
  if(l.sp2 != 0) sp2 <- gam2$sp
  

  ###########
  
  if(Model == "TSS"){ inde2[inde1] <- as.logical(gam2$y); inde2.1 <- inde2[inde1]}
  

  gam3 <- eval(substitute(gam(formula.eq3, binomial(link = margins[3]), gamma=infl.fac, weights=weights, data=data, subset=inde2, knots = knots),list(weights=weights,inde2=inde2))) 


  if(Model %in% c("TSS","TESS")){
  
  X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
  if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParTRIV for more information.")    
    
  }  
    
  X3 <- model.matrix(gam3)
  X3.d2 <- dim(X3)[2]
  l.sp3 <- length(gam3$sp)
  y3 <- gam3$y 
  if(l.sp3 != 0) sp3 <- gam3$sp  


 ##############################################################

if(Model == "T"){

  y1.y2.y3    <- y1*y2*y3 
  y1.y2.cy3   <- y1*y2*(1-y3)
  cy1.y2.y3   <- (1-y1)*y2*y3 
  cy1.y2.cy3  <- (1-y1)*y2*(1-y3) 
  cy1.cy2.cy3 <- (1-y1)*(1-y2)*(1-y3)
  cy1.cy2.y3  <- (1-y1)*(1-y2)*y3
  y1.cy2.cy3  <- y1*(1-y2)*(1-y3)
  y1.cy2.y3   <- y1*(1-y2)*y3
  
} 

if(Model == "TSS"){
 
  cy1       <- (1-y1)
  y1.cy2    <- y1[inde1]*(1-y2)
  y1.y2.cy3 <- y1[inde2]*y2[inde2.1]*(1-y3) 
  y1.y2.y3  <- y1[inde2]*y2[inde2.1]*y3  
  n.sel1 <- table(inde1)[2]
  n.sel2 <- table(inde2)[2]

}

  if(Model == "TESS"){
     
    inde <- inde1
    
    cy1        <- (1-y1)
    y1.y2.y3   <- y1[inde]*y2*y3  
    y1.y2.cy3  <- y1[inde]*y2*(1-y3)  
    y1.cy2.cy3 <- y1[inde]*(1-y2)*(1-y3)
    y1.cy2.y3  <- y1[inde]*(1-y2)*y3 
    
  }


  gp1 <- gam1$nsdf 
  gp2 <- gam2$nsdf
  gp3 <- gam3$nsdf

##############################################################
# Starting values
##############################################################

if(Model == "T"){

#res1 <- residuals(gam1)
#res2 <- residuals(gam2)
#res3 <- residuals(gam3)

#cor1 <- cor(res1, res2)
#cor2 <- cor(res1, res3)
#cor3 <- cor(res2, res3)

tcorrs <- tetrachoric(data[,c(v1[1],v2[1],v3[1])])$rho

cor1 <- tcorrs[1,2]
cor2 <- tcorrs[1,3]
cor3 <- tcorrs[2,3]

cor1 <- sign(cor1)*ifelse(abs(cor1) > 0.85, 0.85, abs(cor1))
cor2 <- sign(cor2)*ifelse(abs(cor2) > 0.85, 0.85, abs(cor2))
cor3 <- sign(cor3)*ifelse(abs(cor3) > 0.85, 0.85, abs(cor3))

}

if(Model == "TSS" || Model == "TESS"){ cor1 <- cor2 <- cor3 <- 0.01; theta12 <- theta13 <- theta23 <- atanh(cor1) } # this can be improved for TESS


if(Chol == FALSE){
      
      theta12 <- atanh(cor1)
      theta13 <- atanh(cor2)
      theta23 <- atanh(cor3)
      
                 }
    
  if(Chol == TRUE) {
    SigmaSt <- matrix(c(1, cor1, cor2, cor1, 1, cor3, cor2, cor3, 1), 3, 3)
    SigmaSt <- PosDefCor(SigmaSt)
    theta12 <- SigmaSt[1, 2]
    theta13 <- SigmaSt[1, 3]
    theta23 <- SigmaSt[2, 3]
    
    theta12.st <- sign(theta12) * sqrt( theta12^2/(1 - theta12^2) )
    th23sol    <- ( (theta23 * sqrt(1 + theta12.st^2) - theta12.st * theta13)/sqrt(1 - theta13^2) )^2
    theta13.st <- sign(theta13) * sqrt( (theta13^2 * (1 + th23sol/(1 - th23sol)))/(1 - theta13^2) )
    theta23.st <- sign(theta23) * sqrt( th23sol/(1 - th23sol) )
    
    theta12 <- theta12.st
    theta13 <- theta13.st
    theta23 <- theta23.st
  }


names(theta12) <- "theta12.st"
names(theta13) <- "theta13.st"
names(theta23) <- "theta23.st"
  
start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, theta12, theta13, theta23)


##############################################################  
# starting values for case of predictors on all parameters
##############################################################  
  
    if(l.flist > 3){
    
    vo <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, theta12 = theta12, theta13 = theta13, theta23 = theta23, n = n )

    overall.svGR <- overall.svG(formula = formula, data = data, ngc = 2, margins = margins, M = M, vo = vo, gam1 = gam1, gam2 = gam2, gam3 = gam3, type = "triv", knots = knots)
        
    start.v = overall.svGR$start.v 
    X4 = overall.svGR$X4; X5 = overall.svGR$X5
    X6 = overall.svGR$X6; X7 = overall.svGR$X7; X8 = overall.svGR$X8
    X4.d2 = overall.svGR$X4.d2; X5.d2 = overall.svGR$X5.d2
    X6.d2 = overall.svGR$X6.d2; X7.d2 = overall.svGR$X7.d2; X8.d2 = overall.svGR$X8.d2
    gp4 = overall.svGR$gp4; gp5 = overall.svGR$gp5
    gp6 = overall.svGR$gp6; gp7 = overall.svGR$gp7; gp8 = overall.svGR$gp8
    gam4 = overall.svGR$gam4; gam5 = overall.svGR$gam5
    gam6 = overall.svGR$gam6; gam7 = overall.svGR$gam7; gam8 = overall.svGR$gam8
    l.sp4 = overall.svGR$l.sp4; l.sp5 = overall.svGR$l.sp5
    l.sp6 = overall.svGR$l.sp6; l.sp7 = overall.svGR$l.sp7; l.sp8 = overall.svGR$l.sp8
    sp4 = overall.svGR$sp4; sp5 = overall.svGR$sp5
    sp6 = overall.svGR$sp6; sp7 = overall.svGR$sp7; sp8 = overall.svGR$sp8
    
    }


##############################################################

if(Model == "T")    func.opt <- triprobgHs
if(Model == "TSS")  func.opt <- triprobgHsSS
if(Model == "TESS") func.opt <- triprobgHsESS

##########################################################
# SPs and penalties
##########################################################
  

GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8)   


if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0) && fp==FALSE ){ 

L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), l.gam4 = length(gam4$coefficients),
              l.gam5 = length(gam5$coefficients), l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), l.gam8 = length(gam8$coefficients))

L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8)

                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
                 qu.mag <- S.m(GAM, L.SP, L.GAM)                             
                                                        } 
 
  
if(!(penCor %in% c("unpen")) && fp==FALSE){
    
    l.sp4 <- 1
    sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, spCor) # all this part is fine when using Chol = T since 
                                                      # we dont' allow for penalty and additive predictors on corrs
  
########################################### 
# no penalty for corrs they are added after  
########################################### 
 
# here I set up rank and off as I can not do it later
# then I will incorporate the penalty for the corrs
 
if(l.sp1 == 0 && l.sp2 == 0 && l.sp3 == 0) qu.mag <- list(rank = 3, off = X1.d2 + X2.d2 + X3.d2 + 1, Ss = NULL) 
if(l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0) qu.mag <- list(rank = c(qu.mag$rank, 3), off = c(qu.mag$off, X1.d2 + X2.d2 + X3.d2 + 1), Ss = qu.mag$Ss)

}


##########################################################


if(missing(parscale)) parscale <- 1   

  respvec <- list(y1 = y1,
                  y2 = y2,
                  y3 = y3,
  		  y1.y2.y3    = y1.y2.y3   , 
  		  y1.y2.cy3   = y1.y2.cy3  ,
  		  cy1.y2.y3   = cy1.y2.y3  , 
  		  cy1.y2.cy3  = cy1.y2.cy3 , 
       		  cy1.cy2.cy3 = cy1.cy2.cy3,
  		  cy1.cy2.y3  = cy1.cy2.y3 , 
                  y1.cy2.cy3  = y1.cy2.cy3 , 
                  y1.cy2.y3   = y1.cy2.y3,                    
                  cy1       = cy1,
                  y1.cy2    = y1.cy2,
                  y1.y2.cy3 = y1.y2.cy3, univ = 0)
 
  lsgam1 <- length(gam1$smooth)
  lsgam2 <- length(gam2$smooth)
  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)
   
  VC <- list(lsgam1 = lsgam1,
             lsgam2 = lsgam2,
             lsgam3 = lsgam3,
             lsgam4 = lsgam4,
             lsgam5 = lsgam5,
             lsgam6 = lsgam6,
             lsgam7 = lsgam7,
             lsgam8 = lsgam8, 
             X1 = X1, inde = inde, inde1 = inde1, inde2 = inde2, inde2.1 = inde2.1,
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,
             X7 = X7,
             X8 = X8,
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,   
             X7.d2 = X7.d2, 
             X8.d2 = X8.d2,
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6,  
             gp7 = gp7,
             gp8 = gp8, 
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6,  
             l.sp7 = l.sp7,
             l.sp8 = l.sp8,
             infl.fac = infl.fac,
             weights = weights, univ.gamls = FALSE,
             fp = fp,
             hess = hess, nCa = nCa,
             Model = Model, gamlssfit = TRUE,
             end = NULL,
             BivD = "N", penCor = penCor,
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "NO", ccss = "no", m2 = m2, m3 = m3, m2d = m2d, m3d = m3d, bl = bl, triv = TRUE,
             X2s = X2s, X3s = X3s,
             approx = approx, gamma = gamma, wc = w.alasso, qu.mag = qu.mag,
             zerov = -10, Chol = Chol, surv.flex = surv.flex, l.flist = l.flist)
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################

  SemiParFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = VC$qu.mag) # there is a reason for calling qu.mag from VC
                                            
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParTRIV.fit.post(SemiParFit = SemiParFit, VC = VC, Model = Model, GAM)
 
  SemiParFit <- SemiParFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.    
    
  ##########################################################################################################################
                                            
                                            
if(gc.l == TRUE) gc()

e.v <- round(min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values), 6)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?SemiParTRIProbit."

if(gradi > 10 && e.v < 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 10 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 10 && e.v < 0)  warning(paste(me2,"\n",me3), call. = FALSE)
  ##########################################################################################################################

gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  # for all.terms
  ##########################################################################################################################

L <- list(fit = SemiParFit$fit, formula = formula, Model = Model,
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8,
          coefficients = SemiParFit$fit$argument, coef.t = NULL, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta12   = SemiParFit.p$theta12, 
          theta12.a = SemiParFit.p$theta12.a,
          theta13   = SemiParFit.p$theta13, 
          theta13.a = SemiParFit.p$theta13.a,
          theta23   = SemiParFit.p$theta23, 
          theta23.a = SemiParFit.p$theta23.a,           
          n = n, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,            
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf11 = SemiParFit.p$edf11,   # what is this for?
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6,
          edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, Vb.t = SemiParFit.p$Vb.t,
          wor.c = SemiParFit$wor.c,
          p111 = SemiParFit$fit$p111,
          p011 = SemiParFit$fit$p011,
          p101 = SemiParFit$fit$p101,
          p110 = SemiParFit$fit$p110,
          p100 = SemiParFit$fit$p100,
          p010 = SemiParFit$fit$p010,
          p001 = SemiParFit$fit$p001,
          p000 = SemiParFit$fit$p000,
          eta1 = SemiParFit$fit$eta1, 
          eta2 = SemiParFit$fit$eta2,
          eta3 = SemiParFit$fit$eta3,
          y1 = y1, y2 = y2, y3 = y3,  n.sel1 = n.sel1, n.sel2 = n.sel2, 
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = hess, 
          respvec = respvec, 
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8,
          VC = VC, magpp = SemiParFit$magpp,
          gamlss = gamlss, gamlssfit = TRUE, Cont = "NO", triv = TRUE,  
          l.flist = l.flist, margins = margins,
          inde2 = inde2, X2s = X2s, X3s = X3s,
          p1n = SemiParFit.p$p1n, p2n = SemiParFit.p$p2n, p3n = SemiParFit.p$p3n, v1 = v1, v2 = v2, v3 = v3, univar.gamlss = FALSE, call = cl,
          surv = FALSE, surv.flex = surv.flex)

class(L) <- c("SemiParTRIV","SemiParBIV")

L

}

