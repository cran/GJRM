copulaSampleSel <- function(formula, data = list(), weights = NULL, subset = NULL,
                             BivD = "N", margins = c("probit","N"), dof = 3,
                             fp = FALSE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t", knots = NULL,
                             drop.unused.levels = TRUE,
                             min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
    BivD2 <- c("C0C90","C0C270","C180C90","C180C270",
               "J0J90","J0J270","J180J90","J180J270",
               "G0G90","G0G270","G180G90","G180G270")
  
  # stop("check next release")
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- qu.mag2 <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- y2m <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  ngc <- 2
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- NULL
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL  
  sp7 <- gp7 <- gam7 <- X7 <- NULL
  sp8 <- gp8 <- gam8 <- X8 <- NULL     
  dataset <- gamlss <- NULL
  
  Sl.sf <- NULL
  sp.method <- "perf"
  
  X2s <- X3s <- X4s <- X5s <- X6s <- X7s <- X8s <- NULL 
    surv.flex <- FALSE

    BivD2 <- c("C0C90","C0C270","C180C90","C180C270",
               "J0J90","J0J270","J180J90","J180J270",
               "G0G90","G0G270","G180G90","G180G270",
               "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270")
  
  
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM","T","PL","HO","GAL0", "GAL90", "GAL180", "GAL270")
  scc  <- c("C0", "C180", "GAL0" , "GAL180", "J0", "J180", "G0", "G180", BivD2)
  sccn <- c("C90", "C270", "GAL90", "GAL270", "J90", "J270", "G90", "G270")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK","GP","GPII","GPo")
  m3   <- c("DAGUM","SM","TW")
  m1d  <- c("PO", "ZTP")
  m2d  <- c("NBI","NBII","PIG","DGP", "DGPII")
  m3d  <- NULL
  bl   <- c("probit", "logit", "cloglog")

  M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, BivD = BivD, 
               opc = opc, extra.regI = extra.regI, margins = margins, bl = bl, dof = dof)

  M$K1 <- NULL

  ct  <- data.frame( c(opc), c(1:14,55,56,57,60,61,62:65) )
  cta <- data.frame( c(opc), c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2,60,61,62:65) )                   




  
  if(BivD %in% BivD2){
  
  if(BivD %in% BivD2[1:4])  BivDt <- "C0" 
  if(BivD %in% BivD2[5:12]) BivDt <- "J0"
  if(BivD %in% BivD2[13 :16]) BivDt <- "C0" # useful for ass dep function but we calculate it differently, so ok like this

  
  nC  <-  ct[which( ct[,1]==BivDt),2]
  nCa <- cta[which(cta[,1]==BivDt),2]     
  
  }
  
  
  if(!(BivD %in% BivD2)){
    
  nC  <-  ct[which( ct[,1]==BivD),2]
  nCa <- cta[which(cta[,1]==BivD),2]     
    
  }
  
  





 #######################################################################################
 if(!is.list(formula)) stop("You must specify a list of equations.")
 M$l.flist <- l.flist <- length(formula) 
 pream.wm(formula, margins, M, l.flist, type = "copSS")
 form.check(formula, l.flist)
 
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
            
  pred.varR <- pred.var(formula, l.flist) 
   
  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  
  
  
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$ordinal <- mf$knots <- mf$BivD <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- NULL                           
  mf$drop.unused.levels <- drop.unused.levels 
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()
        
     data[is.na(data[, v1[1]]), v1[1]] <- 0   
     indS <- data[, v1[1]]    
     indS[is.na(indS)] <- 0   
     indS <- as.logical(indS)  
     data[indS == FALSE, v2[1]] <- 0.01  
     data <- na.omit(data)

  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
    
 ##############################################################  
 # Equation 1
 ##############################################################  
   
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = margins[1]), gamma=infl.fac, weights=weights, 
                              data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights))) 

  X1 <- model.matrix(gam1) ## this may have to be changed to make it more generic using predict maybe
                           ## also in other functions in the package but not really needed maybe
                           ## it depends on the purpose
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp
  gp1 <- gam1$nsdf 
 
 ##############################################################  
 # Equation 2
 ##############################################################   
 
 inde <- as.logical(y1) 
 n.sel <- sum(as.numeric(inde)) 

 form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], m1d, m2d, copSS = TRUE, inde = inde)  
 
 formula.eq2  <- form.eq12R$formula.eq1
 formula.eq2r <- form.eq12R$formula.eq1r
 y2           <- form.eq12R$y1
 y2.test      <- form.eq12R$y1.test 
 y2m          <- form.eq12R$y1m
 
 gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data, subset=inde, knots = knots, drop.unused.levels = drop.unused.levels),list(weights = weights, inde = inde)))
    
    # gam2 not entirely efficient given that we will do gam2.1 but we keep for now
    ######
    # TEST
    ######
    X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(any(class(X2s)=="try-error")) stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######
    
    gam2$formula <- formula.eq2r  
    names(gam2$model)[1] <- as.character(formula.eq2r[2])
    X2.d2 <- length(gam2$coefficients)    
    X2 <- model.matrix(gam2) 
    l.sp2 <- length(gam2$sp)
    if(l.sp2 != 0) sp2 <- gam2$sp     
    gp2 <- gam2$nsdf   
  
    cy <- 1 - y1  
    
    y2 <- y2.test 
    # y2 <- y2[inde] 
    if( margins[2] %in% c("LN") ) y2 <- log(y2) 
  
#######################################################################
# Starting values for dependence parameter and main coefficients
#######################################################################

  form.eq2imr <- update.formula(formula.eq2, ~. + imrGUANN) 
  p.g1 <- predict.gam(gam1)
  imrGUANN <- data$imrGUANN <- dnorm(p.g1)/pnorm(p.g1)
    
  gam2.1 <- eval(substitute(gam(form.eq2imr, gamma=infl.fac, weights=weights, data=data, subset=inde, knots = knots, drop.unused.levels = drop.unused.levels),list(weights = weights, inde = inde)))
  pimr   <- which(names(gam2.1$coefficients)=="imrGUANN")
  c.gam2 <- gam2.1$coefficients[-pimr]
  
  if(l.sp2 != 0) sp2 <- gam2.1$sp
  
  sia <- sqrt(mean(residuals(gam2.1, type = "deviance")^2)+mean(imrGUANN[inde]*(imrGUANN[inde]+p.g1[inde]))*gam2.1$coefficients["imrGUANN"]^2)[[1]]
  co  <- (gam2.1$coefficients["imrGUANN"]/sia)[[1]] 
  
  ass.s <- sign(co)*ifelse(abs(co) > 0.5, 0.5, abs(co))

  i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)
             
##############################################################
# Starting values for whole parameter vector
##############################################################
 
log.nu.2 <- log.sig2.2 <- NULL

if( !(margins[2] %in% c(m1d)) ){

start.snR <- startsn(margins[2], y2)
    
log.sig2.2 <- start.snR$log.sig2.1; names(log.sig2.2) <- "sigma.star"
if( margins[2] %in% c(m3) ){ log.nu.2 <- start.snR$log.nu.1; names(log.nu.2)   <- "nu.star"}     

}

vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2.2 = log.sig2.2, 
           log.nu.2 = log.nu.2, log.nu.1 = NULL, log.sig2.1 = NULL, dof.st = NULL, n = n, drop.unused.levels = drop.unused.levels)
start.v <- overall.sv(margins, M, vo, type = "copSS", c.gam2 = c.gam2)
 
 
##############################################################  
# starting values for case of predictors on all parameters
##############################################################  
  
    if(l.flist > 2){
    
    overall.svGR <- overall.svG(formula, data, ngc = 2, margins, M, vo, gam1, gam2, type = "copSS", inde = inde, c.gam2 = c.gam2, knots = knots)
    
    start.v = overall.svGR$start.v 
    X3 = overall.svGR$X3; X4 = overall.svGR$X4; X5 = overall.svGR$X5
    X3s = overall.svGR$X3s; X4s = overall.svGR$X4s; X5s = overall.svGR$X5s
    X6 = overall.svGR$X6; X7 = overall.svGR$X7; X8 = overall.svGR$X8
    X3.d2 = overall.svGR$X3.d2; X4.d2 = overall.svGR$X4.d2; X5.d2 = overall.svGR$X5.d2
    X6.d2 = overall.svGR$X6.d2; X7.d2 = overall.svGR$X7.d2; X8.d2 = overall.svGR$X8.d2
    gp3 = overall.svGR$gp3; gp4 = overall.svGR$gp4; gp5 = overall.svGR$gp5
    gp6 = overall.svGR$gp6; gp7 = overall.svGR$gp7; gp8 = overall.svGR$gp8
    gam3 = overall.svGR$gam3; gam4 = overall.svGR$gam4; gam5 = overall.svGR$gam5
    gam6 = overall.svGR$gam6; gam7 = overall.svGR$gam7; gam8 = overall.svGR$gam8
    l.sp3 = overall.svGR$l.sp3; l.sp4 = overall.svGR$l.sp4; l.sp5 = overall.svGR$l.sp5
    l.sp6 = overall.svGR$l.sp6; l.sp7 = overall.svGR$l.sp7; l.sp8 = overall.svGR$l.sp8
    sp3 = overall.svGR$sp3; sp4 = overall.svGR$sp4; sp5 = overall.svGR$sp5
    sp6 = overall.svGR$sp6; sp7 = overall.svGR$sp7; sp8 = overall.svGR$sp8
    
    }

##########################################################
  
GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, K1 = NULL)   


if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0) && fp==FALSE ){ 

L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), l.gam4 = length(gam4$coefficients),
              l.gam5 = length(gam5$coefficients), l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), l.gam8 = length(gam8$coefficients))

L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8)

                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
                 qu.mag <- S.m(GAM, L.SP, L.GAM)                             
                                                        }     

  
  
##########################################################


if(missing(parscale)) parscale <- 1   


  respvec <- list(y1 = y1,
                  y2 = y2,
                  y1.y2 = y1.y2, 
                  y1.cy2 = y1.cy2, 
                  cy1.y2 = cy1.y2, 
                  cy1.cy2 = cy1.cy2, 
                  cy1 = cy1,
                  cy = cy, univ = 0)

  lsgam1 <- length(gam1$smooth)
  lsgam2 <- length(gam2$smooth)
  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)

  VC <- list(lsgam1 = lsgam1, robust = FALSE, K1 = NULL,
             lsgam2 = lsgam2, Sl.sf = Sl.sf, sp.method = sp.method,
             lsgam3 = lsgam3, 
             lsgam4 = lsgam4,
             lsgam5 = lsgam5,
             lsgam6 = lsgam6,
             lsgam7 = lsgam7,
             lsgam8 = lsgam8,
             X1 = X1, 
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
             nCa = nCa,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6,    
             l.sp7 = l.sp7, 
             l.sp8 = l.sp8, 
             infl.fac = infl.fac,
             weights = weights,
             fp = fp,
             hess = NULL, univ.gamls = FALSE,
             Model = "B", 
             end = end,
             BivD = BivD, dof.st = log(dof - 2), dof = dof,
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "NO", ccss = "yes", m2 = m2, m3 = m3, 
             m1d = m1d, m2d = m2d, m3d = m3d, bl = bl, inde = inde,
             X2s = X2s, X3s = X3s, X4s = X4s, X5s = X5s, triv = FALSE, y2m = y2m, zerov = -10,
             BivD2 = BivD2, cta = cta, ct = ct, surv.flex = surv.flex, gp2.inf = NULL,
             informative = "no", sp.fixed = NULL,
             zero.tol = 1e-02,
             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr) 
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################
  signind <- 1

  func.opt <- func.OPT(margins, M, type = "copSS") 

  SemiParFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
   
    if(BivD %in% BivD2) rm(signind, envir = .GlobalEnv)
 
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- copulaSampleSel.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)                       
  SemiParFit   <- SemiParFit.p$SemiParFit 
 
  y2.m <- y2  
  if(margins[2] == "LN") y2.m <- exp(y2) # y2.m[inde] <- exp(y2[inde])

  if(gc.l == TRUE) gc()

  ##########################################################################################################################


cov.c(SemiParFit)

  ##########################################################################################################################
gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  # for all.terms
  ##########################################################################################################################


# this is for mice # 

formula.aux <- formula

formula.aux[[1]] <- reformulate( attr(terms(formula.aux[[1]]), "term.labels"), response = "ry" )
environment(formula.aux[[1]]) <- environment(formula[[1]])


formula.aux[[2]] <- reformulate( attr(terms(formula.aux[[2]]), "term.labels"), response = "y" )
environment(formula.aux[[2]]) <- environment(formula[[2]])

###


L <- list(fit = SemiParFit$fit, dataset = dataset, formula = formula, mice.formula = formula.aux, robust = FALSE,
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, 
          gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, coef.t = NULL, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,Vb.t = SemiParFit.p$Vb.t,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,
          OR = SemiParFit.p$OR, 
          GM = SemiParFit.p$GM,    
          n = n, n.sel = n.sel, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
          X7.d2 = X7.d2, X8.d2 = X8.d2,             
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf11 = SemiParFit.p$edf11,   
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
          edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p11 = SemiParFit$fit$p11, p10 = SemiParFit$fit$p10, p01 = SemiParFit$fit$p01, p00 = SemiParFit$fit$p00, 
          p1 = SemiParFit$fit$p1, p2 = SemiParFit$fit$p2,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, etad = SemiParFit$fit$etad,
          etas = SemiParFit$fit$etas, etan = SemiParFit$fit$etan,
          y1 = y1, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = TRUE, 
          respvec = respvec, inde = inde, 
          qu.mag = qu.mag, 
          sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a,
          sigma = SemiParFit.p$sigma2, sigma.a = SemiParFit.p$sigma2.a,
          nu = SemiParFit.p$nu,     nu.a = SemiParFit.p$nu.a,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          X2s = X2s, 
          X3s = X3s,
          X4s = X4s,
          X5s = X5s,
          p1n=SemiParFit.p$p1n , p2n = SemiParFit.p$p2n, 
          VC = VC, Model = NULL, magpp = SemiParFit$magpp,
          gamlss = gamlss, Cont = "NO", tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, 
          l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE, BivD2 = BivD2, dof = dof, dof.a = dof, call = cl,
          surv = FALSE, surv.flex = surv.flex)

class(L) <- c("copulaSampleSel", "SemiParBIV","gjrm")

L

}

