SemiParBIV <- function(formula, data = list(), weights = NULL, subset = NULL,
                             Model = "B", BivD = "N", margins = c("probit","probit"), 
                             dof = 3, gamlssfit = FALSE,
                             fp = FALSE, hess = TRUE, infl.fac = 1, theta.fx = NULL, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t", intf = FALSE, knots = NULL,
                             drop.unused.levels = TRUE){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- y2m <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- i.rho <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- gamlss2 <- dof.st <- NULL
  gamlss2 <- NULL
  sp1 <- sp2 <- c.gam2 <-  X2s <- X3s <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL  
  sp7 <- gp7 <- gam7 <- X7 <- NULL  
  sp8 <- gp8 <- gam8 <- X8 <- NULL  
  log.sig2 <- log.nu <- NULL
  Sl.sf <- NULL
  sp.method <- "perf"
    
  BivD2 <- c("C0C90","C0C270","C180C90","C180C270",
             "J0J90","J0J270","J180J90","J180J270",
             "G0G90","G0G270","G180G90","G180G270")  
  
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM","T","PL","HO")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180", BivD2)
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  mb   <- c("B", "BSS", "BPO", "BPO0")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK","GP","GPII","GPo")
  m3   <- c("DAGUM","SM","TW")
  m1d  <- c("PO", "ZTP") 
  m2d  <- c("NBI", "NBII", "PIG","DGP","DGPII") 
  bl   <- c("probit", "logit", "cloglog") # , "cauchit")   
  M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, BivD = BivD, 
               opc = opc, extra.regI = extra.regI, margins = margins, bl = bl, intf = intf,
               theta.fx = theta.fx, Model = Model, mb = mb, BivD2 = BivD2, dof = dof)  
  surv.flex <- FALSE

  M$K1 <- NULL
  ct  <- data.frame( c(opc),
                    c(1:14,55,56,57,60,61) 
                     )
  cta <- data.frame( c(opc),
                     c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2,60,61) 
                     )                   
if(BivD %in% BivD2){
  
  if(BivD %in% BivD2[1:4])  BivDt <- "C0" 
  if(BivD %in% BivD2[5:12]) BivDt <- "J0"
  
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
  pream.wm(formula, margins, M, l.flist, type = "biv")
  form.check(formula, l.flist)   
  
 #######################################################################################    
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
            
  pred.varR <- pred.var(formula, l.flist) 
   
  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  
  
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$ordinal <- mf$knots <- mf$dof <- mf$intf <- mf$theta.fx <- mf$Model <- mf$BivD <- mf$margins <- mf$fp <- mf$hess <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  mf$drop.unused.levels <- drop.unused.levels
  
  if(Model=="BSS") mf$na.action <- na.pass
  
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
    if(gc.l == TRUE) gc()  

   if(Model=="BSS"){     
   
     data[is.na(data[, v1[1]]), v1[1]] <- 0
     indS <- data[, v1[1]]    
     indS[is.na(indS)] <- 0   
     indS <- as.logical(indS)  
     data[indS == FALSE, v2[1]] <- 0  
     data <- na.omit(data)   
   
                   }
                                      
 
  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  
  
  if(Model=="B"){
    if(v1[1] %in% v2[-1]) end <- 1
    if(v2[1] %in% v1[-1]) end <- 2
                }


 ##############################################################  
 # Equation 1
 ##############################################################  
   
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = margins[1]), gamma=infl.fac, weights=weights, 
                              data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp
  gp1 <- gam1$nsdf 
    
  inde <- rep(TRUE, n) 

 ##############################################################
 # Equation 2 for BPO and binary B
 ##############################################################  

  if((Model=="B" || Model=="BPO" || Model=="BPO0") && margins[2] %in% bl){
  
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link=margins[2]), gamma=infl.fac, weights=weights, 
                           data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights))) # check at later stage the need of eval and substitute

  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 
  if(l.sp2 != 0) sp2 <- gam2$sp
     
  if(Model=="B"){  
  	y1.y2 <- y1*y2
  	y1.cy2 <- y1*(1-y2)
  	cy1.y2 <- (1-y1)*y2
  	cy1.cy2 <- (1-y1)*(1-y2)
                }

  if(Model=="BPO" || Model=="BPO0" ) cy <- 1 - y1
  
  } 
  
 ##############################################################
 # Equation 2 for B and continuous/discrete response 
 ##############################################################  

  if(Model=="B" && !(margins[2] %in% bl) ){
  
    form.eq12R   <- form.eq12(formula.eq2, data, v2, margins[2], m1d, m2d)   
    formula.eq2  <- form.eq12R$formula.eq1
    formula.eq2r <- form.eq12R$formula.eq1r
    y2           <- form.eq12R$y1
    y2.test      <- form.eq12R$y1.test 
    y2m          <- form.eq12R$y1m  
  
    gam2         <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))
    gam2$formula <- formula.eq2r  
    names(gam2$model)[1] <- as.character(formula.eq2r[2])

    y2 <- y2.test  
    if( margins[2] %in% c("LN") ) y2 <- log(y2) 

    X2 <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp)
    if(l.sp2 != 0) sp2 <- gam2$sp
    
    cy <- 1 - y1
    
  } 

 ##############################################################
 # Equation 2 for BSS 
 ##############################################################  

  if(Model=="BSS"){

  inde <- as.logical(y1)

  gam2 <- eval(substitute(gam(formula.eq2, binomial(link = margins[2]), gamma=infl.fac, weights=weights, 
                              data=data, subset=inde, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights,inde=inde)))  
  
######
# TEST
######
X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
if(class(X2s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParBIV for more information.") 
######  
  
  X2.d2 <- length(gam2$coefficients)
  X2 <- model.matrix(gam2) 
  y2 <- gam2$y
  n.sel <- sum(as.numeric(inde))
  
  l.sp2 <- length(gam2$sp)
  if(l.sp2 != 0) sp2 <- gam2$sp 
  
  cy1 <- (1-y1)
  y1.y2 <- y1[inde]*y2
  y1.cy2 <- y1[inde]*(1-y2)
  
  ######################
  # better starting values in BSS case
  
  form.eq2imr <- update.formula(formula.eq2, ~. + imrGUANN) 
  p.g1 <- predict.gam(gam1)
  imrGUANN <- data$imrGUANN <- dnorm(p.g1)/pnorm(p.g1)
    
  gam2.1 <- eval(substitute(gam(form.eq2imr, gamma=infl.fac, binomial(link = margins[2]), weights=weights, data=data, subset=inde, knots = knots, drop.unused.levels = drop.unused.levels),list(weights = weights, inde = inde)))
  pimr   <- which(names(gam2.1$coefficients)=="imrGUANN")
  
  c.gam2 <- gam2.1$coefficients[-pimr]
  
  sia <- sqrt(mean(residuals(gam2.1, type = "deviance")^2)+mean(imrGUANN[inde]*(imrGUANN[inde]+p.g1[inde]))*gam2.1$coefficients["imrGUANN"]^2)[[1]]
  co  <- (gam2.1$coefficients["imrGUANN"]/sia)[[1]] 
  
  ass.s <- co 
  ass.s <- sign(ass.s)*ifelse(abs(ass.s) > 0.2, 0.2, abs(ass.s))
  
  if(l.sp2 != 0) sp2 <- gam2.1$sp
  
  }
  
##############################################################  

  gp2 <- gam2$nsdf
  
##############################################################
# Starting values for dependence parameter
##############################################################

   
if(is.null(theta.fx)){

 
if( !(Model %in% c("BPO0")) ){    
    
if(Model=="B"){ 

  res1 <- residuals(gam1)
  res2 <- residuals(gam2)
  ass.s <- cor(res1, res2, method = "kendall")
  ass.s <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))  
  
              }  

if(Model=="BPO") ass.s <- 0.01
                          
i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)

}

if(Model=="BPO0")  i.rho <- 0

}

names(i.rho) <- "theta.star"  

               
##############################################################
# Starting values for whole parameter vector
##############################################################
           
         
           
if(margins[1] %in% bl && margins[2] %in% c(bl,m1d) && Model %in% c("B","BPO") && is.null(theta.fx)) start.v <- c(gam1$coefficients, gam2$coefficients, i.rho) 

if(Model == "BSS")  start.v <- c(gam1$coefficients, c.gam2, i.rho) 
if(Model == "BPO0") start.v <- c(gam1$coefficients, gam2$coefficients)      
if(!is.null(theta.fx) && Model == "B" && margins[2] %in% bl ) start.v <- c(gam1$coefficients, gam2$coefficients) 

if(margins[1] %in% bl && margins[2] %in% c(m2,m3,m2d) ){

   start.snR <- startsn(margins[2], y2)    
   log.sig2  <- start.snR$log.sig2.1; names(log.sig2) <- "sigma.star"
   if( margins[2] %in% c(m3) ){ log.nu <- start.snR$log.nu.1; names(log.nu) <- "nu.star"}       
     
   if(margins[2] %in% c(m2,m2d)) start.v <- c(gam1$coefficients, gam2$coefficients, log.sig2,         i.rho)        
   if(margins[2] %in% m3)        start.v <- c(gam1$coefficients, gam2$coefficients, log.sig2, log.nu, i.rho)                                  

                                                       } 
    
  
 
##############################################################  
  
if(l.flist > 2){  
 
vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2 = log.sig2, log.nu = log.nu, n = n, drop.unused.levels = drop.unused.levels)  
  
    overall.svGR <- overall.svG(formula, data, ngc = 2, margins, M, vo, gam1, gam2, type = "biv", inde = inde, c.gam2 = c.gam2, knots = knots)
    
    start.v = overall.svGR$start.v 
    X3 = overall.svGR$X3; X4 = overall.svGR$X4; X5 = overall.svGR$X5
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
    X3s = overall.svGR$X3s; X4s = overall.svGR$X4s
  
  
  
}  
  
  
  
  
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

  VC <- list(lsgam1 = lsgam1, robust = FALSE, sp.fixed = NULL,K1 = NULL,
             lsgam2 = lsgam2, Sl.sf = Sl.sf, sp.method = sp.method,
             lsgam3 = lsgam3,
             lsgam4 = lsgam4,
             lsgam5 = lsgam5,
             lsgam6 = lsgam6,
             lsgam7 = lsgam7,
             lsgam8 = lsgam8, 
             X1 = X1, inde = inde,my.env=my.env,
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
             weights = weights,
             fp = fp, univ.gamls = FALSE,
             hess = hess, nCa = nCa,
             Model = Model, gamlssfit = gamlssfit,
             end = end,
             BivD = BivD, dof.st = log(dof - 2), dof = dof, 
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "NO", ccss = "no", m2 = m2, m3 = m3, m2d = m2d, m1d = m1d, m3d = NULL, bl = bl,
             X2s = X2s, X3s = X3s, triv = FALSE, y2m = y2m,
             theta.fx = theta.fx, i.rho = i.rho, 
             BivD2 = BivD2, cta = cta, ct = ct, zerov = -10, surv.flex = surv.flex, gp2.inf = NULL,
             informative = "no") # original n only needed in SemiParBIV.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################



if(gamlssfit == TRUE && !(margins[2] %in% bl)){ 

  form.gamlR <- form.gaml(formula, l.flist, M, type = "biv")

  gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = subset,  
                   margin = margins[2], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, drop.unused.levels = drop.unused.levels), list(weights=weights)))   
                      
  # updated starting values                   
  
  SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
  gamls.upsvR <- gamls.upsv(gamlss1 = NULL, gamlss2, margins, M, l.flist, nstv = NULL, VC, GAM, SP, type = "biv")
  sp <- gamls.upsvR$sp
  start.v <- gamls.upsvR$start.v   
  
  }  
  
  
  ##########################################################################################################################
  ##########################################################################################################################

  if(Model %in% c("BSS", "BPO", "BPO0") || (Model == "B" && margins[2] %in% bl)  )  gamlssfit <- VC$gamlssfit <- TRUE # useful for error massage and joint probs


  func.opt <- func.OPT(margins, M, type = "biv")

  SemiParFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
    
  
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParBIV.fit.post(SemiParFit = SemiParFit, Model = Model, VC = VC, GAM)
                                            
  SemiParFit <- SemiParFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.
 
  y2.m <- y2  
  if(margins[2] == "LN")  y2.m <- exp(y2)

  ##########################################################################################################################


if(gc.l == TRUE) gc()

  ##########################################################################################################################

cov.c(SemiParFit)


  ##########################################################################################################################

gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  # for all.terms
  ##########################################################################################################################


# bit below useful for AT calculations when end is continuous

if( !(Model=="B" && !(margins[2] %in% bl) && end == 2) ) {dataset <- NULL; rm(data) } else { attr(data, "terms") <- NULL; dataset <- data; rm(data) } 


formula.aux <- NULL


if(Model == "BSS"){


# this is for mice # 

formula.aux <- formula

formula.aux[[1]] <- reformulate( attr(terms(formula.aux[[1]]), "term.labels"), response = "ry" )
environment(formula.aux[[1]]) <- environment(formula[[1]])


formula.aux[[2]] <- reformulate( attr(terms(formula.aux[[2]]), "term.labels"), response = "y" )
environment(formula.aux[[2]]) <- environment(formula[[2]])

###



}





L <- list(fit = SemiParFit$fit, dataset = dataset, formula = formula, SemiParFit = SemiParFit, mice.formula = formula.aux,
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, robust = FALSE,
          gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, coef.t = NULL, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,
          OR = SemiParFit.p$OR, 
          GM = SemiParFit.p$GM,    
          n = n, n.sel = n.sel, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X4.d2 = X4.d2, X5.d2 = X5.d2, 
          X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,           
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf11 = SemiParFit.p$edf11,
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6,
          edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8,
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p11 = SemiParFit$fit$p11, p10 = SemiParFit$fit$p10, p01 = SemiParFit$fit$p01, p00 = SemiParFit$fit$p00, 
          p1 = SemiParFit$fit$p1, p2 = SemiParFit$fit$p2,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, etad = SemiParFit$fit$etad,
          etas = SemiParFit$fit$etas, etan = SemiParFit$fit$etan,
          y1 = y1, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = hess, 
          respvec = respvec, inde = inde, 
          qu.mag = qu.mag, 
          sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a,
          sigma  = SemiParFit.p$sigma2, sigma.a  = SemiParFit.p$sigma2.a,
          nu = SemiParFit.p$nu, nu.a = SemiParFit.p$nu.a, Vb.t = SemiParFit.p$Vb.t,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          X2s = X2s, X3s = X3s, p1n=SemiParFit.p$p1n , p2n = SemiParFit.p$p2n, 
          VC = VC, Model = Model, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "NO", tau = SemiParFit.p$tau, 
          tau.a = SemiParFit.p$tau.a, l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE,
          gamlss = gamlss2, BivD2 = BivD2, dof = dof, dof.a = dof, call = cl,
          surv = FALSE, surv.flex = surv.flex)


if(BivD %in% BivD2){       

L$teta1     <- SemiParFit$fit$teta1
L$teta.ind1 <- SemiParFit$fit$teta.ind1   
L$teta2     <- SemiParFit$fit$teta2
L$teta.ind2 <- SemiParFit$fit$teta.ind2   
L$Cop1      <- SemiParFit$fit$Cop1
L$Cop2      <- SemiParFit$fit$Cop2

}


class(L) <- c("SemiParBIV","gjrm")

L

}

