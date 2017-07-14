gamlss <- function(formula, data = list(), weights = NULL, subset = NULL,  
                   margin = "N", surv = FALSE, cens = NULL, 
                   robust = FALSE, rc = 3, lB = NULL, uB = NULL, infl.fac = 1, 
                   rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                   gc.l = FALSE, parscale, extra.regI = "t", gev.par = -0.25,
                   chunk.size = 10000){ 
  
  ##########################################################################################################################
  # preamble 
  ##########################################################################################################################
    
  i.rho <- sp <- qu.mag <- qu.mag1 <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- spgamlss1 <- NULL  
  end <- X2.d2 <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- y1m <- y2m <- NULL
  fp <- FALSE
  gp2 <- gp3 <- 1
  sp1 <- sp2 <- gam2 <- X2 <- sp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL     
  Xd1 <- Xd <- mono.sm.pos <- NULL
  tfc <- NA
  surv.flex <- FALSE
  tempb <- NULL

  m2  <- c("N","N2","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK")
  m3  <- c("DAGUM","SM")
  m1d <- c("PO", "ZTP", "GEVlink")
  m2d <- c("NBI", "NBII","NBIa", "NBIIa","PIG")
  m3d <- c("DEL","SICHEL")
  
  if(margin == "PH" && surv == TRUE) margin <- "cloglog"
  if(margin == "PO" && surv == TRUE) margin <- "logit"  
  
  bl  <- c("probit", "logit", "cloglog")


  ##########################################################################################################################
  if(!is.list(formula)) stop("You must specify a list of one or more equations.")
  l.flist <- length(formula)
  form.check(formula, l.flist, gamlss = TRUE) 
  cl <- match.call()  
  mf <- match.call(expand.dots = FALSE)

  pred.varR <- pred.var(formula, l.flist, gaml = TRUE) 
   
  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  

  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$chunk.size <- mf$gev.par <- mf$surv <- mf$lB <- mf$uB <- mf$robust <- mf$rc <- mf$margin <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- NULL                           
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
 
  n <- dim(data)[1]
          
  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]  
  
  
  if(surv == TRUE && !("(cens)" %in% names(data)) ) stop("You must provide the binary censoring indicator.")
       
       
  if(!("(cens)" %in% names(data))) {cens <- rep(1,dim(data)[1]) 
                        data$cens <- cens
                        names(data)[length(names(data))] <- "(cens)"} else cens <- data[,"(cens)"]                         
  
  M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, robust = robust, extra.regI = extra.regI, margin = margin,
            surv = surv, cens = cens, bl = bl)  
  pream.wm(formula, margins = NULL, M, l.flist, type = "gamls")

  formula.eq1 <- formula[[1]]
  
 ##############################################################  
 ##############################################################  
   
 form.eq12R <- form.eq12(formula.eq1, data, v1, margin, m1d, m2d)   
 
 formula.eq1  <- form.eq12R$formula.eq1
 formula.eq1r <- form.eq12R$formula.eq1r
 y1           <- form.eq12R$y1
 y1.test      <- form.eq12R$y1.test  
 y1m          <- form.eq12R$y1m
   
 if(margin != "GEVlink" && surv == FALSE)                     gam1 <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
 if(margin != "GEVlink" && surv == TRUE && !(margin %in% bl)) gam1 <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights*cens, data=data),list(weights=weights, cens = cens)))
 if(margin == "GEVlink")                                      gam1 <- eval(substitute(gam(formula.eq1, binomial(link = "cloglog"), gamma=infl.fac, weights=weights, data=data),list(weights=weights)))

 #############################################################################################
                                     #############################################################
 if(surv == TRUE && margin %in% bl){ ######### think about this for parametric models ############
  surv.flex <- TRUE                  ########set.seed(1)rm(list=".Random.seed", envir=globalenv())rnorm(dim(data)[1])

  ###################################################################
  ####  cox.ph pre-fit to create response for starting value fit ####
  
  f.eq1 <- form.eq12R$f.eq1
  data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
  tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), data = data, weights = cens),list(cens=cens)))
  data$Sh <- as.vector(mm(predict(tempb, type = "response")))
  
  #data$Sh <- as.vector(mm(rstpm2:::Shat(coxph(Surv(u, delta == 1) ~ 1, data = data, model = TRUE))))
  #does not seem to have an impact
  ###################################################################
  
  cens1 <- ifelse(cens == 0, 1e-07, cens)
  gam1 <- eval(substitute(scam(formula.eq1, gamma=infl.fac, weights=weights*cens1, data=data), list(weights=weights, cens1 = cens1)))
  
  #############################################################
  ### check that we are using the correct smoother of time ####
  
  lsgam1 <- length(gam1$smooth)
  if(lsgam1 == 0) stop("You must use at least a monotonic smooth function of time.")
  
  clsm <- ggr <- NA 
  for(i in 1:lsgam1){ clsm[i] <- class(gam1$smooth[[i]])[1] ### is lsgam1 is and is used for summary when using tensor etc
                      ggr[i]  <- max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn)))
                    }
  
  if( sum(as.numeric(clsm %in% c("mpi.smooth")))==0 ) stop("You must use at least an mpi smooth function of time.")
  if( sum( as.numeric(clsm %in% c("mpi.smooth")) ) != sum( ggr ) ) stop("You must use mpi smooth function(s) of time.")   
  
  ###########################################################

  l.sp1 <- length(gam1$sp)
  if(l.sp1 != 0) sp1 <- gam1$sp
           
  ###########################################################    
  
  if(dim(data)[1] < 2000) sp.c <- 0.2 else sp.c <- 1/sqrt(dim(data)[1])  

  sp1[clsm %in% c("mpi.smooth")] <- sp.c 
  
  gam.call <- gam1$call
  gam.call$sp <- sp1
  gam1 <- eval(gam.call)
  
  ###########################################################

  for(i in 1:lsgam1){ 
  
    if( max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn))) != 0 ) mono.sm.pos <- c(mono.sm.pos, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para) ) 
                                          
  }
     
  
  X1 <- predict(gam1, type = "lpmatrix")
  
  Xd <- Xdpred(gam1, data, v1[1])

  start.v1 <- c( coef(gam1) )

  gam1$y <- data[, v1[1]]
 
 }

 #############################################################################################

 gam1$formula <- formula.eq1r  
 lsgam1 <- length(gam1$smooth)
 
 y1 <- y1.test 
 if( margin %in% c("LN") ) y1 <- log(y1) 
 
 attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
 if( !(surv == TRUE && margin %in% bl) ){
 
     names(gam1$model)[1] <- as.character(formula.eq1r[2])
     X1 <- predict(gam1, type = "lpmatrix")
     l.sp1 <- length(gam1$sp)
     sp1 <- gam1$sp
                                        }
 gp1 <- gam1$nsdf 
 X1.d2 <- dim(X1)[2]
 
 ##############################################################  
 ##############################################################  

log.nu.1 <- log.sig2.1 <- NULL 

if( !(margin %in% c(m1d, bl)) ){

start.snR <- startsn(margin, y1)
    
log.sig2.1 <- start.snR$log.sig2.1; names(log.sig2.1) <- "sigma2.star"
if( margin %in% c(m3) ){ log.nu.1   <- start.snR$log.nu.1;   names(log.nu.1)   <- "nu.star"}     

}

if(margin %in% c(m1d) )    start.v1 <- c( coef(gam1) )
if(margin %in% c(m2,m2d) ) start.v1 <- c( coef(gam1), log.sig2.1           ) 
if(margin %in% c(m3,m3d) ) start.v1 <- c( coef(gam1), log.sig2.1, log.nu.1 ) 


##############################################################  
##############################################################  
  
if(l.flist > 1){ # not used for flexible survival
    
    vo <- list(log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, n = n)
    overall.svGR <- overall.svG(formula, data, ngc = 2, margin, M, vo, gam1, gam2, type = "gaml")
    
    start.v1 <- overall.svGR$start.v 
    X2 <- overall.svGR$X2
    X3 <- overall.svGR$X3  
    X2.d2 <- overall.svGR$X2.d2
    X3.d2 <- overall.svGR$X3.d2
    gp2 <- overall.svGR$gp2
    gp3 <- overall.svGR$gp3
    gam2 <- overall.svGR$gam2
    gam3 <- overall.svGR$gam3
    l.sp2 <- overall.svGR$l.sp2
    l.sp3 <- overall.svGR$l.sp3
    sp2 <- overall.svGR$sp2
    sp3 <- overall.svGR$sp3
    
}
  

##########################################################
# SPs and penalties
##########################################################
##########################################################

spgamlss1 <- c(sp1, sp2, sp3)
GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8) 


if(l.sp1 !=0 || l.sp2 !=0 || l.sp3 !=0) { ##


	L.GAM <- list(l.gam1 = length(coef(gam1)), l.gam2 = length(coef(gam2)), 
	              l.gam3 = length(coef(gam3)), l.gam4 = 0, l.gam5 = 0, 
	              l.gam6 = 0, l.gam7 = 0, l.gam8 = 0)                                 
  
	L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = 0, 
	             l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0)                               
	                      
	qu.mag1 <- S.m(GAM, L.SP, L.GAM)  

                                         } ##

##########################################################
##########################################################

if(missing(parscale)) parscale <- 1   

respvec2 <- list(y1 = y1, univ = 2)               
  
  lsgam2 <- length(gam2$smooth)
  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)


if(robust == TRUE && margin %in% c(m1d, m2d)){

# grid worked out this way seems good enough for the moment
# but would this grid be good for gradient and hessian components?
# maybe we need to look into this again 
# also what about a dinamic grid (that changes with eta and sigma)?

eta.m <- max(predict(gam1, type = "link"))
if( margin %in% c(m2d) ) sigma2.m <- exp(log.sig2.1) else sigma2.m <- 1 # this looks already quite high given that it is from the unconditional fit

if(margin != "ZTP") ygrid <- 0:(max(y1)*100) else ygrid <- 1:(max(y1)*100) 

pdf.test <- distrHsATDiscr(ygrid, eta.m, sigma2.m, 1, margin, y1m, robust = TRUE)$pdf2 > 1e-07

ygrid <- ygrid[1:(table(pdf.test)[2] + 2)]

} else ygrid <- NULL



  VC <- list(lsgam1 = lsgam1, ygrid = ygrid, # why lsgam1? maybe useful outside fitting functions, do not remember, check again
             lsgam2 = lsgam2,
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
             hess = NULL,
             Model = "CC", univ.gamls = TRUE,
             gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = c(margin, margin),
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, m1d = m1d, 
             m2d = m2d, m3d = m3d, bl = bl, triv = FALSE,
             y1m = y1m, y2m = y2m, robust = robust, rc = rc,
             cens = cens, surv = surv,
             lB = lB, uB = uB, gev.par = gev.par,
             chunk.size = chunk.size,
             Xd1 = Xd, 
             mono.sm.pos = mono.sm.pos, 
             surv.flex = surv.flex)

  
  
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################
  if(margin != "GEVlink"){
  
  	if(margin %in% c(m1d, m2d, m2) ) func.opt1 <- bprobgHsContUniv 
  	if(margin %in% c(m3) )           func.opt1 <- bprobgHsContUniv3 
  	if(margin %in% c(bl) )           func.opt1 <- bcontSurvGuniv
  
  }
  
  if(margin == "GEVlink") func.opt1 <- bprobgHsContUnivBIN


    SemiParFit <- SemiParBIV.fit(func.opt = func.opt1, start.v = start.v1, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec2, VC = VC, sp = spgamlss1, qu.mag = qu.mag1) 
                                  
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- gamlss.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)  
                                                          
  y1.m <- y1; if(margin == "LN") y1.m <- exp(y1) 

  SemiParFit <- SemiParFit.p$SemiParFit  

  ##########################################################################################################################

if(gc.l == TRUE) gc()

  ##########################################################################################################################


e.v <- round(min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values), 6)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?gamlss."

if(gradi > 10 && e.v < 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 10 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 10 && e.v < 0)  warning(paste(me2,"\n",me3), call. = FALSE)

  ##########################################################################################################################

gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  # for all.terms when plotting
  ##########################################################################################################################


L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, formula = formula,        
          edf11 = SemiParFit.p$edf11,   ## this is for RE
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, 
          gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, iterlimsp = iterlimsp,
          weights = weights, cens = cens, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
          l.sp7 = l.sp7, l.sp8 = l.sp8,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,  
          sigma2 = SemiParFit.p$sigma2,  
          sigma2.a = SemiParFit.p$sigma2.a, 
          nu = SemiParFit.p$nu, 
          nu.a = SemiParFit.p$nu.a,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
          X6 = X6, X7 = X7, X8 = X8,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
          X7.d2 = X7.d2, X8.d2 = X8.d2,             
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1, Vb.t = SemiParFit.p$Vb.t, coef.t = SemiParFit.p$coef.t,   
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, 
          edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
          edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$etas1, 
          eta3 = SemiParFit$fit$etan1, 
          y1 = y1.m, 
          margins = c(margin, margin),   
          logLik = SemiParFit.p$logLik,
          hess = TRUE,
          qu.mag = qu.mag1, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, 
          gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          VC = VC, magpp = SemiParFit$magpp,
          Cont = "YES",
          l.flist = l.flist, triv = FALSE, univar.gamlss = TRUE, call = cl, gev.par = gev.par,
          ygrid = ygrid,
          r.weights = SemiParFit$fit$d.psi, surv = surv, surv.flex = surv.flex)

class(L) <- c("gamlss","SemiParBIV")


L


}

