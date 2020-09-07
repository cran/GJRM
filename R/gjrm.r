gjrm <- function(formula, data = list(), weights = NULL, subset = NULL,  
                             BivD = "N", margins, Model, dof = 3, ordinal = FALSE,
                             surv = FALSE, cens1 = NULL, cens2 = NULL, cens3 = NULL, dep.cens = FALSE,  
                             gamlssfit = FALSE, fp = FALSE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t", k1.tvc = 0, k2.tvc = 0, 
                             knots = NULL, penCor = "unpen",
                             sp.penCor = 3, Chol = FALSE, gamma = 1, w.alasso = NULL,
                             drop.unused.levels = TRUE, ind.ord = FALSE,
                             min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999){
  
  if(missing(margins)) stop("You must choose the margins' values.")
  if(missing(Model)) stop("You must choose a Model type.")
  
  
  
  if(margins[1] == "PH" && surv == TRUE) margins[1] <- "cloglog"
  if(margins[1] == "PO" && surv == TRUE) margins[1] <- "logit" 
  
  if(margins[2] == "PH" && surv == TRUE) margins[2] <- "cloglog"
  if(margins[2] == "PO" && surv == TRUE) margins[2] <- "logit"   
  
  bl <- c("probit", "logit", "cloglog")  
  
  
  if(surv == FALSE && ordinal == FALSE){
  
  if( (margins[1] %in% bl && margins[2] %in% bl && is.na(margins[3])) || (margins[1] %in% bl && !(margins[2] %in% bl) && Model == "B" && is.na(margins[3]))  ){
            
      L <- eval(substitute(SemiParBIV(formula, data, weights, subset,
                               Model, BivD, margins, dof, gamlssfit, # margins should be probit probit but ok as it should get an error message
                               fp, hess = TRUE, infl.fac, 
                               rinit, rmax, iterlimsp, tolsp,
                               gc.l, parscale, extra.regI, intf = TRUE, 
                               theta.fx = NULL, knots = knots, drop.unused.levels = drop.unused.levels,
                               min.dn = min.dn, min.pr = min.pr, max.pr = max.pr),list(weights=weights)))                               
  
                                                                        }
  
  }
  
  
  
  
  

  
  if(surv == FALSE && ordinal == TRUE){
  
  if( (margins[1] %in% bl && margins[2] %in% bl && is.na(margins[3])) || (margins[1] %in% bl && !(margins[2] %in% bl) && is.na(margins[3]))  ){
  
         
      L <- eval(substitute(CopulaCLM(formula, data, weights, subset,
                               Model, BivD, margins, dof, gamlssfit, 
                               fp, hess = TRUE, infl.fac, 
                               rinit, rmax, iterlimsp, tolsp,
                               gc.l, parscale, extra.regI, intf = TRUE, 
                               theta.fx = NULL, knots = knots, drop.unused.levels = drop.unused.levels, ind.ord = ind.ord,
                               min.dn = min.dn, min.pr = min.pr, max.pr = max.pr),list(weights=weights)))                               
  
                                                                        }
  
  }  
    
  
  
  
  
  
  
  
  
  if( margins[1] %in% bl && !(margins[2] %in% bl) && surv == FALSE && is.na(margins[3]) && Model == "BSS" && ordinal == FALSE){
 

 
      L <- eval(substitute( copulaSampleSel(formula, data, weights, subset,
                              BivD, margins, dof,
                              fp, infl.fac, 
                              rinit, rmax, iterlimsp, tolsp,
                             gc.l, parscale, extra.regI, knots, drop.unused.levels = drop.unused.levels,
                             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr),list(weights=weights)))  
 
                                                                         }





  if(!is.na(margins[3])){
  
  if( margins[1] %in% bl && margins[2] %in% bl && margins[3] %in% bl && surv == FALSE && ordinal == FALSE){

      L <- eval(substitute( SemiParTRIV(formula, data, weights, subset,
                             Model, margins,  
                             penCor, sp.penCor, approx = FALSE, Chol, 
                             infl.fac, gamma, w.alasso, 
                             rinit, rmax, 
                             iterlimsp, tolsp,
                             gc.l, parscale, extra.regI, knots, drop.unused.levels = drop.unused.levels,
                             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr),list(weights=weights))) 

                                                                      }

  }
    

 
  if(  ( !(margins[1] %in% bl) || surv == TRUE ) && ordinal == FALSE ){
    
  ##########################################################################################################################
  # preamble
  ##########################################################################################################################  
  
  robust <- FALSE; t.c = 3
  sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- gam1 <- gam2 <- y1m <- y2m <- indexTeq1 <- indexTeq2 <- NULL  
  i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- NULL    
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL   
  c11 <- c10 <- c01 <- c00 <- NA
  
  Sl.sf <- NULL
  sp.method <- "perf"
  
  Xd1 <- Xd2 <- mono.sm.pos1 <- mono.sm.pos2 <- mono.sm.pos <- NULL
  surv.flex <- FALSE
  
  Deq1 <- pos.pbeq1 <- Deq2 <- pos.pbeq2 <- list()
  
  ###################################
  
  BivD2 <- c("C0C90","C0C270","C180C90","C180C270",
             "J0J90","J0J270","J180J90","J180J270",
             "G0G90","G0G270","G180G90","G180G270")
             
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM","T","PL","HO")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180", BivD2)
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK","GP","GPII","GPo")
  m3   <- c("DAGUM","SM","TW")
  m1d  <- c("PO", "ZTP")
  m2d  <- c("NBI", "NBII","PIG","DGP","DGPII")
  m3d  <- c("DEL","SICHEL")
  
  ct  <- data.frame( c(opc), c(1:14,55,56,57,60,61) )
  cta <- data.frame( c(opc), c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2,60,61) )     
  
  
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
  l.flist <- length(formula)

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
  
  mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$dep.cens <- mf$ordinal <- mf$Model <- mf$knots <- mf$k1.tvc <- mf$k2.tvc <- mf$surv <- mf$BivD <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  mf$drop.unused.levels <- drop.unused.levels 
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
 
  n <- dim(data)[1]
        
  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"] 
 
  if(surv == TRUE && !("(cens1)" %in% names(data)) && margins[1] %in% bl ) stop("You must provide the first binary censoring indicator.")
  if(surv == TRUE && !("(cens2)" %in% names(data)) && margins[2] %in% bl ) stop("You must provide the second binary censoring indicator.")
 
 
  if(!("(cens1)" %in% names(data))) {cens1 <- rep(0,dim(data)[1]) 
                        data$cens1 <- cens1
                        names(data)[length(names(data))] <- "(cens1)"} else cens1 <- data[,"(cens1)"]                         

  if(!("(cens2)" %in% names(data))) {cens2 <- rep(0,dim(data)[1]) 
                        data$cens2 <- cens2
                        names(data)[length(names(data))] <- "(cens2)"} else cens2 <- data[,"(cens2)"]  
                        
  if(!("(cens3)" %in% names(data))) {cens3 <- rep(0,dim(data)[1]) 
                        data$cens3 <- cens3
                        names(data)[length(names(data))] <- "(cens3)"} else cens3 <- data[,"(cens3)"]                          
                        
  M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, BivD = BivD, bl = bl, 
            robust = robust, opc = opc, extra.regI = extra.regI, margins = margins, BivD2 = BivD2, dof = dof,
            surv = surv, c1 = cens1, c2 = cens2, c3 = cens3, dep.cens = dep.cens) 
 
  M$K1 <- NULL
 
  pream.wm(formula, margins, M, l.flist)

  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
    
 ##############################################################  
 # Equation 1
 ##############################################################  
   
 form.eq12R <- form.eq12(formula.eq1, data, v1, margins[1], m1d, m2d)   
 
 formula.eq1  <- form.eq12R$formula.eq1
 formula.eq1r <- form.eq12R$formula.eq1r
 y1           <- form.eq12R$y1
 y1.test      <- form.eq12R$y1.test 
 y1m          <- form.eq12R$y1m

 if(surv == FALSE)                                                   gam1 <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))
 if(surv == TRUE && margins[1] %in% c(m2,m3) && margins[2] %in% bl ) gam1 <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights))) else{
  
         if(surv == TRUE && !(margins[1] %in% bl)) gam1 <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights*cens1, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights, cens1 = cens1)))
         }


if(surv == TRUE && margins[1] %in% bl){ 

  surv.flex <- TRUE                  

  f.eq1 <- form.eq12R$f.eq1
  data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
  tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), data = data, weights = cens1, drop.unused.levels = drop.unused.levels),list(cens1=cens1)))
  data$Sh <- as.vector(mm(predict(tempb, type = "response"), min.pr = min.pr, max.pr = max.pr))
  
  cens11 <- ifelse(cens1 == 0, 1e-07, cens1)
  gam1 <- eval(substitute(scam(formula.eq1, gamma=infl.fac, weights=weights*cens11, data=data), list(weights=weights, cens11 = cens11)))
  
  lsgam1 <- length(gam1$smooth)
  if(lsgam1 == 0) stop("You must use at least a monotonic smooth function of time in the first equation.")
  
  clsm <- ggr <- NA 
  for(i in 1:lsgam1){ clsm[i] <- class(gam1$smooth[[i]])[1] 
                      #ggr[i]  <- max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn)))
                    }
  
  if( sum(as.numeric(clsm[1] %in% c("mpi.smooth")))==0 ) stop("You must have a monotonic smooth of time and it has to be the first to be included.")

  
  #if( sum(as.numeric(clsm %in% c("mpi.smooth")))==0 ) stop("You must use at least an mpi smooth function of time in the first equation.")
  #if( sum( as.numeric(clsm %in% c("mpi.smooth")) ) != sum( ggr ) ) stop("You must use mpi smooth function(s) of time in the first equation.")   
  
  l.sp1 <- length(gam1$sp)
  if(l.sp1 != 0) sp1 <- gam1$sp
           
  ###########################################################    
  
  #if(dim(data)[1] < 2000) sp.c <- 0.2 else sp.c <- 1/sqrt(dim(data)[1])  
  #sp1[clsm %in% c("mpi.smooth")] <- sp.c 

  
  #sp.c <- 1 # there is really no way for making this step better at the moment
  sp1[1] <- 1 #sp.c 

  gam.call <- gam1$call
  gam.call$sp <- sp1
  gam1 <- eval(gam.call)
  
  ###########################################################
  j <- 1
  
  for(i in 1:lsgam1){ 
  
    if( max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$term))) != 0 && clsm[i] == "mpi.smooth" ) mono.sm.pos1 <- c(mono.sm.pos1, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para) ) 
    
    #if( max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn))) != 0 && clsm[i] != "mpi.smooth" ){ 
    #
    #
    #                                                            if( clsm[i] != "pspline.smooth" && k1.tvc !=0) stop("You have to use a ps smooth to allow for doubly penalised tvc terms in the first eq.")
    #                                                                                                                            
    #                                                            if( clsm[i] == "pspline.smooth"){
    #                                                            
    #                                                            pos.pbeq1[[j]] <- c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para)
    #                                                            indexTeq1      <- c(indexTeq1, pos.pbeq1[[j]] ) 
    #                                                            Deq1[[j]]      <- diff(diag(length(pos.pbeq1[[j]])), differences = 1)
    #                                                            j <- j + 1
    #                                                            
    #                                                            }
    #                                                            
    #                                                            
    #                                                                                        }    
    

                    }
       
  
  
  
  X1  <- predict(gam1, type = "lpmatrix")
  if( !is.null(indexTeq1) && k1.tvc !=0){ if(range(X1[, indexTeq1])[1] < 0) stop("Check design matrix for smooth(s) of tvc term(s) in eq. 1.")}

  Xd1 <- Xdpred(gam1, data, v1[1])

  gam1$y <- data[, v1[1]]

  st.v1 <- c( gam1$coefficients )


  if(!is.null(indexTeq1)){
  
     st.v1[mono.sm.pos1] <- exp(st.v1[mono.sm.pos1])
     while( range(Xd1%*%st.v1)[1] < 0 ) st.v1[indexTeq1] <- 0.999*st.v1[indexTeq1]   
           
     gam1$coefficients <- gam1$coefficients.t <- st.v1
     gam1$coefficients.t[mono.sm.pos1] <- exp(gam1$coefficients.t[mono.sm.pos1]) 
  }



 
 }



 gam1$formula <- formula.eq1r  
 lsgam1 <- length(gam1$smooth)
 
 y1 <- y1.test 
 if( margins[1] %in% c("LN") ) y1 <- log(y1) 
 
 attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
 if( !(surv == TRUE && margins[1] %in% bl) ){
 
     names(gam1$model)[1] <- as.character(formula.eq1r[2])
     X1 <- predict(gam1, type = "lpmatrix")
     l.sp1 <- length(gam1$sp)
     sp1 <- gam1$sp
                                        }
 gp1 <- gam1$nsdf 
 X1.d2 <- dim(X1)[2]


 
 ##############################################################
 # Equation 2 
 ##############################################################  

 form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], m1d, m2d)   
 
 formula.eq2  <- form.eq12R$formula.eq1
 formula.eq2r <- form.eq12R$formula.eq1r
 y2           <- form.eq12R$y1
 y2.test      <- form.eq12R$y1.test 
 y2m          <- form.eq12R$y1m

 if(surv == FALSE)                         gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))
 if(surv == TRUE && !(margins[2] %in% bl)) gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights*cens2, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights, cens2 = cens2)))

if(surv == TRUE && margins[2] %in% bl){ 
  surv.flex <- TRUE                  

  f.eq2 <- form.eq12R$f.eq1
  data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
  tempb <- eval(substitute(gam(f.eq2, family = cox.ph(), data = data, weights = cens2, drop.unused.levels = drop.unused.levels),list(cens2=cens2)))
  data$Sh <- as.vector(mm(predict(tempb, type = "response"), min.pr = min.pr, max.pr = max.pr))
  
  cens22 <- ifelse(cens2 == 0, 1e-07, cens2)
  gam2 <- eval(substitute(scam(formula.eq2, gamma=infl.fac, weights=weights*cens22, data=data), list(weights=weights, cens22 = cens22)))
  
  lsgam2 <- length(gam2$smooth)
  if(lsgam2 == 0) stop("You must use at least a monotonic smooth function of time in the second equation.")
  
  clsm <- ggr <- NA 
  for(i in 1:lsgam2){ clsm[i] <- class(gam2$smooth[[i]])[1] 
                      #ggr[i]  <- max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$vn)))
                    }
  
  
  if( sum(as.numeric(clsm[1] %in% c("mpi.smooth")))==0 ) stop("You must have a monotonic smooth of time and it has to be the first to be included.")

  
  #if( sum(as.numeric(clsm %in% c("mpi.smooth")))==0 ) stop("You must use at least an mpi smooth function of time in the second equation.")
  #if( sum( as.numeric(clsm %in% c("mpi.smooth")) ) != sum( ggr ) ) stop("You must use mpi smooth function(s) of time in the second equation.")   
  
  l.sp2 <- length(gam2$sp)
  if(l.sp2 != 0) sp2 <- gam2$sp
           
  ###########################################################    
  
  #if(dim(data)[1] < 2000) sp.c <- 0.2 else sp.c <- 1/sqrt(dim(data)[1])  
  #sp2[clsm %in% c("mpi.smooth")] <- sp.c 

  #sp.c <- 1
  sp2[1] <- 1 #sp.c 
  
  
  gam.call <- gam2$call
  gam.call$sp <- sp2
  gam2 <- eval(gam.call)
  
  ###########################################################

  j <- 1
  for(i in 1:lsgam2){ 
  
  
    if( max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) != 0 && clsm[i] == "mpi.smooth" ) mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para) )   
    
    #if( max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$vn))) != 0 && clsm[i] != "mpi.smooth" ){ 
    #
    #
    #                                                            if( clsm[i] != "pspline.smooth" && k2.tvc !=0) stop("You have to use a ps smooth to allow for doubly penalised tvc terms in the second eq.")
    #                                                                                                                            
    #                                                            if( clsm[i] == "pspline.smooth"){
    #                                                            
    #                                                            pos.pbeq2[[j]] <- c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para)
    #                                                            indexTeq2      <- c(indexTeq2, pos.pbeq2[[j]] ) 
    #                                                            Deq2[[j]]      <- diff(diag(length(pos.pbeq2[[j]])), differences = 1)
    #                                                            j <- j + 1
    #                                                            
    #                                                            }
    #
    #                                                                                            }     

                    }
    
  X2  <- predict(gam2, type = "lpmatrix")
  
  if( !is.null(indexTeq2) && k2.tvc !=0){ if(range(X2[, indexTeq2])[1] < 0) stop("Check design matrix for smooth(s) of tvc term(s) in eq. 2.")}
  
  Xd2 <- Xdpred(gam2, data, v2[1])

  gam2$y <- data[, v2[1]]
 
  st.v2 <- c( gam2$coefficients )


  if(!is.null(indexTeq2)){
  
     st.v2[mono.sm.pos2] <- exp(st.v2[mono.sm.pos2])
     while( range(Xd2%*%st.v2)[1] < 0 ) st.v2[indexTeq2] <- 0.999*st.v2[indexTeq2]   
           
     gam2$coefficients <- gam2$coefficients.t <- st.v2
     gam2$coefficients.t[mono.sm.pos2] <- exp(gam2$coefficients.t[mono.sm.pos2]) 
  }

 
 
 
 
 
 }



 gam2$formula <- formula.eq2r  
 lsgam2 <- length(gam2$smooth)
 
 y2 <- y2.test 
 if( margins[2] %in% c("LN") ) y2 <- log(y2) 
 
 attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
 if( !(surv == TRUE && margins[2] %in% bl) ){
 
     names(gam2$model)[1] <- as.character(formula.eq2r[2])
     X2 <- predict(gam2, type = "lpmatrix")
     l.sp2 <- length(gam2$sp)
     sp2 <- gam2$sp
                                        }
 gp2 <- gam2$nsdf 
 X2.d2 <- dim(X2)[2]

  
#################################################################
# Starting value for dependence parameter (and dof for T if used)
#################################################################

res1 <- residuals(gam1)
res2 <- residuals(gam2)
 
ass.s <- cor(res1, res2, method = "kendall")
ass.s <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))

i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)

dof.st <- log(dof - 2) 
names(dof.st) <- "dof.star"   
                           
##############################################################
# Other starting values + overall
##############################################################
           
if( !(margins[1] %in% c(m1d,bl)) ){

start.snR <- startsn(margins[1], y1)
    
log.sig2.1 <- start.snR$log.sig2.1; names(log.sig2.1) <- "sigma1.star"
if( margins[1] %in% c(m3) ){ log.nu.1   <- start.snR$log.nu.1;   names(log.nu.1)   <- "nu.1.star"}     

}

if( !(margins[2] %in% c(m1d,bl)) ){

start.snR <- startsn(margins[2], y2)
    
log.sig2.2 <- start.snR$log.sig2.1; names(log.sig2.2) <- "sigma2.star"
if( margins[2] %in% c(m3) ){ log.nu.2   <- start.snR$log.nu.1;   names(log.nu.2)   <- "nu.2.star"}     

}

vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2.2 = log.sig2.2, log.nu.2 = log.nu.2, log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, dof.st = dof.st, n = n, drop.unused.levels = drop.unused.levels )

start.v <- overall.sv(margins, M, vo)
   			
##############################################################  
# starting values for case of predictors on all parameters
##############################################################  
  
    if(l.flist > 2){
    
    overall.svGR <- overall.svG(formula, data, ngc = 2, margins, M, vo, gam1, gam2, knots = knots)
                                
    
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
    
    }
    

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
  

##########################################################
# general lists
##########################################################

if(missing(parscale)) parscale <- 1   

  respvec <- respvec2 <- respvec3 <- list(y1 = y1, y2 = y2,
                                          y1.y2 = NULL, y1.cy2 = NULL, 
                                          cy1.y2 = NULL, cy1.cy2 = NULL, 
                                          cy1 = NULL, cy = NULL, univ = 0)
 
  my.env <- new.env()
  my.env$signind <- 1 # this is for mixed copulae

  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)




if(surv == TRUE && dep.cens == FALSE){

if((surv == TRUE && margins[1] %in% bl && margins[2] %in% bl) || (surv == TRUE && margins[1] %in% m2 && margins[2] %in% m2) ){

c11 <- cens1*cens2
c10 <- cens1*(1-cens2)
c01 <- (1-cens1)*cens2
c00 <- (1-cens1)*(1-cens2)

}

if(surv == TRUE && margins[1] %in% c(m2,m3) && margins[2] %in% bl){

c11 <- cens2
c10 <- 1 - cens2
c01 <- NULL
c00 <- NULL

}


}




if(surv == TRUE && dep.cens == TRUE){


c11 <- NULL
c10 <- cens1
c01 <- cens2 # (1-cens1)
c00 <- cens3 # NULL, in case of A cens then 1 - cens1 - cens2


}






#my.env      <- new.env()
my.env$k1   <- k1.tvc
my.env$k2   <- k2.tvc


  VC <- list(lsgam1 = lsgam1, indexTeq1 = indexTeq1, indexTeq2 = indexTeq2, 
             lsgam2 = lsgam2, Deq1 = Deq1, pos.pbeq1 = pos.pbeq1, Deq2 = Deq2, pos.pbeq2 = pos.pbeq2,
             lsgam3 = lsgam3, robust = FALSE, sp.fixed = NULL,
             lsgam4 = lsgam4, Sl.sf = Sl.sf, sp.method = sp.method,
             lsgam5 = lsgam5, K1 = NULL,
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
             l.sp8 = l.sp8, my.env = my.env,
             infl.fac = infl.fac,
             weights = weights,
             fp = fp, 
             gamlssfit = gamlssfit,
             hess = NULL,
             Model = "CC", univ.gamls = FALSE,
             end = end,
             BivD = BivD, nCa = nCa,
             nC = nC, gc.l = gc.l, 
             n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, 
             m1d = m1d, m2d = m2d, m3d = m3d, 
             bl = bl, triv = FALSE,
             y1m = y1m, y2m = y2m, 
             tc = t.c,
             i.rho = i.rho, dof = dof,
             dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct,
             zerov = -10,
             c11 = c11,
             c10 = c10,
             c01 = c01,
             c00 = c00, surv = surv,
             Xd1 = Xd1, Xd2 = Xd2,
             mono.sm.pos1 = mono.sm.pos1, mono.sm.pos2 = mono.sm.pos2, 
             surv.flex = surv.flex,
             mono.sm.pos = mono.sm.pos, gp2.inf = NULL,
             informative = "no",
             zero.tol = 1e-02,
             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr) # original n only needed in SemiParBIV.fit
  
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  ##########################################################################################################################
  # GAMLSS fit
  ##########################################################################################################################
  ##########################################################################################################################

if(gamlssfit == TRUE){ 

  form.gamlR <- form.gaml(formula, l.flist, M)

  surv1 <- surv2 <- surv
  
  if(surv == TRUE && margins[1] %in% c(m2,m3) && margins[2] %in% bl ) surv1 <- FALSE 

  gamlss1 <- eval(substitute(gamlss(form.gamlR$formula.gamlss1, data = data, weights = weights, subset = subset,  
                   margin = margins[1], surv = surv1, cens = cens1, infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, k.tvc = k1.tvc, drop.unused.levels = drop.unused.levels), list(weights=weights,cens1=cens1)))

  gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = subset,  
                   margin = margins[2], surv = surv2, cens = cens2, infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, k.tvc = k2.tvc, drop.unused.levels = drop.unused.levels), list(weights=weights,cens2=cens2)))   
                      
  # updated starting values   

  SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
  gamls.upsvR <- gamls.upsv(gamlss1, gamlss2, margins, M, l.flist, nstv = names(start.v), VC, GAM, SP)
  sp <- gamls.upsvR$sp
  start.v <- gamls.upsvR$start.v 

}

  ##########################################################################################################################
  ##########################################################################################################################
  
  func.opt <- func.OPT(margins, M)  
  
  SemiParFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
    
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)                                     
 
  y1.m <- y1; if(margins[1] == "LN") y1.m <- exp(y1) 
  y2.m <- y2; if(margins[2] == "LN") y2.m <- exp(y2)

  SemiParFit <- SemiParFit.p$SemiParFit  

  if(gc.l == TRUE) gc()

  ##########################################################################################################################


cov.c(SemiParFit)


  ##########################################################################################################################
gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  
  # for all.terms
  ##########################################################################################################################


L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, gamlss1 = gamlss1, gamlss2 = gamlss2, formula = formula, robust = FALSE,     
          edf11 = SemiParFit.p$edf11, surv = surv, 
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, coef.t = SemiParFit.p$coef.t, 
          iterlimsp = iterlimsp,
          weights = weights, cens1 = cens1, cens2 = cens2, cens3 = cens3,
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,  
          sigma21 = SemiParFit.p$sigma21, sigma22 = SemiParFit.p$sigma22, 
          sigma21.a = SemiParFit.p$sigma21.a, sigma22.a = SemiParFit.p$sigma22.a,
          sigma1 = SemiParFit.p$sigma21, sigma2 = SemiParFit.p$sigma22, 
          sigma1.a = SemiParFit.p$sigma21.a, sigma2.a = SemiParFit.p$sigma22.a,                    
          nu1 = SemiParFit.p$nu1, nu2 = SemiParFit.p$nu2, 
          nu1.a = SemiParFit.p$nu1.a, nu2.a = SemiParFit.p$nu2.a,
          dof.a = SemiParFit.p$dof.a, dof = SemiParFit.p$dof,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,            
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, edf7 = SemiParFit.p$edf7,
          edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, edf1.7 = SemiParFit.p$edf1.7, 
          edf1.8 = SemiParFit.p$edf1.8, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
          etad=SemiParFit$fit$etad, etas1 = SemiParFit$fit$etas1, etas2 = SemiParFit$fit$etas2,
          y1 = y1.m, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, 
          respvec = respvec, hess = TRUE,
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          VC = VC, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "YES",
          tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE,
          BivD2 = BivD2, call = cl, surv = surv, surv.flex = surv.flex,
          Vb.t = SemiParFit.p$Vb.t, coef.t = SemiParFit.p$coef.t)
  
if(BivD %in% BivD2){       

L$teta1     <- SemiParFit$fit$teta1
L$teta.ind1 <- SemiParFit$fit$teta.ind1   
L$teta2     <- SemiParFit$fit$teta2
L$teta.ind2 <- SemiParFit$fit$teta.ind2   
L$Cop1      <- SemiParFit$fit$Cop1
L$Cop2      <- SemiParFit$fit$Cop2

}          

class(L) <- c("gjrm","SemiParBIV")


}

L


}

