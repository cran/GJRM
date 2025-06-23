gamlss <- function (formula, data = list(), weights = NULL, subset = NULL, 
          offset =  NULL, family = "N", cens = NULL, type.cens = "R", 
          ub.t = NULL, left.trunc = 0, robust = FALSE, rc = 3, lB = NULL, uB = NULL, 
          infl.fac = 1, rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07, 
          gc.l = FALSE, parscale, gev.par = -0.25, 
          chunk.size = 10000, knots = NULL, informative = "no", 
          inform.cov = NULL, family2 = "-cloglog", fp = FALSE, sp = NULL, 
          drop.unused.levels = TRUE, siginit = NULL, shinit = NULL, 
          sp.method = "perf", hrate = NULL, d.lchrate = NULL, d.rchrate = NULL, 
          d.lchrate.td = NULL, d.rchrate.td = NULL, truncation.time = NULL, 
          min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.9999999, ygrid.tol = 1e-08){
          
          
  if(length(data) == 0) stop("A data frame must be provided.")

  surv <- FALSE        
  extra.regI <- "t" # default
  k.tvc <- 0 # default
  
  mar1surv <- family # for informative cross validation
  mar2surv <- family2 # same here 
  
  upperB <- ub.t
  
        
  if (!is.null(sp)) sp.fixed <- sp
  else sp.fixed <- NULL
  ygrid.tol <- ifelse(ygrid.tol < 1e-160, 1e-160, ygrid.tol)
  v.rB <- upperB
  v.td <- truncation.time
  r.type <- "a"
  i.rho <- sp <- qu.mag <- qu.mag1 <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- spgamlss1 <- indexT <- test.sv.inf <- fgam <- NULL
  end <- X2.d2 <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- gam9 <- y1m <- y2m <- Xi <- X1ni <- X2ni <- gam1TW <- NULL
  gp2 <- gp3 <- 1
  sp1 <- sp2 <- gam2 <- X2 <- sp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- inde.inf2 <- inde.inf1 <- infsetupR <- NULL
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL
  Xd2 <- mono.sm.pos2 <- Sl.sf <- rangeSurv <- NULL
  Xd1 <- Xd <- mono.sm.pos <- gam1TW <- NULL
  indvU <- indvR <- indvL <- indvI <- NULL
  indvUT <- indvRT <- indvLT <- indvIT <- NULL
  Sl.sf1 <- Sl.sf2 <- Sl.sf3 <- list()
  tfc <- no.pb <- NA
  surv.flex <- FALSE
  tempb <- NULL
  D <- pos.pb <- list()
  Gmat12 <- Hmat12 <- NULL
  m2 <- c("tN", "N", "GU", "rGU", "LO", "LN", "WEI", "IG", "GA", "BE", "FISK", "GP", "GPII", "GPo")
  m3 <- c("DAGUM", "SM", "TW")
  m1d <- c("P", "tP", "GEVlink", "DGP0")
  m2d <- c("tNBI", "tNBII", "tPIG", "NBI", "NBII", "PIG", "DGP", "DGPII")
  m3d <- c("DEL", "SICHEL")
  
  
  if(left.trunc > 0 && family %in% c("DGP0", "DGP", "DGPII")) stop("left.trunc option not implemented yet for the chosen distribution. \nGet in touch for more info.")
  
  if( family %in% c("-cloglog", "-logit", "-probit") ) surv <- TRUE
  
  if(family == "-cloglog"     ) family <- "cloglog"
  if(family == "-logit"     ) family <- "logit"
  if(family == "-probit") family <- "probit"

  bl <- c("probit", "logit", "cloglog")
  

  if (surv == TRUE && informative == "yes") {
    if(family2 ==      "-cloglog") family2 <- "cloglog"
    if(family2 ==      "-logit") family2 <- "logit"
    if(family2 == "-probit") family2 <- "probit"
  }
  
  
  if (!is.list(formula)) stop("You must specify a list of one or more equations.")
  l.flist <- length(formula)
  if (l.flist != 3 && family == "TW") stop("You must specify three equations.")
    
  if (family == "TW") {
    formulamgcv <- formula
    formulamgcv[[2]] <- formula[[3]]
    formulamgcv[[3]] <- formula[[2]]
  }
  if (surv == TRUE && family %in% bl && informative == "yes") {
    form.check(formula, l.flist, gamlss = FALSE)
    if (all.equal(formula[[1]], formula[[2]]) != TRUE) stop("The two formulae have to be the same. Get in touch for more info.")
  }
  else form.check(formula, l.flist, gamlss = TRUE)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  if (surv == TRUE && family %in% bl && informative == "yes") pred.varR <- pred.var(formula, l.flist, gaml = TRUE, informative = "yes")
  else pred.varR <- pred.var(formula, l.flist, gaml = TRUE)
  v1 <- pred.varR$v1
  v2 <- pred.varR$v2
  pred.n <- pred.varR$pred.n
  if (!is.null(v.rB)) pred.n <- c(pred.n, v.rB)
  if (!is.null(v.td)) pred.n <- c(pred.n, v.td)
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + "))
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula
  mf$left.trunc <- mf$ub.t <- mf$ygrid.tol <- mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$hrate <- mf$d.lchrate <- mf$d.rchrate <- mf$upperB <- mf$type.cens <- mf$sp.method <- mf$siginit <- mf$shinit <- mf$ordinal <- mf$sp <- mf$fp <- mf$lB <- mf$uB <- mf$family2 <- mf$informative <- mf$inform.cov <- mf$knots <- mf$k.tvc <- mf$chunk.size <- mf$gev.par <- mf$surv <- mf$robust <- mf$rc <- mf$family <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- NULL
  mf$d.lchrate.td <- mf$d.rchrate.td <- mf$truncation.time <- NULL
  mf$drop.unused.levels <- drop.unused.levels
  if (surv == TRUE && type.cens %in% c("I", "mixed")) mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  if (gc.l == TRUE) gc()
  if (surv == TRUE && !("(cens)" %in% names(data))) stop("You must provide the censoring indicator.")
  if (surv == TRUE && !any(data[, "(cens)"] %in% c("I", "IT")) && !is.null(upperB)) stop("Argument upperB is not needed when there are no interval censored \n and/or interval censored left-truncated observations.")
  if (surv == TRUE && !any(data[, "(cens)"] %in% c("UT", "RT", "LT", "IT")) && !is.null(truncation.time)) stop("Argument truncation.time is not needed when there are no left-truncated observations.")
  if (!("(cens)" %in% names(data))) {
    cens <- rep(1, dim(data)[1])
    data$cens <- cens
    names(data)[length(names(data))] <- "(cens)"
  }
  else cens <- data[, "(cens)"]
  if (surv == TRUE && type.cens %in% c("I")) {
    data[cens == 1, v.rB] <- data[cens == 1, v1[1]]
    actual.NAs = as.numeric(which(apply(apply(data, 1, is.na), 2, any)))
    data <- na.omit(data)
    if (length(actual.NAs) > 0) cens = cens[-actual.NAs]
  }
  if (surv == TRUE && type.cens %in% c("mixed")) {
    if (!is.null(truncation.time)) data[!(cens %in% c("UT", "LT", "RT", "IT")), v.td] <- data[!(cens %in% c("UT", "LT", "RT", "IT")), v1[1]]
    if (any(unique(cens) %in% c("I", "IT"))) data[!(cens %in% c("I", "IT")), v.rB] <- data[!(cens %in% c("I", "IT")), v1[1]]
    actual.NAs = as.numeric(which(apply(apply(data, 1, is.na), 2, any)))
    data <- na.omit(data)
    if (length(actual.NAs) > 0) cens = cens[-actual.NAs]
  }
  if (!("(weights)" %in% names(data))) {
    weights <- rep(1, dim(data)[1])
    data$weights <- weights
    names(data)[length(names(data))] <- "(weights)"
  }
  else weights <- data[, "(weights)"]
  
  if (!("(offset)" %in% names(data))) {
    offset <- rep(0, dim(data)[1])
    data$offset <- offset
    names(data)[length(names(data))] <- "(offset)"
  }
  else offset <- data[, "(offset)"]  
  
  
  M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, robust = robust, extra.regI = extra.regI, margin = family, 
            surv = surv, cens = cens, bl = bl, informative = informative, 
            list.inf.cov = inform.cov, sp.method = sp.method, type.cens = type.cens, v.rB = v.rB, v.td = v.td)
  pream.wm(formula, margins = NULL, M, l.flist, type = "gamls")
  
  if(surv == TRUE)                   data[, v1[1]] <- ifelse(data[, v1[1]] < 0.0001, 0.0001, data[, v1[1]])
  if(surv == TRUE && !is.null(v.td)) data[, v.td]  <- ifelse(data[, v.td]  < 0.0001, 0.0001, data[, v.td] )
  if(surv == TRUE && !is.null(v.rB)) data[, v.rB]  <- ifelse(data[, v.rB]  < 0.0001, 0.0001, data[, v.rB] )
  
  
  if (surv == TRUE && !is.null(v.rB) && !is.null(v.td)) rangeSurv <- range(c(data[, v1[1]], data[, v.rB], data[, v.td]))
  if (surv == TRUE && is.null(v.rB) && !is.null(v.td)) rangeSurv <- range(c(data[, v1[1]], data[, v.td]))
  if (surv == TRUE && !is.null(v.rB) && is.null(v.td)) rangeSurv <- range(c(data[, v1[1]], data[, v.rB]))
  if (surv == TRUE && is.null(v.rB) && is.null(v.td)) rangeSurv <- range(data[, v1[1]])
  u1resp <- data[, v1[1]]
  n <- dim(data)[1]
  formula.eq1 <- formula[[1]]
  if (surv == TRUE && informative == "yes") formula.eq2 <- formula[[2]]
  form.eq12R <- form.eq12(formula.eq1, data, v1, family, m1d, m2d)
  formula.eq1 <- form.eq12R$formula.eq1
  formula.eq1r <- form.eq12R$formula.eq1r
  y1 <- form.eq12R$y1
  y1.test <- form.eq12R$y1.test
  y1m <- form.eq12R$y1m
  if (surv == TRUE && family2 %in% bl && informative == "yes") {
    form.eq12R2 <- form.eq12(formula.eq2, data, v2, family2, m1d, m2d)
    formula.eq2 <- form.eq12R2$formula.eq1
    formula.eq2r <- form.eq12R2$formula.eq1r
    y2 <- form.eq12R2$y1
    y2.test <- form.eq12R2$y1.test
    y2m <- form.eq12R2$y1m
  }
  if (family %in% c("GP", "GPII", "GPo", "DGP", "DGPII")) {
    est.ob <- try(gpd.fit(y1, threshold = 0, siglink = exp, show = FALSE), silent = TRUE)
    if (inherits(est.ob, "try-error")) est.ob <- try(gpd.fit(y1, threshold = 0, siglink = exp, show = FALSE, siginit = mean(y1 + mean(y1))/2, shinit = 0.0025), silent = TRUE)
    if (inherits(est.ob, "try-error")) st.v <- c(0.0025, log(mean((y1 + mean(y1))/2)))
    if (!inherits(est.ob, "try-error")) st.v <- c(est.ob$mle[2], est.ob$mle[1])
    esres <- st.v
    if (esres[1] < 0) {
      esresT <- try(resp.check(y1, family, plots = FALSE, print.par = TRUE, loglik = FALSE, os = FALSE, i.f = TRUE, left.trunc = left.trunc), silent = TRUE)
      if (!inherits(esresT, "try-error")) esres <- esresT
    }
    if (esres[1] < 0) esres[1] <- 0.001
    if (family %in% c("GPII", "GPo")) esres[1] <- log(esres[1] + 0.5)
    if (family %in% c("DGPII")) esres[1] <- log(esres[1])
    if (!is.null(siginit)) esres[2] <- siginit
    if (!is.null(shinit)) esres[1] <- shinit
    data$estobXiGP <- esres[1] + y1/(mean(y1) * 100)
  }
  
  
  if (family != "GEVlink" && surv == FALSE) {
    gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, weights = weights, data = data, offset = offset, knots = knots, drop.unused.levels = drop.unused.levels), list(weights = weights, offset = offset)))
    offset <- gam1$offset
    if (sp.method != "perf") {
      gam1ff <- eval(substitute(gam(formula.eq1, gamma = infl.fac, weights = weights, data = data, offset = offset, knots = knots, drop.unused.levels = drop.unused.levels, fit = FALSE), list(weights = weights, offset = offset)))
      Sl.sf1 <- Sl.setup(gam1ff); offset <- gam1ff$offset
      rm(gam1ff)
    }
    
     
    
  }
  if (family != "GEVlink" && surv == TRUE && !(family %in% bl)) {
    gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, weights = weights * cens, data = data, knots = knots, drop.unused.levels = drop.unused.levels), list(weights = weights, cens = cens)))
    if (sp.method != "perf") {
      gam1ff <- eval(substitute(gam(formula.eq1, gamma = infl.fac, weights = weights * cens, data = data, knots = knots, drop.unused.levels = drop.unused.levels, fit = FALSE), list(weights = weights, cens = cens)))
      Sl.sf1 <- Sl.setup(gam1ff)
      rm(gam1ff)
    }
  }
  if (family == "GEVlink") {
    gam1 <- eval(substitute(gam(formula.eq1, binomial(link = "cloglog"), gamma = infl.fac, weights = weights, data = data, knots = knots, drop.unused.levels = drop.unused.levels), list(weights = weights)))
    if (sp.method != "perf") {
      gam1ff <- eval(substitute(gam(formula.eq1, binomial(link = "cloglog"), gamma = infl.fac, weights = weights, data = data, knots = knots, drop.unused.levels = drop.unused.levels, fit = FALSE), list(weights = weights)))
      Sl.sf1 <- Sl.setup(gam1ff)
      rm(gam1ff)
    }
  }
  if (family == "TW") {
    gam1TW <- eval(substitute(gam(formulamgcv, gamma = infl.fac, weights = weights, data = data, offset = offset, knots = knots, drop.unused.levels = drop.unused.levels, family = twlss()), list(weights = weights, offset = offset)))
    if(!is.null(gam1TW$offset[[1]])) offset <- gam1TW$offset[[1]] # this is not ok in mgcv
  }
  if (surv == TRUE && family %in% bl) {
    surv.flex <- TRUE
    f.eq1 <- form.eq12R$f.eq1
    data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
    if (type.cens == "R"){
    
      if(!is.numeric(data[,"(cens)"])) stop("Check that the correct option for argument type.cens has been used.") 
    
      tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), data = data, weights = cens, drop.unused.levels = drop.unused.levels), list(cens = cens)))
      
      }
    if (type.cens %in% c("L", "I")) {
      f.eq1LI <- form.eq12R$f.eq1LI
      data$temp.respV <- data[, v1[1]]
      if (type.cens == "L") data$temp.respV[cens == 0] <- data$temp.respV[cens == 0]/2
      if (type.cens == "I") {
        r.op <- range(data$temp.respV[cens == 1])
        r.upL <- length(data$temp.respV[cens == 0])
        for (i in 1:r.upL) {
          if (data$temp.respV[cens == 0][i] < r.op[1]) {
            r.op <- range(c(r.op, data$temp.respV[cens == 0][i]))
          }
          if (data[cens == 0, v.rB][i] > r.op[2]) {
            data$temp.respV[cens == 0][i] <- data[cens == 0, v.rB][i]
            r.op <- range(c(r.op, data$temp.respV[cens == 0][i]))
          }
          if (data$temp.respV[cens == 0][i] > r.op[1] && data[cens == 0, v.rB][i] < r.op[2]) data$temp.respV[cens == 0][i] <- (data$temp.respV[cens == 0][i] + data[cens == 0, v.rB][i])/2
        }
      }
      tempb <- eval(substitute(gam(f.eq1LI, family = cox.ph(), data = data, weights = rep(1, length(cens)), 
                                   drop.unused.levels = drop.unused.levels), list(cens = cens)))
    }
    if (type.cens %in% c("mixed")) {
      f.eq1LI <- form.eq12R$f.eq1LI
      data$temp.respV <- data[, v1[1]]
      if (any(unique(cens) %in% c("L", "LT"))) 
        data$temp.respV[cens %in% c("L", "LT")] <- data$temp.respV[cens %in% c("L", "LT")]/2
      if (any(unique(cens) %in% c("I", "IT"))) {
        if( any(!(cens %in% c("I", "IT"))) ) r.op  <- range(data$temp.respV[!(cens %in% c("I", "IT"))])
        else r.op <- range(data$temp.respV, data[, v.rB] ) # here we would select the rows such that cens %in% c("I", "IT") but this is redundant here as it is the whole dataset 
        r.upL <- length(data$temp.respV[cens %in% c("I","IT")])
        for (i in 1:r.upL) {
          if (data$temp.respV[cens %in% c("I", "IT")][i] < r.op[1]) r.op <- range(c(r.op, data$temp.respV[cens %in% c("I", "IT")][i]))
          if (data[cens %in% c("I", "IT"), v.rB][i] > r.op[2]) {
            data$temp.respV[cens %in% c("I", "IT")][i] <- data[cens %in% c("I", "IT"), v.rB][i]
            r.op <- range(c(r.op, data$temp.respV[cens %in% c("I", "IT")][i]))
          }
          if (data$temp.respV[cens %in% c("I", "IT")][i] > r.op[1] && data[cens %in% c("I", "IT"), v.rB][i] < r.op[2]) data$temp.respV[cens %in% c("I", "IT")][i] <- (data$temp.respV[cens %in% c("I", "IT")][i] + data[cens %in% c("I", "IT"), v.rB][i])/2
        }
      }
      tempb <- eval(substitute(gam(f.eq1LI, family = cox.ph(), data = data, weights = rep(1, length(cens)), drop.unused.levels = drop.unused.levels), list(cens = cens)))
    }
    data$Sh <- as.vector(mm(predict(tempb, type = "response"), min.pr = min.pr, max.pr = max.pr))
    if (type.cens %in% c("mixed")) {
      cens1 <- cens
      cens1 <- ifelse(cens1 %in% c("L", "I", "R", "LT","IT", "RT"), 0, 1)
      cens1 <- ifelse(cens1 == 0, 1e-07, cens1)
    }
    else cens1 <- ifelse(cens == 0, 1e-07, cens)
    if (type.cens != "R") {
      data[, v1[1]] <- data$temp.respV
      attr(data, "terms") <- NULL
    }
    gam1 <- eval(substitute(scam(formula.eq1, gamma = infl.fac, weights = weights * cens1, data = data), list(weights = weights, cens1 = cens1)))
    lsgam1 <- length(gam1$smooth)
    if (lsgam1 == 0) stop("You must at least use a monotonic smooth function of time.")
    clsm <- ggr <- NA
    for (i in 1:lsgam1) clsm[i] <- class(gam1$smooth[[i]])[1]
    
    if(sum(as.numeric(clsm %in% c("mpi.smooth"))) == 0) stop("You must have a monotonic smooth of time, mpi.")
    pos.mpi <- which(clsm == "mpi.smooth")

    l.sp1 <- length(gam1$sp)
    if (l.sp1 != 0) sp1 <- gam1$sp
    sp1[pos.mpi] <- 1

    gam.call <- gam1$call
    gam.call$sp <- sp1
    gam1 <- eval(gam.call)

    j <- 1
    for (i in 1:lsgam1) if (max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$term))) != 0 && clsm[i] == "mpi.smooth") mono.sm.pos <- c(mono.sm.pos, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para))
    if (type.cens == "mixed") {
      indvU <- indvR <- indvL <- indvI <- rep(0, n)
      if (any(unique(cens) == "U")) indvU[cens == "U"] <- 1
      if (any(unique(cens) == "R")) indvR[cens == "R"] <- 1
      if (any(unique(cens) == "L")) indvL[cens == "L"] <- 1
      if (any(unique(cens) == "I")) indvI[cens == "I"] <- 1
      indvUT <- indvRT <- indvLT <- indvIT <- rep(0, n)
      if (any(unique(cens) == "UT")) indvUT[cens == "UT"] <- 1
      if (any(unique(cens) == "RT")) indvRT[cens == "RT"] <- 1
      if (any(unique(cens) == "LT")) indvLT[cens == "LT"] <- 1
      if (any(unique(cens) == "IT")) indvIT[cens == "IT"] <- 1
    }
    if (type.cens == "R") X1 <- predict(gam1, type = "lpmatrix")
    if (type.cens %in% c("I", "L", "mixed")) {
      data[, v1[1]] <- u1resp
      X1 <- predict(gam1, type = "lpmatrix", newdata = data)
      if (type.cens == "I") {
        data[cens == 0, v1[1]] <- data[cens == 0, v.rB]
        X2 <- predict(gam1, type = "lpmatrix", newdata = data)
      }
      if (type.cens == "mixed") {
        if (any(unique(cens) %in% c("I", "IT"))) {
          data[cens %in% c("I", "IT"), v1[1]] <- data[cens %in% c("I", "IT"), v.rB]
          X2 <- predict(gam1, type = "lpmatrix", newdata = data)
        }
        if (any(unique(cens) %in% c("UT", "LT", "RT","IT"))) {
          data[cens %in% c("UT", "LT", "RT", "IT"), v1[1]] <- data[cens %in% c("UT", "LT", "RT","IT"), v.td]
          X3 <- predict(gam1, type = "lpmatrix", newdata = data)
        }
      }
    }
    if (is.null(X2)) X2 <- matrix(1, dim(X1)[1], dim(X1)[2])
    if (is.null(X3)) X3 <- matrix(1, dim(X1)[1], dim(X1)[2])
    #if (!is.null(indexT) && k.tvc != 0) if (range(X1[, indexT])[1] < 0) stop("Check design matrix for smooth(s) of tvc terms.")
    data[, v1[1]] <- u1resp
    Xd <- Xdpred(gam1, data, v1[1])
    start.v1 <- c(gam1$coefficients)
    gam1$y <- data[, v1[1]]
    if (!is.null(indexT)) {
      start.v2 <- start.v1
      start.v2[mono.sm.pos] <- exp(start.v2[mono.sm.pos])
      while (range(Xd %*% start.v2)[1] < 0) start.v2[indexT] <- 0.999 * start.v2[indexT]
      start.v1[indexT] <- start.v2[indexT]
      gam1$coefficients <- gam1$coefficients.t <- start.v1
      gam1$coefficients.t[mono.sm.pos] <- exp(gam1$coefficients.t[mono.sm.pos])
    }
    if (!is.null(hrate) || !is.null(d.lchrate) || !is.null(d.rchrate) || !is.null(d.lchrate.td) || !is.null(d.rchrate.td)) {
      resExcInd <- survExcInd(n = n, cens = cens, type.cens = type.cens, hrate = hrate, d.lchrate = d.lchrate, d.rchrate = d.rchrate, d.lchrate.td = d.lchrate.td, d.rchrate.td = d.rchrate.td)
      hrate <- resExcInd$hrate
      d.lchrate <- resExcInd$d.lchrate
      d.rchrate <- resExcInd$d.rchrate
      d.lchrate.td <- resExcInd$d.lchrate.td
      d.rchrate.td <- resExcInd$d.rchrate.td
    }
    if (informative == "yes") {
      f.eq1 <- form.eq12R2$f.eq1
      data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
      tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), data = data, weights = 1 - cens, drop.unused.levels = drop.unused.levels), list(cens = cens)))
      data$Sh <- as.vector(mm(predict(tempb, type = "response"),min.pr = min.pr, max.pr = max.pr))
      cens1 <- ifelse((1 - cens) == 0, 1e-07, 1 - cens)
      gam2 <- eval(substitute(scam(formula.eq2, gamma = infl.fac, weights = weights * cens1, data = data), list(weights = weights, cens1 = cens1)))
      lsgam2 <- length(gam2$smooth)
      if (lsgam2 == 0) stop("You must use at least a monotonic smooth function of time.")
      clsm <- ggr <- NA
      for (i in 1:lsgam2) clsm[i] <- class(gam2$smooth[[i]])[1]
      if(sum(as.numeric(clsm %in% c("mpi.smooth"))) == 0) stop("You must have a monotonic smooth of time, mpi, in the informative equation.")
      pos.mpi <- which(clsm == "mpi.smooth")

      if (v1[1] %in% inform.cov) stop("Time can not be an informative covariate.")
      l.sp2 <- length(gam2$sp)
      if (l.sp2 != 0) sp2 <- gam2$sp
      sp2[pos.mpi] <- 1
      gam.call <- gam2$call
      gam.call$sp <- sp2
      gam2 <- eval(gam.call)
      j <- 1
      for (i in 1:lsgam2) if (max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) != 0 && clsm[i] == "mpi.smooth") mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para))
      X2 <- predict(gam2, type = "lpmatrix")
      #if (!is.null(indexT) && k.tvc != 0) if (range(X2[, indexT])[1] < 0) stop("Check design matrix for smooth(s) of tvc terms.")
      Xd2 <- Xdpred(gam2, data, v2[1])
      start.v2 <- c(gam2$coefficients)
      gam2$y <- data[, v2[1]]
    }
  }
  gam1$formula <- formula.eq1r
  lsgam1 <- length(gam1$smooth)
  y1 <- y1.test
  if (family %in% c("LN")) y1 <- log(y1)
  attr(data, "terms") <- NULL
  if (!(surv == TRUE && family %in% bl)) {
    names(gam1$model)[1] <- as.character(formula.eq1r[2])
    X1 <- predict(gam1, type = "lpmatrix")
    l.sp1 <- length(gam1$sp)
    sp1 <- gam1$sp
  }
  gp1 <- gam1$nsdf
  X1.d2 <- dim(X1)[2]
  if (surv == TRUE && family %in% bl && informative == "yes") {
    gam2$formula <- formula.eq2r
    lsgam2 <- length(gam2$smooth)
    y2 <- y2.test
    gp2 <- gam2$nsdf
    X2.d2 <- dim(X2)[2]
  }
  if (family != "TW") {
    log.nu.1 <- log.sig2.1 <- NULL
    if (!(family %in% c(m1d, bl))) {
      start.snR <- startsn(family, y1, left.trunc = left.trunc)
      log.sig2.1 <- start.snR$log.sig2.1
      names(log.sig2.1) <- "sigma.star"
      if (family %in% c("GP", "GPII", "GPo", "DGP", "DGPII")) log.sig2.1 <- esres[2]
      if (family %in% c(m3)) {
        log.nu.1 <- start.snR$log.nu.1
        names(log.nu.1) <- "nu.star"
      }
    }
  }
  if (surv == TRUE && family2 %in% bl && informative == "yes") {
    infsetupR <- inform.setup(gam1, gam2, inform.cov, start.v1, start.v2, lsgam1, lsgam2)
    start.v1 <- infsetupR$start.v1
    start.v2 <- infsetupR$start.v2
    Gmat12 <- matrix(NA, n, length(start.v1))
    Hmat12 <- matrix(NA, length(c(start.v1, start.v2)), length(c(start.v1, start.v2)))
    if (!is.null(infsetupR$par.pos1) && is.null(infsetupR$smo.pos1)) {
      Xi <- Xi <- X1[, c(infsetupR$par.pos1)]
      X1ni <- X1[, -c(infsetupR$par.pos1)]
      X2ni <- X2[, -c(infsetupR$par.pos2)]
      inde.inf1 <- c(infsetupR$par.pos1)
      inde.inf2 <- c(infsetupR$par.pos2)
    }
    if (!is.null(infsetupR$par.pos1) && !is.null(infsetupR$smo.pos1)) {
      Xi <- Xi <- X1[, c(infsetupR$par.pos1, infsetupR$smo.pos1)]
      X1ni <- X1[, -c(infsetupR$par.pos1, infsetupR$smo.pos1)]
      X2ni <- X2[, -c(infsetupR$par.pos2, infsetupR$smo.pos2)]
      inde.inf1 <- c(infsetupR$par.pos1, infsetupR$smo.pos1)
      inde.inf2 <- c(infsetupR$par.pos2, infsetupR$smo.pos2)
    }
    if (is.null(infsetupR$par.pos1) && !is.null(infsetupR$smo.pos1)) {
      Xi <- Xi <- X1[, c(infsetupR$smo.pos1)]
      X1ni <- X1[, -c(infsetupR$smo.pos1)]
      X2ni <- X2[, -c(infsetupR$smo.pos2)]
      inde.inf1 <- c(infsetupR$smo.pos1)
      inde.inf2 <- c(infsetupR$smo.pos2)
    }
  }
  if (family != "TW") {
    if (surv == TRUE && family %in% bl && informative == "yes") start.v1 <- c(start.v1, start.v2)
    else {
      if (family %in% c(m1d)) start.v1 <- c(gam1$coefficients)
      if (family %in% c(m2, m2d)) start.v1 <- c(gam1$coefficients, log.sig2.1)
      if (family %in% c(m3, m3d)) start.v1 <- c(gam1$coefficients, log.sig2.1,log.nu.1)
    }
  }
  if (family == "TW") log.nu.1 <- log.sig2.1 <- 0.1
  if (l.flist > 1 && !(surv == TRUE && family %in% bl)) {
    vo <- list(log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, n = n, drop.unused.levels = drop.unused.levels)
    overall.svGR <- overall.svG(formula, data, ngc = 2, family, M, vo, gam1, gam2, type = "gaml", knots = knots)
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
    if (M$sp.method != "perf") {
      Sl.sf2 <- overall.svGR$Sl.sf2
      Sl.sf3 <- overall.svGR$Sl.sf3
    }
    if (family != "TW") {
      start.v1 <- overall.svGR$start.v
      sp2 <- overall.svGR$sp2
      sp3 <- overall.svGR$sp3
    }
    if (family == "TW") {
      spmgcv <- gam1TW$sp
      l.sp1mgcv <- l.sp1
      l.sp2mgcv <- l.sp3
      l.sp3mgcv <- l.sp2
      if (l.sp1 != 0) {
        sp1 <- spmgcv[1:l.sp1mgcv]
        names(sp1) <- names(gam1$sp)
      }
      if (l.sp2 != 0) {
        sp2 <- spmgcv[(l.sp1mgcv + l.sp2mgcv + 1):(l.sp1mgcv + l.sp2mgcv + l.sp3mgcv)]
        names(sp2) <- names(gam2$sp)
      }
      if (l.sp3 != 0) {
        sp3 <- spmgcv[(l.sp1mgcv + 1):(l.sp1mgcv + l.sp2mgcv)]
        names(sp3) <- names(gam3$sp)
      }
      start.v1TW <- start.v1 <- coef(gam1TW)
      X1.d2mgcv <- X1.d2
      X2.d2mgcv <- X3.d2
      X3.d2mgcv <- X2.d2
      start.v1TW[(X1.d2 + 1):(X1.d2 + X2.d2)] <- start.v1[(X1.d2mgcv + X2.d2mgcv + 1):(X1.d2mgcv + X2.d2mgcv + X3.d2mgcv)]
      start.v1TW[(X1.d2 + X2.d2 + 1):(X1.d2 + X2.d2 + X3.d2)] <- start.v1[(X1.d2mgcv + 1):(X1.d2mgcv + X2.d2mgcv)]
      start.v1 <- start.v1TW
      names(start.v1) <- names(overall.svGR$start.v)
    }
  }
  if (surv == TRUE && family %in% bl && informative == "yes" && l.sp2 > 1 && !is.null(infsetupR$scv) && !is.null(infsetupR$inds2)) sp2 <- sp2[-c(infsetupR$inds2)]
  if (surv == TRUE && family %in% bl && informative == "yes") {
    pic <- NA
    av <- all.vars(formula[[1]])
    for (i in 1:length(inform.cov)) pic[i] <- which(inform.cov[i] == av)
    tfor <- formula(drop.terms(terms(formula[[1]]), pic, keep.response = TRUE))
    fgam <- gam(tfor, data = data, fit = FALSE, drop.unused.levels = drop.unused.levels)
    pfgam <- NA
    for (i in 1:length(fgam$smooth)) pfgam[i] <- fgam$smooth[[i]]$first.para
  }
  spgamlss1 <- c(sp1, sp2, sp3)
  GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
              gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, 
              gam9 = gam9, K1 = NULL)
  if (l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0) {
    L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), 
                  l.gam3 = length(gam3$coefficients), l.gam4 = 0, 
                  l.gam5 = 0, l.gam6 = 0, l.gam7 = 0, l.gam8 = 0, 
                  l.gam9 = 0)
    L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                 l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0, 
                 l.sp9 = 0)
    qu.mag1 <- S.m(GAM, L.SP, L.GAM)
    if (surv == TRUE && family %in% bl && informative == "yes" && l.sp2 > 1 && !is.null(infsetupR$scv)) {
      qu.mag1$Ss <- qu.mag1$Ss[-c(l.sp1 + infsetupR$inds2)]
      qu.mag1$rank <- qu.mag1$rank[-c(l.sp1 + infsetupR$inds2)]
      qu.mag1$off <- qu.mag1$off[-c(l.sp1 + infsetupR$inds2)]
    }
    if (surv == TRUE && family %in% bl && informative == "yes") {
      qu.mag1$off[(l.sp1 + 1):length(qu.mag1$off)] <- pfgam + X1.d2
      test.sv.inf <- start.v1[qu.mag1$off]
    }
  }
  if (missing(parscale)) parscale <- 1
  respvec2 <- list(y1 = y1, univ = 2)
  lsgam2 <- length(gam2$smooth)
  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)
  lsgam9 <- length(gam9$smooth)
  if (robust == TRUE && family %in% c(m1d, m2d)) {
    eta.m <- max(predict(gam1, type = "link"))
    if (family %in% c(m2d)) sigma2.m <- exp(log.sig2.1)
    else sigma2.m <- 1
    if (!(family %in% c("tP","tNBI","tNBII","tPIG"))) ygrid <- 0:(max(y1) * 100)
    else ygrid <- 1:(max(y1) * 100)
    pdf.test <- distrHsATDiscr(ygrid, eta.m, sigma2.m, 1, 
                               family, y1m, robust = FALSE, min.dn = ygrid.tol, 
                               min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)$pdf2 > ygrid.tol
    if (any(pdf.test == FALSE) != TRUE) {
      if (!(family %in% c("tP","tNBI","tNBII","tPIG"))) ygrid <- 0:(max(y1) * 200)
      else ygrid <- 1:(max(y1) * 200)
      pdf.test <- distrHsATDiscr(ygrid, eta.m, sigma2.m, 
                                 1, family, y1m, robust = FALSE, min.dn = ygrid.tol, 
                                 min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)$pdf2 > ygrid.tol
    }
    maxv <- max(which(pdf.test == TRUE))
    ygrid <- ygrid[1:maxv]
  }
  else ygrid <- NULL
  my.env <- new.env()
  my.env$k <- k.tvc
  my.env$indN <- NULL
  my.env$V <- NULL
  if (sp.method == "efs") {
    Sl.sf <- list()
    LSl.sf1 <- length(Sl.sf1)
    LSl.sf2 <- length(Sl.sf2)
    LSl.sf3 <- length(Sl.sf3)
    if (LSl.sf1 == 0) pa1 <- matrix(0, gp1, gp1)
    if (LSl.sf2 == 0) pa2 <- matrix(0, gp2, gp2)
    j <- 1
    if (LSl.sf1 != 0) {
      for (i in 1:LSl.sf1) {
        Sl.sf[[j]] <- Sl.sf1[[i]]
        j <- j + 1
      }
      attr(Sl.sf, "lambda") <- attr(Sl.sf1, "lambda")
      attr(Sl.sf, "E") <- attr(Sl.sf1, "E")
      attr(Sl.sf, "cholesky") <- attr(Sl.sf1, "cholesky")
    }
    if (LSl.sf2 != 0) {
      for (i in 1:LSl.sf2) {
        Sl.sf[[j]] <- Sl.sf2[[i]]
        Sl.sf[[j]]$start <- Sl.sf[[j]]$start + length(gam1$coefficients)
        Sl.sf[[j]]$stop <- Sl.sf[[j]]$stop + length(gam1$coefficients)
        j <- j + 1
      }
      attr(Sl.sf, "lambda") <- c(attr(Sl.sf, "lambda"), attr(Sl.sf2, "lambda"))
      attr(Sl.sf, "cholesky") <- attr(Sl.sf2, "cholesky")
      if (LSl.sf1 != 0) attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), attr(Sl.sf2, "E"))
      if (LSl.sf1 == 0) attr(Sl.sf, "E") <- adiag(pa1, attr(Sl.sf2, "E"))
    }
    if (LSl.sf3 != 0) {
      for (i in 1:LSl.sf3) {
        Sl.sf[[j]] <- Sl.sf3[[i]]
        Sl.sf[[j]]$start <- Sl.sf[[j]]$start + length(gam1$coefficients) + length(gam2$coefficients)
        Sl.sf[[j]]$stop <- Sl.sf[[j]]$stop + length(gam1$coefficients) + length(gam2$coefficients)
        j <- j + 1
      }
      attr(Sl.sf, "lambda") <- c(attr(Sl.sf, "lambda"), attr(Sl.sf3, "lambda"))
      attr(Sl.sf, "cholesky") <- attr(Sl.sf3, "cholesky")
      if (LSl.sf1 != 0 && LSl.sf2 != 0) attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), attr(Sl.sf3, "E"))
      if (LSl.sf1 == 0 && LSl.sf2 != 0) attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), attr(Sl.sf3, "E"))
      if (LSl.sf1 == 0 && LSl.sf2 == 0) attr(Sl.sf, "E") <- adiag(pa1, pa2, attr(Sl.sf3, "E"))
      if (LSl.sf1 != 0 && LSl.sf2 == 0) attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), pa2, attr(Sl.sf3, "E"))
    }
  }
  VC <- list(lsgam1 = lsgam1, ygrid = ygrid, lsgam2 = lsgam2, end.surv = FALSE,
             indexT = indexT, D = D, my.env = my.env, k = k.tvc, 
             pos.pb = pos.pb, lsgam3 = lsgam3, r.type = r.type, Sl.sf = Sl.sf, 
             sp.method = sp.method, lsgam4 = lsgam4, gam1TW = gam1TW, 
             lsgam5 = lsgam5, lsgam6 = lsgam6, lsgam7 = lsgam7, lsgam8 = lsgam8, 
             lsgam9 = lsgam9, fgam = fgam, ad.ind = FALSE, X1 = X1, 
             X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, 
             X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
             X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, 
             X8.d2 = X8.d2, gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, 
             gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, l.sp1 = l.sp1, 
             l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, l.sp5 = l.sp5, 
             l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = 0, 
             infl.fac = infl.fac, weights = weights, offset = offset, fp = fp, hess = NULL, 
             Model = "CC", univ.gamls = TRUE, gc.l = gc.l, n = n, 
             extra.regI = extra.regI, parscale = parscale, margins = c(family, family),
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, 
             m1d = m1d, m2d = m2d, m3d = m3d, bl = bl, triv = FALSE, 
             y1m = y1m, y2m = y2m, robust = robust, rc = rc, cens = cens, 
             surv = surv, lB = lB, uB = uB, gev.par = gev.par, chunk.size = chunk.size, 
             Xd1 = Xd, Xd2 = Xd2, mono.sm.pos = mono.sm.pos, mono.sm.pos2 = mono.sm.pos2, 
             surv.flex = surv.flex, informative = informative, inform.cov = inform.cov, 
             infsetupR = infsetupR, gp2.inf = NULL, Xi = Xi, X1ni = X1ni, 
             X2ni = X2ni, inde.inf2 = inde.inf2, inde.inf1 = inde.inf1, 
             Gmat12 = Gmat12, Hmat12 = Hmat12, v1pred = v1[1], sp.fixed = sp.fixed, 
             indvU = indvU, indvR = indvR, indvL = indvL, indvI = indvI, 
             indvUT = indvUT, indvRT = indvRT, indvLT = indvLT, indvIT = indvIT, 
             hrate = hrate, d.lchrate = d.lchrate, d.rchrate = d.rchrate, 
             d.lchrate.td = d.lchrate.td, d.rchrate.td = d.rchrate.td, 
             zero.tol = 0.01, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, mar1surv = mar1surv, mar2surv = mar2surv,
             end.surv = FALSE, left.trunc = left.trunc)
             
  if (surv == TRUE && family2 %in% bl && informative == "yes") {
    if (!is.null(infsetupR$pcv)) VC$gp2.inf <- gp2 - length(infsetupR$par.pos1)
    else VC$gp2.inf <- gp2
    if (l.sp2 > 1 && !is.null(infsetupR$scv)) VC$l.sp2 <- l.sp2 - length(infsetupR$inds2)
    VC$margins <- c(family, family2)
  }
  if (gc.l == TRUE) gc()
  if (family != "GEVlink") func.opt1 <- chooseFuncOpt1(family, m1d, m2d, m2, m3, bl, informative, type.cens, hrate, truncation.time)
  if (family == "GEVlink") func.opt1 <- bprobgHsContUnivBIN
  y1.m <- y1
  if (family == "LN") y1.m <- exp(y1)
  VC$y1 <- y1.m
  SemiParFit <- SemiParBIV.fit(func.opt = func.opt1, start.v = start.v1, 
                               rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, 
                               tolsp = tolsp, respvec = respvec2, VC = VC, sp = spgamlss1, qu.mag = qu.mag1)
  SemiParFit$robust <- robust
  SemiParFit.p <- gamlss.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)
  SemiParFit <- SemiParFit.p$SemiParFit
  if (gc.l == TRUE) gc()
  cov.c(SemiParFit)
  gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data
  rm(gam1TW)
  L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, formula = formula, 
            robust = robust, edf11 = SemiParFit.p$edf11, gam1 = gam1, 
            gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, 
            gam6 = gam6, gam7 = gam7, gam8 = gam8, coefficients = SemiParFit$fit$argument, 
            iterlimsp = iterlimsp, weights = weights, offset = offset, cens = cens, 
            sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
            l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
            l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, 
            l.sp9 = l.sp9, gam9 = gam9, fp = fp, iter.if = SemiParFit$iter.if, 
            iter.inner = SemiParFit$iter.inner, sigma2 = SemiParFit.p$sigma2, 
            sigma = SemiParFit.p$sigma, sigma2.a = SemiParFit.p$sigma2.a, 
            sigma.a = SemiParFit.p$sigma.a, nu = SemiParFit.p$nu, 
            nu.a = SemiParFit.p$nu.a, X1 = X1, X2 = X2, X3 = X3, 
            X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, 
            X2.d2 = X2.d2, X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, 
            X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, He = SemiParFit.p$He, 
            HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Vb1 = SemiParFit.p$Vb1, 
            Ve = SemiParFit.p$Ve, F = SemiParFit.p$F, F1 = SemiParFit.p$F1, 
            Vb.t = SemiParFit.p$Vb.t, coef.t = SemiParFit.p$coef.t, 
            t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
            edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, 
            edf3 = SemiParFit.p$edf3, edf4 = SemiParFit.p$edf4, 
            edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, 
            edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8, 
            edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, 
            edf1.3 = SemiParFit.p$edf1.3, edf1.4 = SemiParFit.p$edf1.4, 
            edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
            edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, 
            R = SemiParFit.p$R, bs.mgfit = SemiParFit$bs.mgfit, 
            conv.sp = SemiParFit$conv.sp, wor.c = SemiParFit$wor.c, 
            eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$etas1, 
            eta3 = SemiParFit$fit$etan1, y1 = y1.m, margins = c(family, family),
            logLik = SemiParFit.p$logLik, hess = TRUE, 
            qu.mag = qu.mag1, gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, 
            gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, VC = VC, 
            magpp = SemiParFit$magpp, type.cens = type.cens, Cont = "YES", 
            l.flist = l.flist, triv = FALSE, univar.gamlss = TRUE, 
            call = cl, gev.par = gev.par, ygrid = ygrid, r.weights = SemiParFit$fit$d.psi, 
            surv = surv, surv.flex = surv.flex, test.sv.inf = test.sv.inf, 
            rangeSurv = rangeSurv, Model = "CC", end.surv = FALSE, mar1surv = mar1surv, mar2surv = mar2surv, end.surv = FALSE,
            left.trunc = left.trunc)
  class(L) <- c("gamlss", "SemiParBIV", "gjrm")
  L
}
