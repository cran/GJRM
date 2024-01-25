hazsurv <- function (x, eq, newdata, type = "surv", t.range = NULL, t.vec = NULL, 
          intervals = TRUE, n.sim = 100, prob.lev = 0.05, shade = FALSE, 
          bars = FALSE, ylim, ylab, xlab, pch, ls = 100, baseline = FALSE, 
          min.dn = 1e-200, min.pr = 1e-200, max.pr = 1, plot.out = TRUE, 
          print.progress = TRUE, ...) 
{
  pr <- h <- hs <- prs <- CIpr <- CIh <- poe <- poet <- NULL
  pr.avg <- h.avg <- ch.avg <- CIpr.avg <- CIh.avg <- CIch.avg <- NULL
  toleps <- 1e-04
  if (x$univar.gamlss == FALSE && x$surv.flex == TRUE && x$margins[1] %in% 
      c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$bl)) 
    eq <- 2
  if (missing(eq) && x$univar.gamlss == FALSE) 
    stop("You must provide the equation number (either 1 or 2).")
  if (x$univar.gamlss == TRUE) 
    eq <- 1
  if (missing(newdata)) 
    stop("You have to provide a new data frame.")
  if (!is.data.frame(newdata)) 
    stop("You have to provide a new data frame.")
  if (!(type %in% c("surv", "hazard", "cumhaz"))) 
    stop("The type argument can either be surv, hazard or cumhaz")
  if (x$surv.flex == FALSE) 
    stop("This function is only suitable for flexible survival models.")
  if (missing(ylim)) 
    ylim <- NULL
  if (nrow(newdata) < 1) 
    stop("The data frame needs to have at least one row.")
  if (length(t.range) > 2) 
    stop("When using t.range only provide min and max of the time range to be considered.")
  if (!is.null(t.range) & !is.null(t.vec)) 
    stop("You cannot provide both t.range and t.vec. See help for more details.")
  if (shade == TRUE & bars == TRUE) 
    stop("shade and bars are mutually exclusive, they cannot both be TRUE.")
  
  if(eq == 1){
    ntv  <- as.character(x$formula[[1]][2])
    
    #rlb <- range(x$y1)[1]
    
    if(!is.null(t.range)){ 
      rlb <- t.range[1]
    } else if(!is.null(t.vec)){
      rlb <- t.vec #check on rlb works anyway 
    } else{ 
      if(x$univar.gamlss == TRUE){
        rlb <- x$rangeSurv[1]
      } else if(x$gamlssfit == TRUE){
        rlb = x$gamlss1$rangeSurv[1]
      } else {
        rlb <- range(x$y1)[1]
      } 
    }
    
    rlb <- ifelse(rlb < toleps, toleps, rlb)
    
    if(!is.null(t.range)){
      tv <- seq(rlb, t.range[2], length.out = ls)
    } else if(!is.null(t.vec)){
      tv <- rlb
    } else{
      if(x$univar.gamlss == TRUE){
        tv <- seq(rlb, x$rangeSurv[2], length.out = ls)
      } else if(x$gamlssfit == TRUE) {
        tv <- seq(rlb, x$gamlss1$rangeSurv[2], length.out = ls)
      } else {
        tv <- seq(rlb, range(x$y1)[2], length.out = ls)
      } 
    }
    
    
    indp <- 1:x$VC$X1.d2
    gob  <- x$gam1
  }
  
  
  
  if(eq == 2){
    ntv  <- as.character(x$formula[[2]][2])
    
    if(!is.null(t.range)){ 
      rlb <- t.range[1]
    } else if(!is.null(t.vec)){
      rlb <- t.vec
    } else { 
      if(x$gamlssfit == TRUE) rlb <- x$gamlss2$rangeSurv[1] else rlb <- range(x$y2)[1]
    }  
    
    rlb <- ifelse(rlb < toleps, toleps, rlb)
    
    if(!is.null(t.range)){
      tv <- seq(rlb, t.range[2], length.out = ls)
    } else if(!is.null(t.vec)){
      tv <- rlb
    } else{
      if(x$gamlssfit == TRUE) tv = seq(rlb, x$gamlss2$rangeSurv[2], length.out = ls) else tv   <- seq(rlb, range(x$y2)[2], length.out = ls)
    }   
    
    indp <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
    gob  <- x$gam2
  }
  
  ti <- data.frame(tv)
  names(ti) <- ntv
  pr.cumul <- h.cumul <- ch.cumul <- rep(0, length(tv))
  if (intervals == TRUE) 
    CIpr.cumul <- CIh.cumul <- CIch.cumul <- matrix(0, length(tv), 
                                                    n.sim)
  if (!is.null(x$VC$mono.sm.pos)) 
    mono.sm.pos <- x$VC$mono.sm.pos
  else mono.sm.pos <- c(x$VC$mono.sm.pos1, x$VC$mono.sm.pos2 + 
                          x$VC$X1.d2)
  if (intervals == TRUE) {
    bs <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
    bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, 
                                bs[, mono.sm.pos])
  }
  for (obs in 1:dim(newdata)[1]) {
    if (obs%%30 == 0 & print.progress) 
      print(paste(round(obs/dim(newdata)[1] * 100), "% of iterations complete", 
                  sep = ""))
    newdata.temp = as.data.frame(newdata[obs, ])
    row.names(newdata.temp) = NULL
    names(newdata.temp) <- names(newdata)
    newdata.temp <- data.frame(ti, newdata.temp)
    Xpred <- predict(x, newdata.temp, eq = eq, type = "lpmatrix")
    if (baseline == TRUE) {
      Xd <- Xdpred(gob, newdata.temp, ntv)
      ind0 <- (colSums(Xd == 0) == dim(Xpred)[1])
      ind0[1] <- FALSE
      Xpred[, ind0] <- 0
    }
    params1 <- x$coef.t[indp]
    eta1 <- Xpred %*% params1
    pd <- probmS(eta1, x$VC$margins[eq], min.dn = min.dn, 
                 min.pr = min.pr, max.pr = max.pr)
    pr <- pd$pr
    if (intervals == TRUE) {
      eta1s <- Xpred %*% t(bs[, indp])
      pds <- probmS(eta1s, x$VC$margins[eq], min.dn = min.dn, 
                    min.pr = min.pr, max.pr = max.pr)
      prs <- pds$pr
    }
    if (type == "surv") {
      if (intervals == TRUE) {
        for (i in 1:length(tv)) {
          poe <- which(prs[i, ] %in% boxplot.stats(prs[i, 
                                                       ])$out)
          prs[, poe] <- NA
          poe <- union(poe, poet)
          poet <- poe
        }
        CIpr <- prs
      }
      pr.cumul = pr.cumul + pr
      if (intervals == TRUE) 
        CIpr.cumul = CIpr.cumul + CIpr
      if (obs == dim(newdata)[1]) {
        pr.avg = pr.cumul/dim(newdata)[1]
        if (intervals == TRUE) {
          CIpr.avg <- CIpr.cumul/dim(newdata)[1]
          CIpr.avg <- rowQuantiles(CIpr.avg, probs = c(prob.lev/2, 
                                                       1 - prob.lev/2), na.rm = TRUE)
        }
        if (plot.out) {
          if (intervals == TRUE) 
            CIpr.avg.temp = matrix(CIpr.avg, nrow = length(tv))
          if (is.null(ylim) && intervals == TRUE) 
            ylim <- c(min(CIpr.avg.temp[, 1]), max(CIpr.avg.temp[, 
                                                                 2]))
          if (is.null(ylim) && intervals == FALSE) 
            ylim <- c(min(pr.avg), max(pr.avg))
          if (missing(ylab)) 
            ylab <- "Survival function"
          if (missing(xlab)) 
            xlab <- "Time"
          if (missing(pch)) 
            pch <- 19
          if (length(tv) == 1 | bars == TRUE) 
            type.plot = "p"
          else type.plot = "l"
          plot(tv, pr.avg, type = type.plot, ylab = ylab, 
               xlab = xlab, ylim = ylim, pch = pch, ...)
          if (intervals == TRUE) {
            if (length(tv) == 1 | bars == TRUE) {
              for (i in 1:length(tv)) lines(c(tv[i], 
                                              tv[i]), CIpr.avg.temp[i, ])
            }
            else {
              if (shade == FALSE) {
                lines(tv, CIpr.avg.temp[, 1], lty = 2)
                lines(tv, CIpr.avg.temp[, 2], lty = 2)
              }
              else {
                polygon(c(tv, rev(tv)), c(CIpr.avg.temp[, 
                                                        1], rev(CIpr.avg.temp[, 2])), col = "gray80", 
                        border = NA)
                lines(tv, pr.avg, type = "l")
              }
            }
          }
        }
      }
    }
    if (type == "hazard") {
      if (baseline == FALSE) 
        Xd <- Xdpred(gob, newdata.temp, ntv)
      Xthe <- Xd %*% params1
      Gp <- pd$dS
      h <- -Gp/pr * Xthe
      if (intervals == TRUE) {
        Gps <- pds$dS
        Xthes <- Xd %*% t(bs[, indp])
        hs <- -Gps/prs * Xthes
        for (i in 1:length(tv)) {
          poe <- which(hs[i, ] %in% boxplot.stats(hs[i, 
                                                     ])$out)
          hs[, poe] <- NA
          poe <- union(poe, poet)
          poet <- poe
        }
        CIh <- hs
        CIh <- ifelse(CIh < 0, 0, CIh)
      }
      h.cumul = h.cumul + h
      if (intervals == TRUE) 
        CIh.cumul = CIh.cumul + CIh
      if (obs == dim(newdata)[1]) {
        h.avg <- h.cumul/dim(newdata)[1]
        if (intervals == TRUE) {
          CIh.avg <- CIh.cumul/dim(newdata)[1]
          CIh.avg <- rowQuantiles(CIh.avg, probs = c(prob.lev/2, 
                                                     1 - prob.lev/2), na.rm = TRUE)
        }
        if (plot.out) {
          if (intervals == TRUE) 
            CIh.avg.temp = matrix(CIh.avg, nrow = length(tv))
          if (is.null(ylim) && intervals == TRUE) 
            ylim <- c(min(CIh.avg.temp[, 1]), max(CIh.avg.temp[, 
                                                               2]))
          if (is.null(ylim) && intervals == FALSE) 
            ylim <- c(min(h.avg), max(h.avg))
          if (missing(ylab)) 
            ylab <- "Hazard"
          if (missing(xlab)) 
            xlab <- "Time"
          if (missing(pch)) 
            pch = 19
          if (length(tv) == 1 | bars == TRUE) 
            type.plot = "p"
          else type.plot = "l"
          plot(tv, h.avg, type = type.plot, ylab = ylab, 
               xlab = xlab, ylim = ylim, pch = pch, ...)
          if (intervals == TRUE) {
            if (length(tv) == 1 | bars == TRUE) {
              for (i in 1:length(tv)) lines(c(tv[i], 
                                              tv[i]), CIh.avg.temp[i, ])
            }
            else {
              if (shade == FALSE) {
                lines(tv, CIh.avg.temp[, 1], lty = 2)
                lines(tv, CIh.avg.temp[, 2], lty = 2)
              }
              else {
                polygon(c(tv, rev(tv)), c(CIh.avg.temp[, 
                                                       1], rev(CIh.avg.temp[, 2])), col = "gray80", 
                        border = NA)
                lines(tv, h.avg, type = "l")
              }
            }
          }
        }
      }
    }
    if (type == "cumhaz") {
      ch <- -log(pr)
      if (intervals == TRUE) {
        prs <- -log(prs)
        for (i in 1:length(tv)) {
          poe <- which(prs[i, ] %in% boxplot.stats(prs[i, 
                                                       ])$out)
          prs[, poe] <- NA
          poe <- union(poe, poet)
          poet <- poe
        }
        CIch <- prs
      }
      ch.cumul = ch.cumul + ch
      if (intervals == TRUE) 
        CIch.cumul = CIch.cumul + CIch
      if (obs == dim(newdata)[1]) {
        ch.avg = ch.cumul/dim(newdata)[1]
        if (intervals == TRUE) {
          CIch.avg = CIch.cumul/dim(newdata)[1]
          CIch.avg <- rowQuantiles(CIch.avg, probs = c(prob.lev/2, 
                                                       1 - prob.lev/2), na.rm = TRUE)
        }
        if (plot.out) {
          if (intervals == TRUE) 
            CIch.avg.temp = matrix(CIch.avg, nrow = length(tv))
          if (is.null(ylim) && intervals == TRUE) 
            ylim <- c(min(CIch.avg.temp[, 1]), max(CIch.avg.temp[, 
                                                                 2]))
          if (is.null(ylim) && intervals == FALSE) 
            ylim <- c(min(ch.avg), max(ch.avg))
          if (missing(ylab)) 
            ylab <- "Cumulative Hazard"
          if (missing(xlab)) 
            xlab <- "Time"
          if(missing(pch)) pch <- 19            
          if (length(tv) == 1 | bars == TRUE) 
            type.plot = "p"
          else type.plot = "l"
          plot(tv, ch.avg, type = type.plot, ylab = ylab, 
               xlab = xlab, ylim = ylim, pch = pch, ...)
          if (intervals == TRUE) {
            if (length(tv) == 1 | bars == TRUE) {
              for (i in 1:length(tv)) lines(c(tv[i], 
                                              tv[i]), CIch.avg.temp[i, ])
            }
            else {
              if (shade == FALSE) {
                lines(tv, CIch.avg.temp[, 1], lty = 2)
                lines(tv, CIch.avg.temp[, 2], lty = 2)
              }
              else {
                polygon(c(tv, rev(tv)), c(CIch.avg.temp[, 
                                                        1], rev(CIch.avg.temp[, 2])), col = "gray80", 
                        border = NA)
                lines(tv, ch.avg, type = "l")
              }
            }
          }
        }
      }
    }
  }
  out.r <- list(s = pr.avg, h = h.avg, ch = ch.avg, h.sim = hs, 
                s.sim = prs, l.poe = length(poe), CIs = CIpr.avg, CIh = CIh.avg, 
                CIch = CIch.avg)
  invisible(out.r)
}
