hazsurv.plot <- function(x, eq, newdata, type = "surv", t.range = NULL, intervals = TRUE, n.sim = 100, prob.lev = 0.05, 
                         shade = FALSE, ylim, ylab, xlab, ls = 100, baseline = FALSE, 
                         pop.name = NULL, pop.min = NULL,  pop.max = NULL, pop.bin = NULL, pop.build = FALSE, pop.grid = 200, 
                         min.dn = 1e-200, min.pr = 1e-200, max.pr = 1, ...){

pr <- h <- hs <- prs <- CIpr <- CIh <- poe <- poet <- NULL
pr.avg <- h.avg <- ch.avg <- CIpr.avg <- CIh.avg <- CIch.avg <- NULL


toleps <- 1e-04

if(x$univar.gamlss == FALSE && x$surv.flex == TRUE && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$bl) ) eq <- 2
if(missing(eq) && x$univar.gamlss == FALSE) stop("You must provide the equation number (either 1 or 2).")
if(x$univar.gamlss == TRUE) eq <- 1
if(missing(newdata))        stop("You have to provide a new data frame.")   
if(!is.data.frame(newdata)) stop("You have to provide a new data frame.")   

# If we are in individual net survival case or in average net survival with specifically  
# built data, only 1 row of data needs to be provided
if( (is.null(pop.name) && dim(newdata)[1] != 1) || (!is.null(pop.name) && pop.build == TRUE && dim(newdata)[1] != 1) )  stop("Your data frame has to contain one row only.")


if(!(type %in% c("surv", "hazard", "cumhaz"))) stop("The type argument can either be surv, hazard or cumhaz")
if(x$surv.flex == FALSE) stop("This function is only suitable for flexible survival models.")
if(missing(ylim)) ylim <- NULL   


if(!is.null(pop.name) && is.null(pop.min) && is.null(pop.max) && is.null(pop.bin)) stop("You indicated the population variable with respect to which you want to \n stratify but didn\'t give the range/value you want to consider.")


# In the setting in which user has to provide his/her own data, raise error if the data isn't   
# coherent (i.e. it doesn't contain the variable wrt which he/she wants to stratify)
if(!is.null(pop.name) && pop.build == FALSE  && !(pop.name %in% names(newdata)) ) stop("The population variable you have chosen doesn\'t appear in the data provided.")

# In the setting in which the data has to be built, raise error if the user provided the stratification
# variable among the other variables given in (1-row!) newdata
if(!is.null(pop.name) && pop.build == TRUE  && (pop.name %in% names(newdata)) ) stop("As you are in a build setting, you need not provide \n the variable with respect to which you want to stratify in newdata.")

# Assuming it doesn't make sense to want a binary stratification variable in the pop.build = TRUE setting,
# check whether the data provided is actually binary
if(!is.null(pop.name) && !is.null(pop.bin) && pop.build) stop("pop.build setting isn\'t allowed when stratification variable is binary.")
if(!is.null(pop.name) && !is.null(pop.bin) && !all(unique(newdata[, pop.name]) %in% c(0,1)) ) stop(paste("The data you provided for variable", pop.name, "isn't binary." ))




#################################################################################################


if(eq == 1){
  ntv  <- as.character(x$formula[[1]][2])
  
  #rlb <- range(x$y1)[1]
 
  if(is.null(t.range)){  if(x$univar.gamlss == TRUE) rlb <- x$rangeSurv[1] else rlb <- range(x$y1)[1] } else rlb <- t.range[1]   
 
  rlb <- ifelse(rlb < toleps, toleps, rlb)
  
  
  
  if(is.null(t.range)){ 
  
  tv <- seq(rlb, range(x$y1)[2], length.out = ls)  
  
  if(x$univar.gamlss == TRUE) tv <- seq(rlb, x$rangeSurv[2], length.out = ls) else tv <- seq(rlb, range(x$y1)[2], length.out = ls)

  
  } else tv <- seq(rlb, t.range[2], length.out = ls)
  
  
  indp <- 1:x$VC$X1.d2
  gob  <- x$gam1
}

if(eq == 2){
  ntv  <- as.character(x$formula[[2]][2])
  
  if(is.null(t.range)){ rlb <- range(x$y2)[1] } else rlb <- t.range[1] 
  rlb <- ifelse(rlb < toleps, toleps, rlb)
  if(is.null(t.range)){ tv   <- seq(rlb, range(x$y2)[2], length.out = ls) } else tv <- seq(rlb, t.range[2], length.out = ls)  
  
  indp <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
  gob  <- x$gam2
}


ti <- data.frame(tv)
names(ti) <- ntv



#################################################################################################
#################################################################################################

# Extract data needed
if(!is.null(pop.name) && !pop.build){ # Case A: when data is provided by user, activate selection process
  
  
  if(is.null(pop.bin)){ # Case A.1: non-binary stratification variable
    if(!is.null(pop.min) && !is.null(pop.max)) ind.pop = newdata[, pop.name] >= pop.min & newdata[, pop.name] <= pop.max
    if(!is.null(pop.min) &&  is.null(pop.max)) ind.pop = newdata[, pop.name] >= pop.min
    if( is.null(pop.min) && !is.null(pop.max)) ind.pop = newdata[, pop.name] <= pop.max
  } else { # Case A.2: binary stratification variable
    if(pop.bin == 0) ind.pop = newdata[, pop.name] == 0
    if(pop.bin == 1) ind.pop = newdata[, pop.name] == 1
  }
  
  
  newdata.complete = newdata[ind.pop, ]
  row.names(newdata.complete) = NULL # this is to reset row indices
  
  if(dim(newdata)[1] == 0 ) stop("There are 0 observations in the (sub-)population you have chosen.")
  
} else if(!is.null(pop.name) && pop.build) { # Case B: when pop.build == TRUE, build the data needed
  
  strat.var = seq(from = pop.min, to = pop.max, length.out = pop.grid)
  newdata.complete = data.frame(strat.var, newdata)
  names(newdata.complete) = c(pop.name, names(newdata))
  
} else { # Case C: traditional individual case, so just save newdata as given
  
  newdata.complete = newdata 
  
}

#################################################################################################
#################################################################################################

#############################################
# Cumulative quantities which we will use at the end of the 'for cycle' for the plot
# Note1: if we are plotting individual surv/hazard plot this will just be the individual curve (so we can think of these
# cumulative quantities as more general versions of the base ones)
# Note2: we define three seperate couples of quantities for clarity. Only one out of three couples will actually be used
# (depending on type chosen)

pr.cumul <- h.cumul <- ch.cumul <- rep(0, ls)
if(intervals == TRUE) CIpr.cumul <- CIh.cumul <- CIch.cumul <- matrix(0, ls, n.sim)

#############################################


  if(!is.null(x$VC$mono.sm.pos)) mono.sm.pos <- x$VC$mono.sm.pos else mono.sm.pos <- c(x$VC$mono.sm.pos1, x$VC$mono.sm.pos2 + x$VC$X1.d2)  
  
  # we only sample once
  
  if(intervals == TRUE){
  
      bs                <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
      bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos]) 
      
  }    

#############################################


for(obs in 1:dim(newdata.complete)[1]){ 
  
  if(obs %% 30 == 0) print(paste(round(obs/dim(newdata.complete)[1]*100), '% of iterations complete', sep = ''))
  
  newdata = as.data.frame(newdata.complete[obs,])
  row.names(newdata) = NULL
  names(newdata) <- names(newdata.complete)
  
  newdata <- data.frame(ti, newdata)
  
  Xpred <- predict(x, newdata, eq = eq, type = "lpmatrix")
  
  if(baseline == TRUE){                                   # should be general enough but it may need checking in future ...
  
     Xd            <- Xdpred(gob, newdata, ntv)
     ind0          <- (colSums(Xd == 0) == dim(Xpred)[1])
     ind0[1]       <- FALSE                               # intercept must stay to have surv prob properly scaled
     Xpred[, ind0] <- 0
  
  } 
  
  
  params1 <- x$coef.t[indp]
  eta1 <- Xpred%*%params1
  pd <- probmS(eta1, x$VC$margins[eq], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)    
  pr <- pd$pr 
  
  
  
  
  if(intervals == TRUE){
    
      eta1s <- Xpred%*%t(bs[,indp])
      pds <- probmS(eta1s, x$VC$margins[eq], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr) 
      prs <- pds$pr 
  
  }
  
  
  #################################################################################################
  
  if(type == "surv"){
    
   if(intervals == TRUE){
    
    # safety check
    for(i in 1:ls){ 
      poe <- which(prs[i,] %in% boxplot.stats(prs[i,])$out) 
      prs[, poe] <- NA 
      poe  <- union(poe, poet) 
      poet <- poe
                  }
     
    CIpr <- prs 
   }
    
    ###########################################################
    # Keeping track of confidence intervals and survival curves
    
    pr.cumul = pr.cumul + pr
    if(intervals == TRUE) CIpr.cumul = CIpr.cumul + CIpr
    
    ###########################################################
    
    
    if(obs == dim(newdata.complete)[1]){ # i.e. if we have arrived to the end of the for cycle, then produce plot
      
      # Obtain average quantities from the cumulative ones
      
      pr.avg = pr.cumul/dim(newdata.complete)[1]
      
      if(intervals == TRUE){ CIpr.avg <- CIpr.cumul/dim(newdata.complete)[1]
                             CIpr.avg <- rowQuantiles(CIpr.avg, probs = c(prob.lev/2, 1-prob.lev/2), na.rm = TRUE)
                           }
      
      
      if(is.null(ylim) && intervals == TRUE)  ylim <- c(min(CIpr.avg[,1]),max(CIpr.avg[,2]) )  
      if(is.null(ylim) && intervals == FALSE) ylim <- c(min(pr.avg),      max(pr.avg) )  
      if(missing(ylab))                       ylab <- "Survival function"  
      if(missing(xlab))                       xlab <- "Time"  
      
      plot(tv, pr.avg, type = "l", ylab = ylab, xlab = xlab, ylim = ylim, ...)
      
      if(intervals == TRUE){
      
        if(shade == FALSE){
          lines(tv, CIpr.avg[,1], lty = 2)
          lines(tv, CIpr.avg[,2], lty = 2)
        }else{
          polygon(c(tv, rev(tv)), c(CIpr.avg[,1],rev(CIpr.avg[,2])), col= "gray80", border = NA)
          lines(tv, pr.avg, type = "l")
        }
        
      }
      
    }
      
      
    }
    

  
  
  
  #################################################################################################
  
  if(type == "hazard"){
    
    if(baseline == FALSE) Xd <- Xdpred(gob, newdata, ntv) # if true, already calculated above
    
    Xthe <- Xd%*%params1   
    
    Gp <- pd$dS
    h  <- -Gp/pr*Xthe
    
    
    if(intervals == TRUE){
    
      Gps   <- pds$dS
      Xthes <- Xd%*%t(bs[,indp]) 
      hs    <- -Gps/prs*Xthes
    
      # safety check
      for(i in 1:ls){ 
         poe <- which(hs[i,] %in% boxplot.stats(hs[i,])$out)
         hs[, poe] <- NA 
         poe  <- union(poe, poet)
         poet <- poe 
                    }
    
      CIh <- hs 
      CIh <- ifelse(CIh < 0, 0, CIh)
    
    }
    
    ###########################################################
    # Keeping track of confidence intervals and hazard curves
    
    h.cumul = h.cumul + h
    if(intervals == TRUE) CIh.cumul = CIh.cumul + CIh
    
    ###########################################################
    
    if(obs == dim(newdata.complete)[1]){ # i.e. if we have arrived to the end of the for cycle, then produce plot
    
      h.avg <- h.cumul/dim(newdata.complete)[1]
      
      if(intervals == TRUE){ CIh.avg <- CIh.cumul/dim(newdata.complete)[1]
                             CIh.avg <- rowQuantiles(CIh.avg, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
                            }
                            
      if(is.null(ylim) && intervals == TRUE)  ylim <- c(min(CIh.avg[,1]), max(CIh.avg[,2]) )  
      if(is.null(ylim) && intervals == FALSE) ylim <- c(min(h.avg),       max(h.avg) )  
      if(missing(ylab))                       ylab <- "Hazard"  
      if(missing(xlab))                       xlab <- "Time"  
      
      plot(tv, h.avg, type = "l", ylab = ylab, xlab = xlab, ylim = ylim, ...)
      
      if(intervals == TRUE){
      
        if(shade == FALSE){
          lines(tv, CIh.avg[,1], lty = 2)
          lines(tv, CIh.avg[,2], lty = 2)
        }else{
          polygon(c(tv,rev(tv)), c(CIh.avg[,1],rev(CIh.avg[,2])), col= "gray80", border = NA)
          lines(tv, h.avg, type = "l")
        }
      }
    
    }
    
    
    
    
  }
  
  
  #################################################################################################
  
  
  if(type == "cumhaz"){
    
    ch <- -log(pr)
    
    
    if(intervals == TRUE){ 
    
        prs <- -log(prs)
    
        for(i in 1:ls){ 
          poe <- which(prs[i,] %in% boxplot.stats(prs[i,])$out)
          prs[, poe] <- NA 
          poe  <- union(poe, poet)
          poet <- poe 
                      }  
        
    CIch <- prs
   
    }
    
    ###########################################################
    # Keeping track of confidence intervals and hazard curves
    
    ch.cumul = ch.cumul + ch
    if(intervals == TRUE) CIch.cumul = CIch.cumul + CIch
    
    ###########################################################
    
    if(obs == dim(newdata.complete)[1]){ # i.e. if we have arrived to the end of the for cycle, then produce plot
    
      ch.avg = ch.cumul/dim(newdata.complete)[1]
      if(intervals == TRUE){ CIch.avg = CIch.cumul/dim(newdata.complete)[1]
                             CIch.avg <- rowQuantiles(CIch.avg, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
                           }
      
      
      
      # if(is.null(ylim) && intervals == TRUE)  ylim <- c(min(CIch.cumul[,1]),max(CIch.cumul[,2]) )
      
      if(is.null(ylim) && intervals == TRUE)  ylim <- c(min(CIch.avg[,1]), max(CIch.avg[,2]) )
      if(is.null(ylim) && intervals == FALSE) ylim <- c(min(ch.avg),max(ch.avg) )  
      if(missing(ylab))                       ylab <- "Cumulative Hazard"  
      if(missing(xlab))                       xlab <- "Time"  
      
      plot(tv, ch.avg, type = "l", ylab = ylab, xlab = xlab, ylim = ylim, ...)
      
      if(intervals == TRUE){
        if(shade == FALSE){
          lines(tv, CIch.avg[,1], lty = 2)
          lines(tv, CIch.avg[,2], lty = 2)
        }else{
          polygon(c(tv,rev(tv)), c(CIch.avg[,1],rev(CIch.avg[,2])), col= "gray80", border = NA)
          lines(tv, ch.avg, type = "l")
        }
      }
      
      
      }
    
  
  }
  
 
  
   
}

#################################################################################################
#################################################################################################

out.r <- list(s = pr.avg, h = h.avg, ch = ch.avg, h.sim = hs, s.sim = prs, l.poe = length(poe), 
              CIs = CIpr.avg, CIh = CIh.avg, CIch = CIch.avg)
invisible(out.r)

}

