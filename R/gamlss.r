gamlss <- function(formula, data = list(), weights = NULL, subset = NULL,  
                   margin = "N", surv = FALSE, cens = NULL, type.cens = "R", upperB = NULL,  
                   robust = FALSE, rc = 3, lB = NULL, uB = NULL, infl.fac = 1, 
                   rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                   gc.l = FALSE, parscale, extra.regI = "t", gev.par = -0.25,
                   chunk.size = 10000, k.tvc = 0, knots = NULL,
                   informative = "no", inform.cov = NULL, margin2 = "PH", fp = FALSE,
                   sp = NULL, drop.unused.levels = TRUE,
                   siginit = NULL, shinit = NULL,
                   sp.method = "perf"){ 
  
  ##########################################################################################################################
  # preamble 
  ##########################################################################################################################
    
    
  if(!is.null(sp)) sp.fixed <- sp else sp.fixed <- NULL  
  #if(!is.null(sp) && robust == FALSE) stop("Smoothing parameters can not be fixed.")  

  
  v.rB <- upperB
  
  r.type <- "a"  
  i.rho <- sp <- qu.mag <- qu.mag1 <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- spgamlss1 <- indexT <- test.sv.inf <- fgam <- NULL  
  end <- X2.d2 <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- y1m <- y2m <- Xi <- X1ni <- X2ni <- gam1TW <- NULL
  gp2 <- gp3 <- 1
  sp1 <- sp2 <- gam2 <- X2 <- sp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- inde.inf2 <- inde.inf1 <- infsetupR <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL
  Xd2 <- mono.sm.pos2 <- Sl.sf <- rangeSurv <- NULL
  Xd1 <- Xd <- mono.sm.pos <- gam1TW <- NULL
  indvU <- indvR <- indvL <- indvI <- NULL
  Sl.sf1 <- Sl.sf2 <- Sl.sf3 <- list()
  tfc <- no.pb <- NA
  surv.flex <- FALSE
  tempb <- NULL
  D <- pos.pb <- list()
  Gmat12 <- Hmat12 <- NULL
  
  
  #if(robust == TRUE) sp.method = "efs" # maybe drop option above?

  m2  <- c("N","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK","GP","GPII","GPo")
  m3  <- c("DAGUM","SM","TW")
  m1d <- c("PO", "ZTP", "GEVlink")
  m2d <- c("NBI", "NBII","PIG","DGP","DGPII")
  m3d <- c("DEL","SICHEL")
  
  if(margin == "PH" && surv == TRUE) margin <- "cloglog"
  if(margin == "PO" && surv == TRUE) margin <- "logit"  
  
  bl  <- c("probit", "logit", "cloglog")
  
  
  if(surv == TRUE && informative == "yes"){
  
   if(margin2 == "PH") margin2 <- "cloglog"
   if(margin2 == "PO") margin2 <- "logit"   
  
  }
  
  
  

  ##########################################################################################################################
  if(!is.list(formula)) stop("You must specify a list of one or more equations.")
  l.flist <- length(formula)
  
  if(l.flist != 3 && margin == "TW") stop("You must specify three equations.")
  
  
  
  if(surv == TRUE && margin %in% bl && informative == "yes"){
  
  form.check(formula, l.flist, gamlss = FALSE) 
  
  if(all.equal(formula[[1]], formula[[2]]) != TRUE) stop("The two formulae have to be the same. Get in touch for more info.")
  
  } else form.check(formula, l.flist, gamlss = TRUE)
  
  
  cl <- match.call()  
  mf <- match.call(expand.dots = FALSE)

  if(surv == TRUE && margin %in% bl && informative == "yes") pred.varR <- pred.var(formula, l.flist, gaml = TRUE, informative = "yes") else pred.varR <- pred.var(formula, l.flist, gaml = TRUE) 
   
  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  
  
  
  ###
  if(!is.null(v.rB)) pred.n <- c(pred.n, v.rB)
  ###
  

  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  
  mf$upperB <- mf$type.cens <- mf$sp.method <- mf$siginit <- mf$shinit <- mf$ordinal <- mf$sp <- mf$fp <- mf$lB <- mf$uB <- mf$margin2 <- mf$informative <- mf$inform.cov <- mf$knots <- mf$k.tvc <- mf$chunk.size <- mf$gev.par <- mf$surv <- mf$robust <- mf$rc <- mf$margin <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- NULL                           
  mf$drop.unused.levels <- drop.unused.levels
  
  #########
  if( surv == TRUE && type.cens %in% c("I","mixed") ) mf$na.action <- na.pass # this might be a bit weak
  #########
  
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
  

          
  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]  
  
  
  if(surv == TRUE && !("(cens)" %in% names(data)) ) stop("You must provide the binary censoring indicator.")
       
       
  if(!("(cens)" %in% names(data))) {cens <- rep(1,dim(data)[1]) 
                        data$cens <- cens
                        names(data)[length(names(data))] <- "(cens)"} else cens <- data[,"(cens)"]                         
  
  M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, robust = robust, extra.regI = extra.regI, margin = margin,
            surv = surv, cens = cens, bl = bl, informative = informative, list.inf.cov = inform.cov, sp.method = sp.method,
            type.cens = type.cens, v.rB = v.rB)  
            
  pream.wm(formula, margins = NULL, M, l.flist, type = "gamls")
  
  
  
  if( surv == TRUE && type.cens %in% c("I") ){
  
     data[cens == 1, v.rB] <- data[cens == 1, v1[1] ]    
     data <- na.omit(data)  
  
  }  
 
  if( surv == TRUE && type.cens %in% c("mixed") && any(unique(cens) == "I") ){

     data[cens != "I", v.rB] <- data[cens != "I", v1[1] ]    
     data <- na.omit(data)  
               
  } 
  
  
  if( surv == TRUE && !is.null(v.rB) )  rangeSurv <- range( c( data[, v1[1]], data[, v.rB] )  ) # for hazard plot
  if( surv == TRUE &&  is.null(v.rB) )  rangeSurv <- range(    data[, v1[1]]  ) 
  
  
  
  u1resp <- data[, v1[1]] # this is for survival, but generally useful
  n <- dim(data)[1]


  formula.eq1 <- formula[[1]]
  if(surv == TRUE && informative == "yes") formula.eq2 <- formula[[2]]
  
 ##############################################################  
 ##############################################################  
   
 form.eq12R <- form.eq12(formula.eq1, data, v1, margin, m1d, m2d)   # this will protect ZTP from having zeros in the response
 
 formula.eq1  <- form.eq12R$formula.eq1
 formula.eq1r <- form.eq12R$formula.eq1r
 y1           <- form.eq12R$y1
 y1.test      <- form.eq12R$y1.test  
 y1m          <- form.eq12R$y1m
 
 
if(surv == TRUE && margin2 %in% bl && informative == "yes"){

 form.eq12R2 <- form.eq12(formula.eq2, data, v2, margin2, m1d, m2d)   
 
 formula.eq2  <- form.eq12R2$formula.eq1
 formula.eq2r <- form.eq12R2$formula.eq1r
 y2           <- form.eq12R2$y1
 y2.test      <- form.eq12R2$y1.test  
 y2m          <- form.eq12R2$y1m

}
 
 ### #esres <- resp.check(y1, margin, plots = FALSE, print.par = TRUE, loglik = FALSE, os = FALSE, i.f = TRUE)
 ## est.ob <- try( gpd.fit(y1, threshold = 0, siglink = exp, show = FALSE, siginit = mean(y1 + mean(y1))/2, shinit = 0.0025), silent = TRUE)  
 #est.ob <- try(  gpd.fit(y1, threshold = 0, siglink = exp, show = FALSE), silent = TRUE) 
 #if(class(est.ob) == "try-error") esres <- resp.check(y1, margin, plots = FALSE, print.par = TRUE, loglik = FALSE, os = FALSE, i.f = TRUE)
 #if(class(est.ob) != "try-error") esres <- c( est.ob$mle[2], est.ob$mle[1]  )
 #esres <- abs(esres)
 #                                 esres <- try( resp.check(y1, margin, plots = FALSE, print.par = TRUE, loglik = FALSE, os = FALSE, i.f = TRUE), silent = TRUE)
 #if(class(esres) == "try-error") esres <- c( 0.0025, log( mean((y1 + mean(y1))/2)  )  ) 
 
 
 
 if( margin %in% c("GP","GPII","GPo","DGP","DGPII") ){ 
 
   est.ob <- try(  gpd.fit(y1, threshold = 0, siglink = exp, show = FALSE), silent = TRUE) 
 
   if(class(est.ob) == "try-error") est.ob <- try( gpd.fit(y1, threshold = 0, siglink = exp, show = FALSE, siginit = mean(y1 + mean(y1))/2, shinit = 0.0025), silent = TRUE)  
   if(class(est.ob) == "try-error") st.v <- c( 0.0025, log( mean((y1 + mean(y1))/2)  )  ) 
   if(class(est.ob) != "try-error") st.v <- c( est.ob$mle[2], est.ob$mle[1]  ) 
 
   esres <- st.v # abs(st.v)
      
   if(esres[1] < 0){
 
     esresT <- try( resp.check(y1, margin, plots = FALSE, print.par = TRUE, loglik = FALSE, os = FALSE, i.f = TRUE), silent = TRUE)
     if(class(esresT) != "try-error") esres <- esresT  
 
    }
    
   if(esres[1] < 0) esres[1] <- 0.001
 
   if( margin %in% c("GPII","GPo") ) esres[1] <- log(esres[1] + 0.5) 
   if( margin %in% c("DGPII") )      esres[1] <- sqrt(esres[1])
   
   if( !is.null(siginit) ) esres[2] <- siginit 
   if( !is.null(shinit) )  esres[1] <- shinit  

   data$estobXiGP <- esres[1] + y1/(mean(y1)*100) 
 
 }
 
 ### remember to check sigma2 para when covariates are in it, at the moment we have rnorm and not runif as above        
         
     
 if(margin != "GEVlink" && surv == FALSE){ gam1   <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))
                  
                  if(sp.method != "perf"){ gam1ff <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels, fit = FALSE),list(weights=weights)))
                                           Sl.sf1 <- Sl.setup(gam1ff)
                                           rm(gam1ff)}
                  
                                          }
 
 
 if(margin != "GEVlink" && surv == TRUE && !(margin %in% bl)){ gam1   <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights*cens, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights, cens = cens)))
                                      
                                      if(sp.method != "perf"){ gam1ff <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights*cens, data=data, knots = knots, drop.unused.levels = drop.unused.levels, fit = FALSE),list(weights=weights, cens = cens)))
                                                               Sl.sf1 <- Sl.setup(gam1ff)
                                                               rm(gam1ff)}
                                      
                                                              }
                                                             
 if(margin == "GEVlink"){ gam1 <- eval(substitute(gam(formula.eq1, binomial(link = "cloglog"), gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))
  
  if(sp.method != "perf"){ gam1ff <- eval(substitute(gam(formula.eq1, binomial(link = "cloglog"), gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels, fit = FALSE),list(weights=weights)))
                             Sl.sf1 <- Sl.setup(gam1ff)
                             rm(gam1ff)}
             
                           }






  if(margin == "TW"){
  
  gam1TW <- eval(substitute(gam(formula, gamma = infl.fac, weights = weights, data = data, knots = knots, drop.unused.levels = drop.unused.levels, family = twlss()), list(weights = weights)))
 
  # test below, Sl.sf can be set up directly or in pieces
  # gam1TWff <- eval(substitute(gam(formula, gamma = infl.fac, weights = weights, data = data, knots = knots, drop.unused.levels = drop.unused.levels, family = twlss(), fit = FALSE), list(weights = weights)))
  # Sl.sf <- Sl.setup(gam1TWff)
  # rm(gam1TWff)
  
  }


 #############################################################################################
                                     #############################################################
 if(surv == TRUE && margin %in% bl){ ######### think about this for parametric models ############
  surv.flex <- TRUE                  ######### set.seed(1)rm(list=".Random.seed", envir=globalenv())rnorm(dim(data)[1])
                                     ### no schall-f avalable here ##
  ###################################################################
  ####  cox.ph pre-fit to create response for starting value fit ####
  
  f.eq1             <- form.eq12R$f.eq1
  data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
  
  
  if(type.cens == "R") tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), data = data, weights = cens, drop.unused.levels = drop.unused.levels),list(cens=cens)))
  
  
  if( type.cens %in% c("L","I") ){
  
    f.eq1LI         <- form.eq12R$f.eq1LI
    data$temp.respV <- data[, v1[1]] # this is u1
  
  
      if(type.cens == "L") data$temp.respV[cens == 0] <-  data$temp.respV[cens == 0]/2
      
      
      if(type.cens == "I"){
            
           r.op  <- range(data$temp.respV[cens == 1])
           r.upL <- length(data$temp.respV[cens == 0])
           
           for(i in 1:r.upL){ 
           
                  if( data$temp.respV[cens == 0][i] < r.op[1] ) {r.op <- range(c(r.op, data$temp.respV[cens == 0][i]))}    
                  
                  if( data[cens == 0, v.rB][i] > r.op[2] )      { data$temp.respV[cens == 0][i] <- data[cens == 0, v.rB][i]; r.op <- range(c(r.op, data$temp.respV[cens == 0][i]))} 
                  
                  if( data$temp.respV[cens == 0][i] > r.op[1] && data[cens == 0, v.rB][i] < r.op[2]) data$temp.respV[cens == 0][i] <- (data$temp.respV[cens == 0][i] + data[cens == 0, v.rB][i])/2       
                             }
                             
                           }
      
      
      # careful here as I need to space the whole range of outcome(s) values
      # With left censoring I space more than I should, but maybe this is just fine
      # or I could include some zeros as well? We may revisit this but looks ok for the time being
            

    tempb <- eval(substitute(gam(f.eq1LI, family = cox.ph(), data = data, weights = rep(1, length(cens)), drop.unused.levels = drop.unused.levels),list(cens=cens)))

                                 }  
                                 
                                 
                                 
                                 
      
  if( type.cens %in% c("mixed") ){
  
    f.eq1LI         <- form.eq12R$f.eq1LI
    data$temp.respV <- data[, v1[1]] # this is u1
    
      # "R" is fine as it is
    
      if( any(unique(cens) == "L") ) data$temp.respV[cens == "L"] <-  data$temp.respV[cens == "L"]/2
      
      if( any(unique(cens) == "I") ){
      
      
           r.op  <-  range(data$temp.respV[cens != "I"])
           r.upL <- length(data$temp.respV[cens == "I"])
           
           for(i in 1:r.upL){ 
           
                  if( data$temp.respV[cens == "I"][i] < r.op[1] ) {r.op <- range(c(r.op, data$temp.respV[cens == "I"][i]))}    
                  
                  if( data[cens == "I", v.rB][i] > r.op[2] )      { data$temp.respV[cens == "I"][i] <- data[cens == "I", v.rB][i]; r.op <- range(c(r.op, data$temp.respV[cens == "I"][i]))} 
                  
                  if( data$temp.respV[cens == "I"][i] > r.op[1] && data[cens == "I", v.rB][i] < r.op[2]) data$temp.respV[cens == "I"][i] <- (data$temp.respV[cens == "I"][i] + data[cens == "I", v.rB][i])/2       

                             }

                                    }      

    # this will also apply when we have only two categories with U and L, or U and R for instance
    # in this case the fit may change a bit (for R) since I do not fit below the same model as if I used type.cens = R but it is ok
    
    tempb <- eval(substitute(gam(f.eq1LI, family = cox.ph(), data = data, weights = rep(1, length(cens)), drop.unused.levels = drop.unused.levels),list(cens=cens)))

  }
  
  ########
  ### the above part does not affect original u1 and u2 but uses a new temp. variable ####
  ##############
  
  
  
  data$Sh <- as.vector(mm(predict(tempb, type = "response")))
  
  #data$Sh <- as.vector(mm(rstpm2:::Shat(coxph(Surv(u, delta == 1) ~ 1, data = data, model = TRUE))))
  #does not seem to have an impact
  ###################################################################
  
  cens1 <- ifelse(cens == 0, 1e-07, cens)
  
  
  if(type.cens != "R"){
  
  data[, v1[1] ] <- data$temp.respV # new u1 is now temp.respV which covers the whole range
  attr(data, "terms") <- NULL
  
                      }
      
      
  gam1 <- eval(substitute(scam(formula.eq1, gamma = infl.fac, weights = weights*cens1, data = data), list(weights = weights, cens1 = cens1)))
       
  #############################################################
  ### check that we are using the correct smoother of time ####
  
  lsgam1 <- length(gam1$smooth)
  if(lsgam1 == 0) stop("You must use at least a monotonic smooth function of time.")
  
  clsm <- ggr <- NA 
  for(i in 1:lsgam1){ clsm[i] <- class(gam1$smooth[[i]])[1] ### is lsgam1 is and is used for summary when using tensor etc
                      #ggr[i]  <- max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn)))
                    }
  
  
  
  if( sum(as.numeric(clsm[1] %in% c("mpi.smooth")))==0 ) stop("You must have a monotonic smooth of time and it has to be the first to be included.")
  
    
  #if(clsm[1] != "mpi.smooth") stop("The first smoother must be a monotonic smooth of time.")  
  #if( sum( as.numeric(clsm %in% c("mpi.smooth")) ) != sum( ggr ) ) stop("You must use mpi smooth function(s) of time.")   
  
  
  
  
  ###########################################################

  l.sp1 <- length(gam1$sp)
  if(l.sp1 != 0) sp1 <- gam1$sp
           
  ###########################################################    
  
  #if(dim(data)[1] < 2000) sp.c <- 0.2 else sp.c <- 1/sqrt(dim(data)[1])  
  #sp1[clsm %in% c("mpi.smooth")] <- sp.c 

  #sp.c <- 1 # this will have to be improved but not sure how at the moment
            # it may not affect results though, only maybe in peculiar situations
  sp1[1] <- 1 

  gam.call <- gam1$call
  gam.call$sp <- sp1
  gam1 <- eval(gam.call)
  
  ###########################################################

 
  j <- 1
  for(i in 1:lsgam1){ 
  
    if( max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn))) != 0 && clsm[i] == "mpi.smooth" ) mono.sm.pos <- c(mono.sm.pos, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para) ) 
   
    if( max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$vn))) != 0 && clsm[i] != "mpi.smooth" ){ 
    
    # these checks are not entirely general and the user may
    # still specify the model in the wrong way
    
                                                                if( clsm[i] != "pspline.smooth" && k.tvc !=0) stop("You have to use a ps smooth to allow for doubly penalised tvc terms.")
                                                                                                                                
                                                                if( clsm[i] == "pspline.smooth"){
                                                                
                                                                pos.pb[[j]] <- c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para)
                                                                indexT      <- c(indexT, pos.pb[[j]] ) # this is fine, do not touch it, good for starting values
                                                                D[[j]]      <- diff(diag(length(pos.pb[[j]])), differences = 1)
                                                                j <- j + 1
                                                                
                                                                # scaling factor should not matter since lambda is fixed here
                                                                #D[[j]] <- diff(diag(rep(0,length(pos.pb))), differences = 1)
                                                                
                                                                }
                                                                
                                                                
                                                                                            }
  
  
  }  
  
  
 
 
 
 
if( type.cens == "mixed"){ 
  
indvU <- indvR <- indvL <- indvI <- rep(0, n)

if( any(unique(cens) == "U") ) indvU[cens == "U"] <- 1 
if( any(unique(cens) == "R") ) indvR[cens == "R"] <- 1 
if( any(unique(cens) == "L") ) indvL[cens == "L"] <- 1  
if( any(unique(cens) == "I") ) indvI[cens == "I"] <- 1   
  
  
} 
  
  
if(type.cens == "R") X1 <- predict(gam1, type = "lpmatrix")


if(type.cens %in% c("I", "L", "mixed")){

        data[, v1[1]] <- u1resp                                 # 'cause I need to evaluate at the correct u1 values and not temp.respV 
        X1 <- predict(gam1, type = "lpmatrix", newdata = data)
        
      	     if( type.cens == "I" ){ 
      	     
      		data[cens == 0, v1[1]] <- data[cens == 0, v.rB] 
                X2 <- predict(gam1, type = "lpmatrix", newdata = data)  
                
                                   }    
 
      	     if( type.cens == "mixed" ){ 
      	     
                          if( any(unique(cens) == "I") ){ 
                          
                              data[cens == "I", v1[1]] <- data[cens == "I", v.rB]  
                              X2 <- predict(gam1, type = "lpmatrix", newdata = data)  
                            
                                                        }
                                       }     
}

                                                                     
  if(is.null(X2)) X2 <- matrix(1, dim(X1)[1], dim(X1)[2])       
  
  
  
  if( !is.null(indexT) && k.tvc !=0){ if(range(X1[, indexT])[1] < 0) stop("Check design matrix for smooth(s) of tvc terms.")}
  
  data[, v1[1]] <- u1resp # to go back to normal
  
  Xd <- Xdpred(gam1, data, v1[1]) # this part can be analytical given that we only use B-splines
                                  # but it can also be a tprs or any other smoother so we may want to leave it generic

  start.v1 <- c( gam1$coefficients )

  gam1$y <- data[, v1[1]] 
  
  
  
  

  if(!is.null(indexT)){
  
     start.v2 <- start.v1
     start.v2[mono.sm.pos] <- exp(start.v2[mono.sm.pos])
     while( range(Xd%*%start.v2)[1] < 0 ) start.v2[indexT] <- 0.999*start.v2[indexT]   
     
     # this is to make sure that the condition is not violated from the start
      
     start.v1[indexT]  <- start.v2[indexT]  
     gam1$coefficients <- gam1$coefficients.t <- start.v1
     gam1$coefficients.t[mono.sm.pos] <- exp(gam1$coefficients.t[mono.sm.pos]) 
  }

  
  
  
  
  
  
  
  
  
  
  
  
if(informative == "yes"){  
  
  
  f.eq1 <- form.eq12R2$f.eq1
  data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
  tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), data = data, weights = 1 - cens, drop.unused.levels = drop.unused.levels),list(cens=cens))) 
  data$Sh <- as.vector(mm(predict(tempb, type = "response")))
  

  ###################################################################
  
  cens1 <- ifelse( (1 - cens) == 0, 1e-07, 1 - cens)
  gam2 <- eval(substitute(scam(formula.eq2, gamma=infl.fac, weights=weights*cens1, data=data), list(weights=weights, cens1 = cens1)))
  
  #############################################################
  ### check that we are using the correct smoother of time ####
  
  lsgam2 <- length(gam2$smooth)
  if(lsgam2 == 0) stop("You must use at least a monotonic smooth function of time.")
  
  clsm <- ggr <- NA 
  for(i in 1:lsgam2){ clsm[i] <- class(gam2$smooth[[i]])[1] ### is lsgam1 is and is used for summary when using tensor etc
                      #ggr[i]  <- max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$vn)))
                    }
  
  if( sum(as.numeric(clsm[1] %in% c("mpi.smooth")))==0 ) stop("You must have a monotonic smooth of time and it has to be the first to be included.")
  if(v1[1] %in% inform.cov) stop("Time can not be an informative covariate.") # loose check as some cases may escape...grepl etc

  
  #if(clsm[1] != "mpi.smooth") stop("the first smoother must be a monotonic smooth of time.")   
  
  ###########################################################

  l.sp2 <- length(gam2$sp)
  if(l.sp2 != 0) sp2 <- gam2$sp
           
  ###########################################################    
  
  #sp.c <- 1 
  sp2[1] <- 1

  gam.call <- gam2$call
  gam.call$sp <- sp2
  gam2 <- eval(gam.call)
  
  ###########################################################

 
  j <- 1
  for(i in 1:lsgam2){ 
  
    if( max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$vn))) != 0 && clsm[i] == "mpi.smooth" ) mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para) ) 
   
    if( max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$vn))) != 0 && clsm[i] != "mpi.smooth" ){ # this part is not important and was for experimentation
    
    
                                                                if( clsm[i] != "pspline.smooth" && k.tvc !=0) stop("You have to use a ps smooth to allow for doubly penalised tvc terms.")
                                                                                                                                
                                                                if( clsm[i] == "pspline.smooth"){
                                                                
                                                                pos.pb[[j]] <- c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para)
                                                                indexT      <- c(indexT, pos.pb[[j]] ) # this is fine, do not touch it, good for starting values
                                                                D[[j]]      <- diff(diag(length(pos.pb[[j]])), differences = 1)
                                                                j <- j + 1
                                                                
                                                                # scaling factor should not matter since lambda is fixed here
                                                                #D[[j]] <- diff(diag(rep(0,length(pos.pb))), differences = 1)
                                                                
                                                                }
                                                                
                                                                
                                                                                            }
  
  
  }  
  
     
  
  X2 <- predict(gam2, type = "lpmatrix")
  
  if( !is.null(indexT) && k.tvc !=0){ if(range(X2[, indexT])[1] < 0) stop("Check design matrix for smooth(s) of tvc terms.")}
  
  Xd2 <- Xdpred(gam2, data, v2[1]) # this part can be analytical given that we only use B-splines
                                   # but it cn also be a tprs or any other smoother so we may want to leave it generic

  start.v2 <- c( gam2$coefficients )

  gam2$y <- data[, v2[1]]  
  
 
}  
  
  
  

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
 



if(margin == "TW") start.v1 <- coef(gam1TW) 
 
 
 
 
 if(surv == TRUE && margin %in% bl && informative == "yes"){
 
  gam2$formula <- formula.eq2r  
  lsgam2 <- length(gam2$smooth)
  
  y2 <- y2.test 
  
  gp2 <- gam2$nsdf 
  X2.d2 <- dim(X2)[2]
  
 
 }
 
 
 
 
 
 
 
 
 
 
 
 
 ##############################################################  
 ##############################################################  


if(margin != "TW"){ 

  log.nu.1 <- log.sig2.1 <- NULL 

if( !(margin %in% c(m1d, bl)) ){

start.snR <- startsn(margin, y1)
    
log.sig2.1 <- start.snR$log.sig2.1; names(log.sig2.1) <- "sigma.star"

if( margin %in% c("GP","GPII","GPo","DGP","DGPII") ) log.sig2.1 <- esres[2] # this is redundant as it is done already in startsn if we employ resp.check

if( margin %in% c(m3) ){ log.nu.1 <- start.snR$log.nu.1; names(log.nu.1) <- "nu.star"}     

                               }
}


##########################################################
##########################################################

if(surv == TRUE && margin2 %in% bl && informative == "yes"){

infsetupR <- inform.setup(gam1, gam2, inform.cov, start.v1, start.v2, lsgam1, lsgam2)

start.v1 <- infsetupR$start.v1
start.v2 <- infsetupR$start.v2

Gmat12 <- matrix(NA, n, length(start.v1))
Hmat12 <- matrix(NA, length( c( start.v1, start.v2 ) ), length( c( start.v1, start.v2 ) ) )


if( !is.null(infsetupR$par.pos1) &&  is.null(infsetupR$smo.pos1) ) { Xi <- Xi <- X1[, c(infsetupR$par.pos1)]
                                                                     X1ni <- X1[, -c(infsetupR$par.pos1)] 
                                                                     X2ni <- X2[, -c(infsetupR$par.pos2)] 
                                                                     inde.inf1 <- c(infsetupR$par.pos1)
                                                                     inde.inf2 <- c(infsetupR$par.pos2)
                                                                   }

if( !is.null(infsetupR$par.pos1) && !is.null(infsetupR$smo.pos1) ){ Xi <- Xi <- X1[, c(infsetupR$par.pos1, infsetupR$smo.pos1)]
                                                                    X1ni <- X1[, -c(infsetupR$par.pos1, infsetupR$smo.pos1)] 
                                                                    X2ni <- X2[, -c(infsetupR$par.pos2, infsetupR$smo.pos2)] 
                                                                    inde.inf1 <- c(infsetupR$par.pos1, infsetupR$smo.pos1)
                                                                    inde.inf2 <- c(infsetupR$par.pos2, infsetupR$smo.pos2)                                                                    
                                                                    
                                                                  }

if(  is.null(infsetupR$par.pos1) && !is.null(infsetupR$smo.pos1) ){ Xi <- Xi <- X1[, c(infsetupR$smo.pos1)]
                                                                    X1ni <- X1[, -c(infsetupR$smo.pos1)] 
                                                                    X2ni <- X2[, -c(infsetupR$smo.pos2)] 
                                                                    inde.inf1 <- c(infsetupR$smo.pos1)
                                                                    inde.inf2 <- c(infsetupR$smo.pos2)                                                                      
                                                                    
                                                                  }


}

# for the spline part we need to be careful with the matrix design set up.
# Specifically the above is assuming that the two equations
# are the same

##########################################################
##########################################################


if(margin != "TW"){ 

if(surv == TRUE && margin %in% bl && informative == "yes") start.v1 <- c( start.v1, start.v2 ) else{ 

if(margin %in% c(m1d) )    start.v1 <- c( gam1$coefficients )
if(margin %in% c(m2,m2d) ) start.v1 <- c( gam1$coefficients, log.sig2.1           ) 
if(margin %in% c(m3,m3d) ) start.v1 <- c( gam1$coefficients, log.sig2.1, log.nu.1 ) 
                                                                                                    }

}


##############################################################  
##############################################################  
  

if(margin == "TW") log.nu.1 <- log.sig2.1 <- 0.1
  
if( l.flist > 1 && !(surv == TRUE && margin %in% bl) ){ # not used for flexible survival
    
    vo <- list(log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, n = n, drop.unused.levels = drop.unused.levels)
    overall.svGR <- overall.svG(formula, data, ngc = 2, margin, M, vo, gam1, gam2, type = "gaml", knots = knots) # in M I will have to introduce robust = TRUE?
    
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
    
    
    if(M$sp.method != "perf"){
    
    Sl.sf2 <- overall.svGR$Sl.sf2 
    Sl.sf3 <- overall.svGR$Sl.sf3 
    
    }
    
    
    
    
        
    if(margin != "TW"){ start.v1 <- overall.svGR$start.v 
                        sp2 <- overall.svGR$sp2
                        sp3 <- overall.svGR$sp3
                      }
                      
    if(margin == "TW"){ 
                        if(l.sp1 != 0){ sp1 <- gam1TW$sp[1:l.sp1];                                     names(sp1) <- names(gam1$sp)}
                        if(l.sp2 != 0){ sp2 <- gam1TW$sp[(l.sp1 + 1):(l.sp1 + l.sp2)];                 names(sp2) <- names(gam2$sp)}
                        if(l.sp3 != 0){ sp3 <- gam1TW$sp[(l.sp1 + l.sp2 + 1):(l.sp1 + l.sp2 + l.sp3)]; names(sp3) <- names(gam3$sp)}
                        
                        names(start.v1) <- names(overall.svGR$start.v)
                        
                      }
    

}





##########################################################
# SPs and penalties
##########################################################
##########################################################

if(surv == TRUE && margin %in% bl && informative == "yes" && l.sp2 > 1 && !is.null(infsetupR$scv) && !is.null(infsetupR$inds2)) sp2 <- sp2[-c(infsetupR$inds2)] 

if(surv == TRUE && margin %in% bl && informative == "yes"){ # there is always a smooth in this model's case!

pic <- NA

av <- all.vars(formula[[1]])

for(i in 1:length(inform.cov)) pic[i] <- which(inform.cov[i] == av ) 

tfor <- formula(drop.terms(terms(formula[[1]]), pic, keep.response = TRUE))

fgam <- gam(tfor, data = data, fit = FALSE, drop.unused.levels = drop.unused.levels)

pfgam <- NA

for(i in 1:length(fgam$smooth)) pfgam[i] <- fgam$smooth[[i]]$first.para  

}



spgamlss1 <- c(sp1, sp2, sp3)
GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, K1 = NULL) 
            
            
if(l.sp1 !=0 || l.sp2 !=0 || l.sp3 !=0) { ## this will always work with informative case as there is always a smooth


	L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), 
	              l.gam3 = length(gam3$coefficients), l.gam4 = 0, l.gam5 = 0, 
	              l.gam6 = 0, l.gam7 = 0, l.gam8 = 0)                                 
  
	L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = 0, 
	             l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0)                               

                
	qu.mag1 <- S.m(GAM, L.SP, L.GAM)  


        if(surv == TRUE && margin %in% bl && informative == "yes" && l.sp2 > 1 && !is.null(infsetupR$scv)){ 
        
             qu.mag1$Ss   <-   qu.mag1$Ss[-c(l.sp1 + infsetupR$inds2)]
             qu.mag1$rank <- qu.mag1$rank[-c(l.sp1 + infsetupR$inds2)] 
             qu.mag1$off  <-  qu.mag1$off[-c(l.sp1 + infsetupR$inds2)] 
                                                                                                          }   
         
        if(surv == TRUE && margin %in% bl && informative == "yes") {qu.mag1$off[(l.sp1 + 1):length(qu.mag1$off)] <- pfgam + X1.d2
                                                                    
                                                                    test.sv.inf <- start.v1[qu.mag1$off]}
         
         
        
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

pdf.test <- distrHsATDiscr(ygrid, eta.m, sigma2.m, 1, margin, y1m, robust = TRUE)$pdf2 > sqrt(.Machine$double.eps)

maxv <- max( which(pdf.test == TRUE) ) + 10 # 10 is a sort of safety thing but it could be increased

ygrid <- ygrid[1:maxv]

} else ygrid <- NULL


my.env      <- new.env()
my.env$k    <- k.tvc 
my.env$indN <- NULL
my.env$V    <- NULL



# if(margin != "TW" || robust == TRUE) Sl.sf <- list(Sl.sf1, Sl.sf2, Sl.sf3)






if(sp.method == "efs"){########


Sl.sf <- list()

LSl.sf1 <- length(Sl.sf1)
LSl.sf2 <- length(Sl.sf2)
LSl.sf3 <- length(Sl.sf3)

if( LSl.sf1 == 0 ) pa1 <- matrix(0, gp1, gp1)
if( LSl.sf2 == 0 ) pa2 <- matrix(0, gp2, gp2)

j <- 1

if( LSl.sf1 != 0 ){
        for(i in 1:LSl.sf1 ){ Sl.sf[[j]] <- Sl.sf1[[i]]
                              j <- j + 1 } 
        attr(Sl.sf, "lambda") <- attr(Sl.sf1, "lambda")      
        attr(Sl.sf, "E")      <- attr(Sl.sf1, "E")                       
        
        
}

if( LSl.sf2 != 0 ){
        for(i in 1:LSl.sf2 ){ Sl.sf[[j]] <- Sl.sf2[[i]]
                              Sl.sf[[j]]$start <- Sl.sf[[j]]$start + length(gam1$coefficients)
                              Sl.sf[[j]]$stop  <- Sl.sf[[j]]$stop  + length(gam1$coefficients)   
                              j <- j + 1 } 
        attr(Sl.sf, "lambda") <- c(attr(Sl.sf, "lambda"), attr(Sl.sf2, "lambda")) 
  
        if(LSl.sf1 != 0) attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), attr(Sl.sf2, "E"))                      
        if(LSl.sf1 == 0) attr(Sl.sf, "E") <- adiag(pa1,              attr(Sl.sf2, "E"))                      
        
        
}

if( LSl.sf3 != 0 ){
        for(i in 1:LSl.sf3 ){ Sl.sf[[j]] <- Sl.sf3[[i]]
                              Sl.sf[[j]]$start <- Sl.sf[[j]]$start + length(gam1$coefficients) + length(gam2$coefficients)
                              Sl.sf[[j]]$stop  <- Sl.sf[[j]]$stop  + length(gam1$coefficients) + length(gam2$coefficients)
                              j <- j + 1 } 
        attr(Sl.sf, "lambda") <- c(attr(Sl.sf, "lambda"), attr(Sl.sf3, "lambda"))
        
        
        if(LSl.sf1 != 0 && LSl.sf2 != 0)  attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), attr(Sl.sf3, "E"))
        if(LSl.sf1 == 0 && LSl.sf2 != 0)  attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), attr(Sl.sf3, "E"))

        if(LSl.sf1 == 0 && LSl.sf2 == 0)  attr(Sl.sf, "E") <- adiag(pa1, pa2, attr(Sl.sf3, "E"))
        if(LSl.sf1 != 0 && LSl.sf2 == 0)  attr(Sl.sf, "E") <- adiag(attr(Sl.sf, "E"), pa2, attr(Sl.sf3, "E"))
                    
}
     
     
     
}##############     


                           

  VC <- list(lsgam1 = lsgam1, ygrid = ygrid, # why lsgam1? maybe useful outside fitting functions, do not remember, check again
             lsgam2 = lsgam2, indexT = indexT, D = D, my.env = my.env, k = k.tvc, pos.pb = pos.pb,
             lsgam3 = lsgam3, r.type = r.type, Sl.sf = Sl.sf, sp.method = sp.method,
             lsgam4 = lsgam4, gam1TW = gam1TW,
             lsgam5 = lsgam5,
             lsgam6 = lsgam6,
             lsgam7 = lsgam7,
             lsgam8 = lsgam8, fgam = fgam, ad.ind = FALSE,
             X1 = X1, 
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,  
             X7 = X7, #tol.rc = tol.rc,
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
             Xd1 = Xd, Xd2 = Xd2,
             mono.sm.pos = mono.sm.pos, mono.sm.pos2 = mono.sm.pos2, 
             surv.flex = surv.flex,
             informative = informative, inform.cov = inform.cov,
             infsetupR = infsetupR, gp2.inf = NULL,
             Xi = Xi, X1ni = X1ni, X2ni = X2ni,
             inde.inf2 = inde.inf2, inde.inf1 = inde.inf1,
             Gmat12 = Gmat12, Hmat12 = Hmat12, v1pred = v1[1],
             sp.fixed = sp.fixed,
             indvU = indvU,
             indvR = indvR,
             indvL = indvL,
             indvI = indvI)





if(surv == TRUE && margin2 %in% bl && informative == "yes"){

if( !is.null(infsetupR$pcv) ) VC$gp2.inf <- gp2 - length(infsetupR$par.pos1) else VC$gp2.inf <- gp2


if(l.sp2 > 1 && !is.null(infsetupR$scv)) VC$l.sp2 <- l.sp2 - length(infsetupR$inds2)


VC$margins <- c(margin, margin2)


}




  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################
  if(margin != "GEVlink"){
  
  	if(margin %in% c(m1d, m2d, m2) )                                     func.opt1 <- bprobgHsContUniv 
  	if(margin %in% c(m3) )                                               func.opt1 <- bprobgHsContUniv3 
  	if(margin %in% c(bl) && informative == "no" && type.cens == "R")     func.opt1 <- bcontSurvGuniv
  	if(margin %in% c(bl) && informative == "no" && type.cens == "L")     func.opt1 <- bcontSurvGunivL
  	if(margin %in% c(bl) && informative == "no" && type.cens == "I")     func.opt1 <- bcontSurvGunivI
  	if(margin %in% c(bl) && informative == "no" && type.cens == "mixed") func.opt1 <- bcontSurvGunivMIXED
  	if(margin %in% c(bl) && informative == "yes")                        func.opt1 <- bcontSurvGunivInform
  
  }
  
  if(margin == "GEVlink") func.opt1 <- bprobgHsContUnivBIN

  y1.m <- y1; if(margin == "LN") y1.m <- exp(y1) 
  VC$y1 <- y1.m 

    SemiParFit <- SemiParBIV.fit(func.opt = func.opt1, start.v = start.v1, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec2, VC = VC, sp = spgamlss1, qu.mag = qu.mag1) 
    SemiParFit$robust <- robust  
      
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- gamlss.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)  
                                                          

  SemiParFit <- SemiParFit.p$SemiParFit  

  ##########################################################################################################################

if(gc.l == TRUE) gc()

  ##########################################################################################################################


cov.c(SemiParFit)


  ##########################################################################################################################

gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  # for all.terms when plotting
  ##########################################################################################################################


rm(gam1TW) 
 

L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, formula = formula, robust = robust,       
          edf11 = SemiParFit.p$edf11,   ## this is for RE
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, 
          gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, 
          iterlimsp = iterlimsp,
          weights = weights, cens = cens, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
          l.sp7 = l.sp7, l.sp8 = l.sp8,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,  
          sigma2 = SemiParFit.p$sigma2, sigma = SemiParFit.p$sigma2,  
          sigma2.a = SemiParFit.p$sigma2.a, sigma.a = SemiParFit.p$sigma2.a,  
          nu = SemiParFit.p$nu, 
          nu.a = SemiParFit.p$nu.a,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
          X6 = X6, X7 = X7, X8 = X8,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
          X7.d2 = X7.d2, X8.d2 = X8.d2,             
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Vb1 = SemiParFit.p$Vb1, Ve = SemiParFit.p$Ve, 
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
          VC = VC, magpp = SemiParFit$magpp, type.cens = type.cens,
          Cont = "YES",
          l.flist = l.flist, triv = FALSE, univar.gamlss = TRUE, call = cl, gev.par = gev.par,
          ygrid = ygrid,
          r.weights = SemiParFit$fit$d.psi, surv = surv, surv.flex = surv.flex, test.sv.inf = test.sv.inf,
          rangeSurv = rangeSurv)

class(L) <- c("gamlss","SemiParBIV","gjrm")


L


}

