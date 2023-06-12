LM.bpm <- function(formula, data = list(), weights = NULL, subset = NULL, model, hess = TRUE){

  
  Model <- model
  sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- NULL  
  
  end <- data <- var <- G <- var.eig <- params <- resf <- VC <- respvec <- qu.mag <- X1 <- X2 <- gam1 <- gam2 <- 0
  
  BivD <- "N"
  fp <- FALSE
  
  if(!(Model %in% c("B", "BSS")) || missing(Model)) stop("Error in parameter Model value. It should be one of: B or BSS.")
  if(length(formula) > 2) stop("This test is not designed for varying correlation coefficient models.")


  ig <- interpret.gam(formula)
  mf <- match.call(expand.dots = FALSE)
  
  pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names,ig[[2]]$response))
  fake.formula <- paste(ig[[1]]$response, "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(ig$fake.formula)
  mf$formula <- fake.formula  
  mf$Model <- mf$hess <- NULL  
  mf$drop.unused.levels <- TRUE 
  if(Model=="BSS") mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(Model=="BSS"){ 
     indS <- as.logical(data[,ig[[1]]$response])==FALSE 
     indS <- ifelse( is.na(indS), FALSE, indS) 
     data[indS, ig[[2]]$response] <- ifelse( is.na(data[indS, ig[[2]]$response]), 0, data[indS, ig[[2]]$response]) 
     data <- na.omit(data)
     }
  
  if(is.null(weights)) weights <- rep(1,dim(data)[1]) else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  
  if(Model=="B"){
  if(ig[[1]]$response %in% ig[[2]]$pred.names ) end <- 1
  if(ig[[2]]$response %in% ig[[1]]$pred.names ) end <- 2
  }
  

  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), weights=weights, 
                              data=data),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp else sp1 <- NULL 
  
  
  inde <- rep(TRUE, n)


  if(Model=="B"){
  
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link="probit"), weights=weights, 
                           data=data),list(weights=weights))) 
  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 

  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)
  cy1.y2 <- (1-y1)*y2
  cy1.cy2 <- (1-y1)*(1-y2)

  func.opt <- bprobgHs                       
  
  }
  

  
  if(Model=="BSS"){

  inde <- y1 > 0
  gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), weights=weights, 
                              data=data, subset=inde),list(weights=weights,inde=inde)))                              
  X2.d2 <- length(gam2$coefficients)
  X2 <- model.matrix(gam2) 
  y2 <- gam2$y # rep(0,length(inde)); y2[inde] <- gam2$y
  l.sp2 <- length(gam2$sp)
  
  cy1 <- (1-y1)
  y1.y2 <- y1[inde]*y2
  y1.cy2 <- y1[inde]*(1-y2)
  
  func.opt <- bprobgHsSS 

  }

  if(l.sp2 != 0) sp2 <- gam2$sp else sp2 <- NULL 
  gp1 <- gam1$nsdf
  gp2 <- gam2$nsdf   
  
  
if( l.sp1!=0 || l.sp2!=0){ 
  
 sp <- c(sp1, sp2)
 GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = NULL, gam4 = NULL, 
             gam5 = NULL, gam6 = NULL, gam7 = NULL, gam8 = NULL, gam9 = NULL)  
             
 L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = 0, l.sp4 = 0, 
             l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0, l.sp9 = 0)    
             
 L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = 0, l.gam3 = 0, l.gam4 = 0,
              l.gam5 = 0, l.gam6 = 0, l.gam7 = 0, l.gam8 = 0, l.gam9 = 0)             
 
 qu.mag <- S.m(GAM, L.SP, L.GAM)               
                

                           }


  respvec <- list(y1 = y1,
                  y2 = y2,
                  y1.y2 = y1.y2, 
                  y1.cy2 = y1.cy2, 
                  cy1.y2 = cy1.y2, 
                  cy1.cy2 = cy1.cy2, 
                  cy1 = cy1)
  
  VC <- list(X1 = X1, 
             X2 = X2, X3 = NULL, inde = inde,
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             gp1 = gp1, 
             gp2 = gp2, gp3 = NULL,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2, l.sp3 = 0,
             weights = weights,
             hess = hess,
             Model = Model,
             end = end, fp = fp,
             BivD = BivD, nC = 1, extra.regI = FALSE, margins = c("probit","probit"),
             bl = c("probit", "logit", "cloglog", "cauchit"), triv = FALSE, univ.gamls = FALSE , n = n,
             min.dn = 1e-323, min.pr = 1e-32, max.pr = 0.9999999)


params <- c(gam1$coefficients, gam2$coefficients,0)


l.splist <- list( l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = 0, 
                  l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, 
                  l.sp7 = 0, l.sp8 = 0 , l.sp9 = 0)


if( l.sp1==0 && l.sp2==0 ) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(qu.mag, sp, VC, univ = 0, l.splist)


resf <- func.opt(params, respvec, VC, ps)

G   <- resf$gradient
var <- resf$hessian

tolH <- sqrt(.Machine$double.eps)

var.eig <- eigen(var, symmetric=TRUE)   
if(min(var.eig$values) < tolH) var.eig$values[which(var.eig$values < tolH)] <- tolH
var <- var.eig$vec%*%tcrossprod(diag(1/var.eig$val),var.eig$vec)  

ev <- as.numeric(t(G)%*%var%*%G)

rm(data, var, G, var.eig, params, resf, VC, respvec, qu.mag, X1, X2, gam1, gam2 )

return(pchisq(ev,1,lower.tail=FALSE))

}






