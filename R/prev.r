prev <- function(x, sw = NULL, joint = TRUE, n.sim = 100, prob.lev = 0.05){


if(joint == TRUE) type <- "joint" 
if(joint == FALSE) type <- "univariate" 


if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous/discrete margins.")
if(x$Cont == "NO" && x$VC$ccss == "yes" ) stop("This function is not suitable for selection models with continuous/discrete margin.")


# if(any(x$margins != "probit")) stop("You have to use probit margins for this calculation to make sense.")


lb <- wm <- ub <- qz <- sv <- Vv <- G <- X2sg <- Xsg <- 1
wms <- NA


if((x$Model=="B" || x$Model=="BPO") && x$triv == FALSE) stop("This function is suitable for sample selection models only.")

if(x$triv == TRUE && x$Model != "TSS") stop("This function is suitable for sample selection models only.")


# if( !( type %in% c("naive","univariate","joint") ) ) stop("Error in parameter type value. It should be one of: naive, univariate or joint.")


if(x$Model == "BSS") Xsg <- x$X2s
if(x$Model == "TSS") Xsg <- x$X3s


if(type == "univariate"){

if(x$Model == "BSS"){
        etasg <- Xsg%*%x$gam2$coefficients
        Vv    <- x$gam2$Vp
                    }

if(x$Model == "TSS"){
        etasg <- Xsg%*%x$gam3$coefficients
        Vv    <- x$gam3$Vp
                    }

}  



if(type == "joint" || type == "naive"){ # naive useful for sw below

if(x$Model == "BSS") etasg <- x$eta2
if(x$Model == "TSS") etasg <- x$eta3

Vv <- x$Vb

}




if(!is.null(sw)) { if( length(sw)!=length(etasg) ) stop("sw must have the same length as the number of observations used in fitting.")  }
if(is.null(sw)) sw <- rep(1,length(etasg)) 




#######

wm <- weighted.mean(probm(etasg, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr, w=sw)
names(wm) <- "prev"

#######



if(type != "naive"){



  if(type == "joint") coefm <- x$coefficients    
  
  
  if(type == "univariate"){ 
  
    if(x$Model == "BSS") coefm <- x$gam2$coefficients
    if(x$Model == "TSS") coefm <- x$gam3$coefficients
  
  }
  
  
   bs <- rMVN(n.sim, mean = coefm, sigma=Vv)
  
  
  if(type == "joint"){
  
    if(x$Model == "BSS") bs <- bs[, x$X1.d2 + (1 : x$X2.d2) ]
    if(x$Model == "TSS") bs <- bs[, x$X1.d2 + x$X2.d2 + (1 : x$X3.d2) ]

  }
  
  
  
  
  
 
  ps  <- probm( Xsg%*%t(bs) , x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 
  wms <- colWeightedMeans( ps, w = sw, na.rm = FALSE)
  bb <- quantile(wms, probs = c(prob.lev/2,1-prob.lev/2), na.rm=TRUE )

  lb <- bb[1]
  ub <- bb[2] 
  



} # end type







if( type == "naive"){ 

if(x$Model == "BSS") inde <- x$inde
if(x$Model == "TSS") inde <- x$inde2
sw <- sw[inde]
if(x$Model == "BSS") resp <- x$y2
if(x$Model == "TSS") resp <- x$y3




qz <- qnorm(prob.lev/2, lower.tail = FALSE)
wm <- weighted.mean(resp, w = sw)
sv <- sqrt( (wm*(1 - wm))/length(resp) )
lb <- wm - qz*sv 
ub <- wm + qz*sv 
  
}
  

  res <- c(lb, wm, ub)
  

  rm(lb, ub, qz, sv, Vv, G, Xsg)

  out <- list(res = res, prob.lev = prob.lev, sim.prev = wms, prev = wm)
 
  class(out) <- "prev"

  out

}

























