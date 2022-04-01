SemiParBIV.fit  <- function(func.opt, start.v, 
                                   rinit, rmax, iterlim, iterlimsp, tolsp,
                                   respvec, VC,
                                   sp = NULL, qu.mag = NULL){ 


l.sp1 <- VC$l.sp1 
l.sp2 <- VC$l.sp2 
l.sp3 <- VC$l.sp3 
l.sp4 <- VC$l.sp4 
l.sp5 <- VC$l.sp5 
l.sp6 <- VC$l.sp6 
l.sp7 <- VC$l.sp7 
l.sp8 <- VC$l.sp8 
l.sp9 <- VC$l.sp9 


score.hist <- rp <- D <- L <- Sl.sfTemp <- St <- NULL

#if(VC$sp.method == "efs") iterlim <- 1

l.splist <- list( l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                  l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                  l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9 )


if(!is.null(VC$sp.fixed)) sp <- VC$sp.fixed  


if( ( l.sp1==0 && l.sp2==0 && l.sp3==0 && l.sp4==0 && l.sp5==0 && l.sp6==0 && l.sp7==0 && l.sp8==0 && l.sp9==0) || VC$fp==TRUE) ps <- ps1 <- list(S.h = 0, S.h1 = 0, S.h2 = 0, qu.mag = NULL) else ps <- ps1 <- pen(qu.mag, sp, VC, univ = respvec$univ, l.splist)



if(VC$triv == TRUE){

if( VC$penCor == "ridge" ) qu.mag <- ps$qu.mag # need to do this as in the ridge case I expand S in pen
                                               # not the same for lasso and alasso as these penalties
                                               # change based on parameter estimates
                                               
# can't do lasso as it changes at each iteration  
# and we need th stuff below

#if(VC$robust == TRUE) fit$bcorPost <- bcorrecPost(VC, c(fit$argument), fit)   


if( VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp


}



###########################################################

  parsc <- rep(VC$parscale, length(start.v) ); sc <- TRUE

  fit  <- fit1 <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                     respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                     iterlim = iterlim), silent = sc)   

  if(inherits(fit, "try-error") || is.null(fit$l)){
  
    fit  <- fit1 <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                       respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                       iterlim = iterlim/4), silent = sc)  
                     
        if(inherits(fit, "try-error")|| is.null(fit$l)){
        
            fit  <- fit1 <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                               respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                               iterlim = iterlim/10), silent = sc)   


if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
 
             
                                     }
                     
  }


  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL
  conv.sp <- TRUE
  #####################################################################

  #VC$fp <- TRUE
  
  
  
  
  
  
  
    if((VC$fp==FALSE && is.null(VC$sp.fixed) && (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0 || l.sp9!=0)) ){






         if(VC$sp.method == "perf"){ #########################


       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  


	  while( stoprule.SP > tolsp ){ 

             fito   <- fit$l
             o.ests <- c(fit$argument) 
             spo  <- sp             
            
             wor.c <- working.comp(fit) 
             
             if( VC$triv == TRUE && VC$penCor %in% c("lasso", "alasso") ) qu.mag <- fit$qu.mag   
               
                	bs.mgfit <- try(magic(y = wor.c$Z,  
                	                      X = wor.c$X,
                	                      sp= sp, S = qu.mag$Ss,                
                        	              off = qu.mag$off, rank = qu.mag$rank, 
                                	      gcv = FALSE,
                                	      gamma = VC$infl.fac), silent = sc)
                		if(inherits(bs.mgfit, "try-error"))      {conv.sp <- FALSE; break} 
                		if(any( is.na(bs.mgfit$sp) ) == TRUE) {conv.sp <- FALSE; break} 
                		
                		
                	sp <- bs.mgfit$sp; if(!is.null(VC$sp.fixed)) sp <- VC$sp.fixed
                	
                	iter.sp <- iter.sp + 1; names(sp) <- names(spo)
                	
                	
                	
                	
             if( VC$triv == TRUE && VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp	

             ps <- pen(qu.mag, sp, VC, univ = respvec$univ, l.splist) # in triv case, I reset ridge penalty but ok
               
             
             #if(any( is.na(o.ests) ) == TRUE) {conv.sp <- FALSE; break} 

             # VC$my.env$k <- VC$k, does not like it  
               
             fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc, respvec = respvec, VC = VC, ps = ps, blather = TRUE, iterlim = iterlim), silent = sc)                           
                          
                          if(inherits(fit, "try-error") || is.null(fit$l)){conv.sp <- FALSE

                                                         ps <- ps1 # pen(qu.mag, spo, VC, univ = respvec$univ, l.splist)

                                                         fit <- try( trust(func.opt, c(fit1$argument), rinit=rinit, rmax = rmax,  parscale = parsc,  
			                                              respvec = respvec, VC = VC, 
			                                              ps = ps, 
                                                                      blather = TRUE, iterlim = iterlim), silent = sc)  
                                                                      

if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
 
    
                                                        # break
                          } 
                          
                                                             
             iter.inner <- iter.inner + fit$iterations   	                                    
              
             if(iter.sp >= iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- abs(fit$l - fito)/(0.1 + abs(fit$l))  # max(abs(o.ests - c(fit$argument)))
                  
           
          } # end smoothing fitting loop
          

       if(VC$gc.l == TRUE) gc()   

       
       #if(!is.null(VC$sp.fixed)) bs.mgfit$sp <- bs.mgfit$sp.full <- VC$sp.fixed  

       magpp <- magic.post.proc(wor.c$X, bs.mgfit) 
       


} #########################



        #Hp <- PDef( fit$hessian - fit$S.h )$res # safety check
        #Hp <- PDef( fit$hessian )$res 
        #Hp <- Sl.initial.repara(Sl.sfTemp, Hp, inverse = TRUE) + crossprod(rp$E)
        #Hp <- Sl.initial.repara(Sl.sfTemp, Hp, inverse = TRUE) 
        #LD <- LDfun(Hp, eigen.fix)        
        #if (sum(fit$bdrop)) {
        #    q <- length(fit$bdrop)
        #    ibd <- !fit$bdrop
        #    Vtemp <- Vb
        #    Vb <- matrix(0, q, q)
        #    Vb[ibd, ibd] <- Vtemp
        #}




if(VC$sp.method == "efs"){ #########################



LDfun <- function(Hp, eigen.fix){

       rank <- dim(Hp)[1]
       D <- diag(Hp)
        if(sum(!is.finite(D)) > 0) stop("non finite values in Hessian")
        
        if (min(D) < 0) {
            Dthresh <- max(D) * sqrt(.Machine$double.eps)
            if (-min(D) < Dthresh) {
                indefinite <- FALSE
                D[D < Dthresh] <- Dthresh
                                   } else indefinite <- TRUE
                        } else indefinite <- FALSE
                        
        if(indefinite){
            if (eigen.fix) {
                eh <- eigen(Hp, symmetric = TRUE)
                ev <- abs(eh$values)
                Hp <- eh$vectors %*% (ev * t(eh$vectors))
                           }
                           else{
                Ib <- diag(rank) * abs(min(D))
                Ip <- diag(rank) * abs(max(D) * .Machine$double.eps^0.5)
                Hp <- Hp + Ip + Ib
                               }
            D <- rep(1, ncol(Hp))
            indefinite <- TRUE
        }
        else {
            D <- D^-0.5
            Hp <- D * t(D * Hp)
            Ip <- diag(rank) * .Machine$double.eps^0.5
        }    
       
    
         L <- suppressWarnings(chol(Hp, pivot = TRUE))
 
         mult <- 1
         while (attr(L, "rank") < rank) {
             if (eigen.fix) {
                 eh <- eigen(Hp, symmetric = TRUE)
                 ev <- eh$values
                 thresh <- max(min(ev[ev > 0]), max(ev) * 1e-06) * 
                   mult
                 mult <- mult * 10
                 ev[ev < thresh] <- thresh
                 Hp <- eh$vectors %*% (ev * t(eh$vectors))
                 L <- suppressWarnings(chol(Hp, pivot = TRUE))
             }
             else {
                 L <- suppressWarnings(chol(Hp + Ip, pivot = TRUE))
                 Ip <- Ip * 100
             }
             indefinite <- TRUE
        }

list(L = L, D = D)

}



stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  


controlEFS <- list(efs.lspmax = 15, eps = 1e-07, tol = 1e-06, tiny = .Machine$double.eps^0.5, efs.tol = 0.1)
score.hist <- rep(0, 200) 
mult       <- 1
lsp        <- log(sp)
gamma      <- 1
Mp         <- -1
eigen.fix  <- FALSE

Sl.termMult <- getFromNamespace("Sl.termMult", "mgcv")
ldetS       <- getFromNamespace("ldetS", "mgcv")

Mp <- ncol(totalPenaltySpace(qu.mag$Ss, NULL, qu.mag$off, length(fit$argument))$Z) # check this when using TE, should be ok

Sl.sfTemp <- VC$Sl.sf # check this for te and ti, should be fine

for(i in 1:length(Sl.sfTemp)) Sl.sfTemp[[i]]$D <- solve(Sl.sfTemp[[i]]$D) # PDef(Sl.sfTemp[[i]]$D)$res.inv not good here. There are i numbers in eigen vals
                                                                          # With TE or TI all ok since there is only one D

for(iter in 1:200){    
    
        o.ests <- c(fit$argument) 
        
        rp <- ldetS(VC$Sl.sf, rho = lsp, fixed = rep(FALSE, length(lsp)), np = length(fit$argument), root = TRUE) 
        o.estsStar <- Sl.initial.repara(Sl.sfTemp, o.ests, inverse = TRUE)

        Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit$hessian)$res.inv, inverse = TRUE) # this works to obtain the correctly parametrised Vb
                                                                                      # but I need L and D for REML score, so I still have to run the stuff below

        LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
        
        L <- LD$L
        D <- LD$D
      
        ipiv <- piv <- attr(L, "pivot")
        p <- length(piv)
        ipiv[piv] <- 1:p
        
        Vb <- crossprod(forwardsolve(t(L), diag(D, nrow = p)[piv, , drop = FALSE])[ipiv, , drop = FALSE]) # this matches with the stuff above
        Vb <- Sl.repara(rp$rp, Vb, inverse = TRUE)  
        
        SVb   <- Sl.termMult(VC$Sl.sf, Vb)
        trVS  <- rep(0, length(SVb))
        
        for (i in 1:length(SVb)){
            ind     <- attr(SVb[[i]], "ind")
            trVS[i] <- sum(diag(SVb[[i]][, ind]))
                                 }

        start <- Sl.repara(rp$rp, o.estsStar)
        
        Sb  <- Sl.termMult(VC$Sl.sf, start, full = TRUE)
        bSb <- rep(0, length(Sb))
        
        for(i in 1:length(Sb)) bSb[i] <- sum(start * Sb[[i]])
        
        S1 <- rp$ldet1 # rank of pens
        
        a <- pmax(controlEFS$tiny, S1 * exp(-lsp) - trVS) 
        r <- a/pmax(controlEFS$tiny, bSb)
        
        r[a == 0 & bSb == 0] <- 1
        
        r[!is.finite(r)] <- 1e+06
        
        lsp1 <- pmin(lsp + log(r) * mult, controlEFS$efs.lspmax)
        
        max.step <- max(abs(lsp1 - lsp))
        
        
          ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        old.reml <- -as.numeric((-fit$l - drop(t(fit$argument) %*% fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)   
        #old.reml <- -as.numeric((-fit$l - drop(t(start) %*% St %*% start)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)   
        #verified, difference in St, ok
        
        
        sp1 <- exp(lsp1); names(sp1) <- names(lsp1)
        if(!is.null(VC$sp.fixed)) sp1 <- VC$sp.fixed 	
        if( VC$triv == TRUE && VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp1	
        ps <- pen(qu.mag, sp1, VC, univ = respvec$univ, l.splist)       
             fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc,  
                          respvec = respvec, VC = VC, 
                          ps = ps, 
                          blather = TRUE, iterlim = iterlim), silent = sc)                                                    
                          if(inherits(fit, "try-error") || is.null(fit$l)){conv.sp <- FALSE
                                                         ps <- ps1 # pen(qu.mag, sp, VC, univ = respvec$univ, l.splist)
                                                         fit <- try( trust(func.opt, c(fit1$argument), rinit=rinit, rmax = rmax,  parscale = parsc,  
			                                              respvec = respvec, VC = VC, 
			                                              ps = ps, 
                                                                      blather = TRUE, iterlim = iterlim), silent = sc)                                                                       
                         if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
                         if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
                                                                         }  
        iter.inner <- iter.inner + fit$iterations   	                                        
        rp <- ldetS(VC$Sl.sf, rho = lsp1, fixed = rep(FALSE, length(lsp1)), np = length(fit$argument), root = TRUE)         
        Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit$hessian)$res.inv, inverse = TRUE) 
        LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D        

        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        fit$REML <- -as.numeric((-fit$l - drop(t(fit$argument) %*% fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)   




            if (fit$REML <= old.reml) {
            
                if (max.step < 0.05) {
                
                  lsp2 <- pmin(lsp + log(r) * mult * 2, 12)
 
        sp2 <- exp(lsp2); names(sp2) <- names(lsp2)
        if(!is.null(VC$sp.fixed)) sp2 <- VC$sp.fixed 	
        if( VC$triv == TRUE && VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp2	
        ps <- pen(qu.mag, sp2, VC, univ = respvec$univ, l.splist)       
             fit2 <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax, parscale = parsc,  
                          respvec = respvec, VC = VC, 
                          ps = ps, 
                          blather = TRUE, iterlim = iterlim), silent = sc)                                                    
                          if(inherits(fit2, "try-error") || is.null(fit2$l)){conv.sp <- FALSE
                                                         ps <- ps1 # pen(qu.mag, spo, VC, univ = respvec$univ, l.splist)
                                                         fit2 <- try( trust(func.opt, c(fit1$argument), rinit=rinit, rmax = rmax,  parscale = parsc,  
			                                              respvec = respvec, VC = VC, 
			                                              ps = ps, 
                                                                      blather = TRUE, iterlim = iterlim), silent = sc)                                                                       
                         if((inherits(fit2, "try-error") || is.null(fit2$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
                         if((inherits(fit2, "try-error") || is.null(fit2$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
                                                                         }  
        iter.inner <- iter.inner + fit2$iterations  
        rp <- ldetS(VC$Sl.sf, rho = lsp2, fixed = rep(FALSE, length(lsp2)), np = length(fit$argument), root = TRUE)         
        Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit2$hessian)$res.inv, inverse = TRUE) 
        LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D           
        
        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        fit2$REML <- -as.numeric((-fit2$l - drop(t(fit2$argument) %*% fit2$S.h %*% fit2$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)   

         
         
                  if(fit2$REML < fit$REML) {
                    fit <- fit2
                    lsp <- lsp2
                    mult <- mult * 2
                  }
                  else {
                    lsp <- lsp1
                  }
                }
                else lsp <- lsp1
            }
            else {
                while (fit$REML > old.reml && mult > 1) {
                  mult <- mult/2
                  lsp1 <- pmin(lsp + log(r) * mult, controlEFS$efs.lspmax)
                  
 
        sp1 <- exp(lsp1); names(sp1) <- names(lsp1)
        if(!is.null(VC$sp.fixed)) sp1 <- VC$sp.fixed 	
        if( VC$triv == TRUE && VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp1	
        ps <- pen(qu.mag, sp1, VC, univ = respvec$univ, l.splist)       
             fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc,  
                          respvec = respvec, VC = VC, 
                          ps = ps, 
                          blather = TRUE, iterlim = iterlim), silent = sc)                                                    
                          if(inherits(fit, "try-error") || is.null(fit$l)){conv.sp <- FALSE
                                                         ps <- ps1 # pen(qu.mag, spo, VC, univ = respvec$univ, l.splist)
                                                         fit <- try( trust(func.opt, c(fit1$argument), rinit=rinit, rmax = rmax,  parscale = parsc,  
			                                              respvec = respvec, VC = VC, 
			                                              ps = ps, 
                                                                      blather = TRUE, iterlim = iterlim), silent = sc)                                                                       
                         if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
                         if((inherits(fit, "try-error") || is.null(fit$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
                                                                         }  
        iter.inner <- iter.inner + fit$iterations  
        rp <- ldetS(VC$Sl.sf, rho = lsp1, fixed = rep(FALSE, length(lsp1)), np = length(fit$argument), root = TRUE)                   
        Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit$hessian)$res.inv, inverse = TRUE) 
        LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D           
        
        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        fit$REML <- -as.numeric((-fit$l - drop(t(fit$argument) %*% fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)   
              
              
                }
                lsp <- lsp1
                if (mult < 1) 
                  mult <- 1
            }
            score.hist[iter] <- fit$REML
            if (iter > 3 && max.step < 0.05 && max(abs(diff(score.hist[(iter - 3):iter]))) < controlEFS$efs.tol) break
            if (iter == 1) 
                old.ll <- fit$l
            else {
                if (abs(old.ll - fit$l) < 100 * controlEFS$eps *abs(fit$l)) 
                  break
                old.ll <- fit$l
            }



      } # end smoothing fitting loop
          
   
       sp      <- exp(lsp)
       iter.sp <- iter
       if(iter > 200 ) conv.sp <- FALSE else conv.sp <- TRUE 
       if(VC$gc.l == TRUE) gc()
       St <- crossprod(rp$E)

    
}



        


        

}else{
    
    wor.c <- working.comp(fit) 
        
    # if(is.null(VC$sp.fixed))  
    # not working properly in classical case although not needed
    
    
    bs.mgfit <- magic(wor.c$Z, wor.c$X, numeric(0), list(), numeric(0))  
    
    
    #if(!is.null(VC$sp.fixed)) bs.mgfit <- magic(y = wor.c$Z, X = wor.c$X, sp = VC$sp.fixed, S = qu.mag$Ss, off = qu.mag$off, rank = qu.mag$rank, gcv = FALSE, gamma = 1000000, control = list(tol = 1e+6)) 
    
    magpp <- magic.post.proc(wor.c$X, bs.mgfit)
        
    }




rm(fit1, ps1)

                  list(fit = fit, score.hist = score.hist,
                       iter.if = iter.if, sp.method = VC$sp.method,
                       conv.sp = conv.sp, 
                       iter.sp = iter.sp, 
                       iter.inner = iter.inner, 
                       bs.mgfit = bs.mgfit, 
                       wor.c = wor.c, 
                       sp = sp, magpp = magpp, rp = rp$rp, Sl = VC$Sl.sf, D = D, L = L, Sl.sfTemp = Sl.sfTemp, St = St)
      
 
}




