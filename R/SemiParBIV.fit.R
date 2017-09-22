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

l.splist <- list( l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                  l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                  l.sp7 = l.sp7, l.sp8 = l.sp8 )


if( ( l.sp1==0 && l.sp2==0 && l.sp3==0 && l.sp4==0 && l.sp5==0 && l.sp6==0 && l.sp7==0 && l.sp8==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0, qu.mag = NULL) else ps <- pen(qu.mag, sp, VC, univ = respvec$univ, l.splist)



if(VC$triv == TRUE){

if( VC$penCor == "ridge" ) qu.mag <- ps$qu.mag # need to do this as in the ridge case I expand S in pen
                                               # not the same for lasso and alasso as these penalties
                                               # change based on parameter estimates
                                               
# can't do lasso as it changes at each iteration  
# and we need th stuff below

if( VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp


}



###########################################################

  parsc <- rep(VC$parscale, length(start.v) ); sc <- TRUE

  fit  <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                     respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                     iterlim = iterlim), silent = sc)   

  if(class(fit) == "try-error" || is.null(fit$l)){
  
    fit  <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                       respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                       iterlim = iterlim/4), silent = sc)  
                     
        if(class(fit) == "try-error"|| is.null(fit$l)){
        
            fit  <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                               respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                               iterlim = iterlim/10), silent = sc)   


if((class(fit) == "try-error" || is.null(fit$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
 
             
                                     }
                     
  }


  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL
  
  #####################################################################

  #VC$fp <- TRUE
    if((VC$fp==FALSE && (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0)) ){

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
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                		
                	sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo)
                	
             if( VC$triv == TRUE && VC$penCor %in% c("lasso", "alasso") ) VC$sp <- sp	

             ps <- pen(qu.mag, sp, VC, univ = respvec$univ, l.splist) # in triv case, I reset ridge penalty but ok
               
             # VC$my.env$k <- VC$k, does not like it  
               
             fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc,  
                          respvec = respvec, VC = VC, 
                          ps = ps, 
                          blather = TRUE, iterlim = iterlim), silent = sc)                           
                          
                          if(class(fit) == "try-error" || is.null(fit$l)){conv.sp <- FALSE

                                                         ps <- pen(qu.mag, spo, VC, univ = respvec$univ, l.splist)

                                                         fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc,  
			                                              respvec = respvec, VC = VC, 
			                                              ps = ps, 
                                                                      blather = TRUE, iterlim = iterlim), silent = sc)  
                                                                      

if((class(fit) == "try-error" || is.null(fit$l)) && VC$gamlssfit == FALSE ) stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$gamlssfit == TRUE  ) stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")         
 
    
                                                        # break
                          } 
                          
                                                             
             iter.inner <- iter.inner + fit$iterations   	                                    
              
             if(iter.sp >= iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- abs(fit$l - fito)/(0.1 + abs(fit$l))  # max(abs(o.ests - c(fit$argument)))
                  
           
          } # end smoothing fitting loop
          

       if(VC$gc.l == TRUE) gc()   

       magpp <- magic.post.proc(wor.c$X, bs.mgfit) 
       

    }else{
    
    wor.c <- working.comp(fit) 
    
    bs.mgfit <- magic(wor.c$Z, wor.c$X, numeric(0), list(), numeric(0))    
    magpp    <- magic.post.proc(wor.c$X, bs.mgfit)
        
    }


                  list(fit = fit, 
                       iter.if = iter.if, 
                       conv.sp = conv.sp, 
                       iter.sp = iter.sp, 
                       iter.inner = iter.inner, 
                       bs.mgfit = bs.mgfit, 
                       wor.c = wor.c, 
                       sp = sp, magpp = magpp)
      
 
}




